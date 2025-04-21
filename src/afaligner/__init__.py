from datetime import timedelta
import json
import math
import os.path
import shutil
import subprocess

from aeneas.audiofilemfcc import AudioFileMFCC
from aeneas.language import Language
from aeneas.synthesizer import Synthesizer
from aeneas.textfile import TextFile, TextFileFormat
from aeneas.exacttiming import TimeValue
import numpy as np
import jinja2

from afaligner.c_dtwbd_wrapper import c_FastDTWBD

BASE_DIR = os.path.dirname(os.path.realpath(__file__))


def align(
        text_dir, audio_dir, output_dir=None, output_format='smil',
        sync_map_text_path_prefix='', sync_map_audio_path_prefix='',
        skip_penalty=None, radius=None,
        times_as_timedelta=False, language=Language.ENG,
):

    if skip_penalty is None:
        skip_penalty = 0.75

    if radius is None:
        radius = 100

    if output_dir is not None:
        tmp_dir = os.path.join(output_dir, 'tmp')
    else:
        parent_dir = os.path.dirname(os.path.dirname(os.path.join(text_dir, '')))
        tmp_dir = os.path.join(parent_dir, 'tmp')

    os.makedirs(tmp_dir, exist_ok=True)

    text_paths = (os.path.join(text_dir, f) for f in sorted(os.listdir(text_dir)))
    audio_paths = (os.path.join(audio_dir, f) for f in sorted(os.listdir(audio_dir)))

    sync_map = build_sync_map(
        text_paths, audio_paths, tmp_dir,
        sync_map_text_path_prefix=sync_map_text_path_prefix,
        sync_map_audio_path_prefix=sync_map_audio_path_prefix,
        skip_penalty=skip_penalty,
        radius=radius,
        times_as_timedelta=times_as_timedelta,
        language=language,
    )

    if output_dir is not None:
        if output_format == 'smil':
            output_smil(sync_map, output_dir)
        elif output_format == 'json':
            output_json(sync_map, output_dir)

    shutil.rmtree(tmp_dir)

    return sync_map


def build_sync_map(
        text_paths, audio_paths, tmp_dir,
        sync_map_text_path_prefix, sync_map_audio_path_prefix,
        skip_penalty, radius,
        times_as_timedelta,
        language,
):


    synthesizer = Synthesizer()
    parse_parameters = {'is_text_unparsed_id_regex': 'f[0-9]+'}

    sync_map = {}
    process_next_text = True
    process_next_audio = True

    while True:
        if process_next_text:
            try:
                text_path = next(text_paths)
            except StopIteration:
                break

            text_name = get_name_from_path(text_path)
            output_text_name = os.path.join(sync_map_text_path_prefix, text_name)
            textfile = TextFile(text_path, file_format=TextFileFormat.UNPARSED, parameters=parse_parameters)
            textfile.set_language(language)
            text_wav_path = os.path.join(tmp_dir, f'{drop_extension(text_name)}_text.wav')
            sync_map[output_text_name] = {}

            # Produce synthesized audio, get anchors
            anchors, _, _ = synthesizer.synthesize(textfile, text_wav_path)

            # Get fragments, convert anchors timings to the frames indicies
            fragments = [a[1] for a in anchors]
            anchors = np.array([int(a[0] / TimeValue('0.040')) for a in anchors])

            # MFCC frames sequence memory layout is a n x l 2D array,
            # where n - number of frames and l - number of MFFCs
            # i.e it is c-contiguous, but after dropping the first coefficient it siezes to be c-contiguous.
            # Should decide whether to make a copy or to work around the first coefficient.
            text_mfcc_sequence = np.ascontiguousarray(
                AudioFileMFCC(text_wav_path).all_mfcc.T[:, 1:]
            )

        if process_next_audio:
            try:
                audio_path = next(audio_paths)
            except StopIteration:
                break

            audio_name = get_name_from_path(audio_path)
            output_audio_name = os.path.join(sync_map_audio_path_prefix, audio_name)
            audio_wav_path = os.path.join(tmp_dir, f'{drop_extension(audio_name)}_audio.wav')
            subprocess.run(['ffmpeg', '-n', '-i', audio_path, '-rf64', 'auto', audio_wav_path])

            audio_mfcc_sequence = np.ascontiguousarray(
                AudioFileMFCC(audio_wav_path).all_mfcc.T[:, 1:]
            )

            # Keep track to calculate frames timings
            audio_start_frame = 0

        n = len(text_mfcc_sequence)
        m = len(audio_mfcc_sequence)

        _, path = c_FastDTWBD(text_mfcc_sequence, audio_mfcc_sequence, skip_penalty, radius=radius)

        if len(path) == 0:
            print(
                f'No match between {text_name} and {audio_name}. '
                f'Alignment is terminated. '
                f'Adjust skip_penalty or input files.'
            )
            return {}

        # Project path to the text and audio sequences
        text_path_frames = path[:, 0]
        audio_path_frames = path[:, 1]

        last_matched_audio_frame = audio_path_frames[-1]

        # Find first and last matched frames
        first_matched_text_frame = text_path_frames[0]
        last_matched_text_frame = text_path_frames[-1]

        # Map only those fragments that intersect matched frames
        anchors_boundary_indices = np.searchsorted(
            anchors, [first_matched_text_frame, last_matched_text_frame]
        )
        map_anchors_from = max(anchors_boundary_indices[0] - 1, 0)
        map_anchors_to = anchors_boundary_indices[1]
        anchors_to_map = anchors[map_anchors_from:map_anchors_to]
        fragments_to_map = fragments[map_anchors_from:map_anchors_to]

        # Get anchors indicies in the path projection to the text sequence
        text_path_anchor_indices = np.searchsorted(text_path_frames, anchors_to_map)

        # Get anchors' frames in audio sequence, calculate their timings
        anchors_matched_frames = audio_path_frames[text_path_anchor_indices]
        timings = (np.append(anchors_matched_frames, audio_path_frames[-1]) + audio_start_frame) * 0.040

        # Map fragment_ids to timings, update mapping of the current text file
        fragment_map = {
            f: {
                'audio_file': output_audio_name,
                'begin_time': format_time(bt, times_as_timedelta),
                'end_time': format_time(et, times_as_timedelta),
            }
            for f, bt, et in zip(fragments_to_map, timings[:-1], timings[1:])
        }

        sync_map[output_text_name].update(fragment_map)

        # Decide whether to process next file or to align the tail of the current one

        if map_anchors_to == len(anchors):
            # Process next text if no fragments are left
            process_next_text = True
        else:
            # Otherwise align tail of the current text
            process_next_text = False
            text_mfcc_sequence = text_mfcc_sequence[last_matched_text_frame:]
            fragments = fragments[map_anchors_to:]
            anchors = anchors[map_anchors_to:] - last_matched_text_frame

        if last_matched_audio_frame == m - 1 or not process_next_text:
            # Process next audio if there are no unmatched audio frames in the tail
            # or there are more text fragments to map, i.e.
            # we choose to process next audio if we cannot decide.
            # This strategy is correct if there are no extra fragments in the end.
            process_next_audio = True
        else:
            # Otherwise align tail of the current audio
            process_next_audio = False
            audio_mfcc_sequence = audio_mfcc_sequence[last_matched_audio_frame:]
            audio_start_frame += last_matched_audio_frame

    return sync_map


def get_name_from_path(path):
    return os.path.split(path)[1]


def drop_extension(path):
    return os.path.splitext(path)[0]


def format_time(t, as_timedelta=False):
    tdelta = timedelta(seconds=t)
    if as_timedelta:
        return tdelta
    return timedelta_to_str(tdelta)


def timedelta_to_str(tdelta):
    hours = int(tdelta.total_seconds()) // 3600
    minutes = int(tdelta.total_seconds() % 3600) // 60
    seconds = int(tdelta.total_seconds()) % 60
    ms = int(tdelta.microseconds) // 1000
    return f'{hours:d}:{minutes:0>2d}:{seconds:0>2d}.{ms:0>3d}'


def output_smil(sync_map, output_dir):
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(os.path.join(BASE_DIR, 'templates/')),
        autoescape=True
    )
    template = env.get_template('template.smil')

    for text_path, fragments in sync_map.items():
        parallels = []
        n = get_number_of_digits_to_name(len(fragments))
        for i, t in enumerate(fragments.items(), start=1):
            fragment_id, info = t
            # EPUB3 standard requires clipBegin < clipEnd
            if info['begin_time'] != info['end_time']:
                parallels.append({
                    'id': f'par{i:0>{n}}',
                    'fragment_id': fragment_id,
                    'audio_path': info['audio_file'],
                    'begin_time': info['begin_time'],
                    'end_time': info['end_time'],
                })

        smil = template.render(sequentials=[{
            'id': 'seq1',
            'text_path': text_path,
            'parallels': parallels
        }])

        text_name = get_name_from_path(text_path)
        file_path = os.path.join(output_dir, f'{drop_extension(text_name)}.smil')
        with open(file_path, 'w') as f:
            f.write(smil)


def get_number_of_digits_to_name(num):
    if num <= 0:
        return 0

    return math.floor(math.log10(num)) + 1


def output_json(sync_map, output_dir):
    for text_path, fragments in sync_map.items():
        text_name = get_name_from_path(text_path)
        file_path = os.path.join(output_dir, f'{drop_extension(text_name)}.json')
        with open(file_path, 'w') as f:
            json.dump(fragments, f, indent=2)


def print_sync_map(sync_map):
    for text, fragment_map in sync_map.items():
        print(text)
        for fragment, info in fragment_map.items():
            print(fragment, info['audio_file'], info['begin_time'], info['end_time'])