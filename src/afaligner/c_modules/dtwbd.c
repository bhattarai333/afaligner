#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>  // For logging
#include "dtwbd.h"
#include "sparse_matrix.h"
#include "uthash.h"
#include <stdbool.h>
#include "fastdtwbd.h"


#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

// Log file pointer
FILE *log_file = NULL;

// Open log file for writing
void open_log_file() {
    log_file = fopen("./output/afaligner.log", "a");
    if (!log_file) {
        perror("Error opening log file");
        exit(1);
    }
}

// Close log file
void close_log_file() {
    if (log_file) {
        fclose(log_file);
    }
}

// Log a message to the log file
void log_message(const char *message) {
    if (log_file) {
        fprintf(log_file, "%s\n", message);
        fflush(log_file);  // Ensure it's written immediately
    }
}

ssize_t DTWBD(
    double *s, size_t n,
    double *t, size_t m,
    size_t dim,
    double skip_penalty,
    size_t *window,
    size_t *path_buffer,
    double *path_distance
) {
    log_message("Entering DTWBD function");

    SparseMatrix *D_sparse = create_sparse_matrix();
    if (!D_sparse) {
        log_message("Error: Failed to create sparse matrix");
        return -1;
    }

    double d;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            if (window && (j < window[2*i] || j >= window[2*i+1])) {
                continue;
            }

            d = euclid_distance(&s[i*dim], &t[j*dim], dim);

            double min_prev_distance = DBL_MAX;
            ssize_t prev_i = -1, prev_j = -1;

            D_matrix_element *e;
            e = get_element(D_sparse, i - 1, j);
            if (e && e->distance < min_prev_distance) {
                min_prev_distance = e->distance;
                prev_i = i - 1;
                prev_j = j;
            }

            e = get_element(D_sparse, i, j - 1);
            if (e && e->distance < min_prev_distance) {
                min_prev_distance = e->distance;
                prev_i = i;
                prev_j = j - 1;
            }

            e = get_element(D_sparse, i - 1, j - 1);
            if (e && e->distance < min_prev_distance) {
                min_prev_distance = e->distance;
                prev_i = i - 1;
                prev_j = j - 1;
            }

            D_matrix_element cur;
            cur.distance = d + (min_prev_distance == DBL_MAX ? 0 : min_prev_distance);
            cur.prev_i = prev_i;
            cur.prev_j = prev_j;

            set_element(D_sparse, i, j, cur);
        }
    }

    double min_path_distance = DBL_MAX;
    size_t end_i = 0, end_j = 0;
    bool match = false;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            D_matrix_element *e = get_element(D_sparse, i, j);
            if (!e) continue;

            double cur_path_distance = e->distance + skip_penalty * (n - i + m - j - 2);

            if (cur_path_distance < min_path_distance) {
                min_path_distance = cur_path_distance;
                end_i = i;
                end_j = j;
                match = true;
            }
        }
    }

    ssize_t path_len = 0;

    if (match) {
        D_matrix_element *e;
        *path_distance = min_path_distance;
        for (ssize_t i = end_i, j = end_j; i != -1; ) {
            e = get_element(D_sparse, i, j);
            if (!e) break;
            path_buffer[2 * path_len] = i;
            path_buffer[2 * path_len + 1] = j;
            path_len++;
            ssize_t next_i = e->prev_i;
            ssize_t next_j = e->prev_j;
            i = next_i;
            j = next_j;
        }
        reverse_path(path_buffer, path_len);
    }

    free_sparse_matrix(D_sparse);

    log_message("Exiting DTWBD function");
    return path_len;
}

double euclid_distance(double *x, double *y, size_t l) {
    log_message("Entering euclid_distance function");

    double sum = 0;
    for (size_t i = 0; i < l; i++) {
        double v = x[i] - y[i];
        sum += v * v;
    }

    log_message("Exiting euclid_distance function");
    return sqrt(sum);
}

void reverse_path(size_t *path, ssize_t path_len) {
    log_message("Entering reverse_path function");

    for (size_t i = 0, j = path_len - 1; i < j; i++, j--) {
        size_t tmp_s = path[2 * i];
        size_t tmp_t = path[2 * i + 1];
        path[2 * i] = path[2 * j];
        path[2 * i + 1] = path[2 * j + 1];
        path[2 * j] = tmp_s;
        path[2 * j + 1] = tmp_t;
    }

    log_message("Exiting reverse_path function");
}

ssize_t FastDTWBD(
    double *s,  // first sequence of MFCC frames – n x l contiguous array
    double *t,  // second sequence of MFCC frames – m x l contiguous array
    size_t n,   // number of frames in first sequence
    size_t  m,   // number of frames in second sequence
    size_t  l,   // number of MFCCs per frame
    double skip_penalty,    // penalty for skipping one frame
    int radius,             // radius of path projection
    double *path_distance,  // place to store warping path distance
    size_t *path_buffer     // buffer to store resulting warping path – (n+m) x 2 contiguous array
) {
    log_message("Entering FastDTWBD function");

    ssize_t path_len;
    size_t min_sequence_len = 2 * (radius + 1) + 1;

    if (n < min_sequence_len || m < min_sequence_len) {
        log_message("Sequence length too short, using DTWBD directly");
        return DTWBD(s, t, n, m, l, skip_penalty, NULL, path_distance, path_buffer);
    }

    double *coarsed_s = get_coarsed_sequence(s, n, l);
    double *coarsed_t = get_coarsed_sequence(t, m, l);

    path_len = FastDTWBD(coarsed_s, coarsed_t, n/2, m/2, l, skip_penalty, radius, path_distance, path_buffer);

    size_t *window = get_window(n, m, path_buffer, path_len, radius);

    path_len = DTWBD(s, t, n, m, l, skip_penalty, window, path_distance, path_buffer);

    free(coarsed_s);
    free(coarsed_t);
    free(window);

    log_message("Exiting FastDTWBD function");
    return path_len;
}

double *get_coarsed_sequence(double *s, size_t n, size_t l) {
    log_message("Entering get_coarsed_sequence: n = %zu, l = %zu\n", n, l);

    if (n % 2 != 0) {
        log_message("Warning: n is not even, which may cause issues in coarsing.\n");
    }

    size_t coarsed_sequence_len = n / 2;
    log_message("Coarsed sequence length: %zu\n", coarsed_sequence_len);

    double *coarsed_sequence = malloc(coarsed_sequence_len * l * sizeof(double));
    if (!coarsed_sequence) {
        log_message("Error: Memory allocation for coarsed_sequence failed.\n");
        return NULL;
    }

    log_message("Memory allocated for coarsed_sequence.\n");

    for (size_t i = 0; 2 * i + 1 < n; i++) {
        for (size_t j = 0; j < l; j++) {
            coarsed_sequence[l * i + j] = (s[l * (2 * i) + j] + s[l * (2 * i + 1) + j]) / 2;
        }
    }

    log_message("Exiting get_coarsed_sequence.\n");
    return coarsed_sequence;
}


size_t *get_window(size_t n, size_t m, size_t *path_buffer, size_t path_len, int radius) {
    log_message("Entering get_window: n = %zu, m = %zu, path_len = %zu, radius = %d\n", n, m, path_len, radius);

    size_t *window = malloc(2 * n * sizeof(size_t));
    if (!window) {
        log_message("Error: Memory allocation for window failed.\n");
        return NULL;
    }
    log_message("Memory allocated for window.\n");

    for (size_t i = 0; i < n; i++) {
        window[2 * i] = m;    // maximum value for lower limit
        window[2 * i + 1] = 0;  // minimum value for upper limit
    }

    log_message("Initial window setup complete. Starting to update window with path buffer.\n");

    for (size_t k = 0; k < path_len; k++) {
        size_t i = path_buffer[2 * k];
        size_t j = path_buffer[2 * k + 1];

        log_message("Processing path_buffer[%zu]: i = %zu, j = %zu\n", k, i, j);

        for (ssize_t x = -radius; x < radius + 1; x++) {
            log_message("Updating window limits for x = %zd\n", x);

            // update lower window limit
            update_window(window, n, m, 2 * (i + x), 2 * (j - radius));
            update_window(window, n, m, 2 * (i + x) + 1, 2 * (j - radius));

            // update upper window limit
            update_window(window, n, m, 2 * (i + x), 2 * (j + radius + 1) + 1);
            update_window(window, n, m, 2 * (i + x) + 1, 2 * (j + radius + 1) + 1);
        }
    }

    log_message("Exiting get_window.\n");
    return window;
}


void update_window(size_t *window, size_t n, size_t m, ssize_t i, ssize_t j) {
    if (i < 0 || i >= n) {
        log_message("update_window: Index i = %zd is out of bounds (0 to %zu).\n", i, n - 1);
        return;
    }

    if (j < 0) {
        j = 0;
    }
    if (j > m - 1) {
        j = m - 1;
    }

    log_message("update_window: i = %zd, j = %zd, current window limits: [%zu, %zu]\n", i, j, window[2 * i], window[2 * i + 1]);

    if (j < window[2 * i]) {
        window[2 * i] = j;
        log_message("update_window: Lower limit updated for i = %zd: new window[2 * i] = %zu\n", i, j);
    }

    if (j >= window[2 * i + 1]) {
        window[2 * i + 1] = j + 1;
        log_message("update_window: Upper limit updated for i = %zd: new window[2 * i + 1] = %zu\n", i, j + 1);
    }
}
