#ifndef FASTDTWBD_H
#define FASTDTWBD_H

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include "dtwbd.h"  // Including dtwbd.h for shared structures and functions

#if defined(_WIN32) || defined(__WIN32__)
    #ifdef BUILDING_FASTDTWBD
        #define EXPORT __declspec(dllexport)
    #else
        #define EXPORT __declspec(dllimport)
    #endif
#else
    #ifdef BUILDING_FASTDTWBD
        #define EXPORT __attribute__((visibility("default")))
    #else
        #define EXPORT
    #endif
#endif

// FastDTWBD function prototype (make sure it's marked with EXPORT)
EXPORT ssize_t FastDTWBD(
    double *s,  // first sequence of MFCC frames – n x l contiguous array
    double *t,  // second sequence of MFCC frames – m x l contiguous array
    size_t n,   // number of frames in first sequence
    size_t m,   // number of frames in second sequence
    size_t l,   // number of MFCCs per frame
    double skip_penalty,    // penalty for skipping one frame
    int radius,             // radius of path projection
    double *path_distance,  // place to store warping path distance
    size_t *path_buffer     // buffer to store resulting warping path – (n+m) x 2 contiguous array
);

// Additional helper function prototypes if needed for FastDTWBD implementation
EXPORT double *get_coarsed_sequence(double *s, size_t n, size_t l);
EXPORT size_t *get_window(size_t n, size_t m, size_t *path_buffer, size_t path_len, int radius);
EXPORT void update_window(size_t *window, size_t n, size_t m, ssize_t i, ssize_t j);

#endif // FASTDTWBD_H
