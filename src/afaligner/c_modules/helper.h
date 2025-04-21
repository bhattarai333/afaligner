#ifndef HELPER_H
#define HELPER_H

#include <stdlib.h>

#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

// Struct definition for D_matrix element
typedef struct {
    double distance;
    ssize_t prev_i;
    ssize_t prev_j;
} D_matrix_element;

// Function declarations for sequence processing
double *get_coarsed_sequence(double *s, size_t n, size_t l);

// Window management functions
size_t *get_window(size_t n, size_t m, size_t *path_buffer, size_t path_len, int radius);
void update_window(size_t *window, size_t n, size_t m, ssize_t i, ssize_t j);

// Distance calculation functions
double euclid_distance(double *x, double *y, size_t l);
double get_distance(D_matrix_element *D_matrix, size_t n, size_t m, size_t *window, size_t i, size_t j);

// Path processing functions
D_matrix_element get_best_candidate(D_matrix_element *candidates, size_t n);
void reverse_path(size_t *path, ssize_t path_len);

// Main algorithm functions
ssize_t DTWBD(
    double *s,          // first sequence of MFCC frames – n x l contiguous array
    double *t,          // second sequence of MFCC frames – m x l contiguous array
    size_t n,           // number of frames in first sequence
    size_t m,           // number of frames in second sequence
    size_t l,           // number of MFCCs per frame
    double skip_penalty,// penalty for skipping one frame
    size_t *window,     // n x 2 contiguous array for frame ranges
    double *path_distance,  // place to store warping path distance
    size_t *path_buffer    // buffer to store resulting warping path
);

ssize_t FastDTWBD(
    double *s,          // first sequence of MFCC frames – n x l contiguous array
    double *t,          // second sequence of MFCC frames – m x l contiguous array
    size_t n,           // number of frames in first sequence
    size_t m,           // number of frames in second sequence
    size_t l,           // number of MFCCs per frame
    double skip_penalty,// penalty for skipping one frame
    int radius,         // radius of path projection
    double *path_distance,  // place to store warping path distance
    size_t *path_buffer    // buffer to store resulting warping path
);

#endif // HELPER_H