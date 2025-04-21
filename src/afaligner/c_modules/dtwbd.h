#ifndef DTWBD_H
#define DTWBD_H

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include "uthash.h"

// Matrix element used for DTW/BD computation
typedef struct {
    double distance;
    ssize_t prev_i;
    ssize_t prev_j;
} D_matrix_element;

// Main DTWBD function
ssize_t DTWBD(
    double *x, size_t n,
    double *y, size_t m,
    size_t dim,
    double skip_penalty,
    size_t *window,
    size_t *path_buffer,
    double *path_distance
);

double *get_coarsed_sequence(double *s, size_t n, size_t l);

size_t *get_window(size_t n, size_t m, size_t *path_buffer, size_t path_len, int radius);

void update_window(size_t *window, size_t n, size_t m, ssize_t i, ssize_t j);

// Helper functions
double euclid_distance(double *x, double *y, size_t l);

double get_distance(
    void *matrix, size_t n, size_t m,
    size_t *window, size_t i, size_t j,
    int use_sparse
);

D_matrix_element get_best_candidate(D_matrix_element *candidates, size_t n);

void reverse_path(size_t *path, ssize_t path_len);

#endif // DTWBD_H
