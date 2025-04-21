//helper.h
#ifndef HELPER_H
#define HELPER_H

#include <stdlib.h>

#if defined(_MSC_VER)
#include <BaseTsd.h>
#include <stdbool.h>
typedef SSIZE_T ssize_t;
#endif

// First, define the D_matrix_element structure
typedef struct {
    double distance;
    ssize_t prev_i;
    ssize_t prev_j;
} D_matrix_element;

// Then define the sparse matrix structures
typedef struct sparse_element {
    size_t i;
    size_t j;
    D_matrix_element value;
    struct sparse_element* next;
} sparse_element;

typedef struct {
    size_t rows;
    size_t cols;
    sparse_element* elements;
} sparse_matrix;

// Sparse matrix operations
sparse_matrix* create_sparse_matrix(size_t rows, size_t cols);
void set_sparse_value(sparse_matrix* matrix, size_t i, size_t j, D_matrix_element value);
D_matrix_element get_sparse_value(sparse_matrix* matrix, size_t i, size_t j);
void free_sparse_matrix(sparse_matrix* matrix);

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
    double *s,
    double *t,
    size_t n,
    size_t m,
    size_t l,
    double skip_penalty,
    size_t *window,
    double *path_distance,
    size_t *path_buffer
);

ssize_t FastDTWBD(
    double *s,
    double *t,
    size_t n,
    size_t m,
    size_t l,
    double skip_penalty,
    int radius,
    double *path_distance,
    size_t *path_buffer
);

#endif // HELPER_H