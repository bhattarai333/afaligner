#ifndef HELPER_H
#define HELPER_H

#include <stddef.h>
#include <sys/types.h>

typedef struct D_matrix_element {
    double distance;
    ssize_t prev_i;
    ssize_t prev_j;
} D_matrix_element;

typedef struct sparse_element {
    size_t i;
    size_t j;
    D_matrix_element value;
    struct sparse_element* next;
} sparse_element;

typedef struct sparse_matrix {
    size_t rows;
    size_t cols;
    sparse_element* elements;
} sparse_matrix;

// Sparse matrix operations
sparse_matrix* create_sparse_matrix(size_t rows, size_t cols);
void set_sparse_value(sparse_matrix* matrix, size_t i, size_t j, D_matrix_element value);
D_matrix_element get_sparse_value(sparse_matrix* matrix, size_t i, size_t j);
void free_sparse_matrix(sparse_matrix* matrix);

// Helper functions
double* get_coarsed_sequence(double* s, size_t n, size_t l);
size_t* get_window(size_t n, size_t m, size_t* path_buffer, size_t path_len, int radius);
void update_window(size_t* window, size_t n, size_t m, ssize_t i, ssize_t j);
double euclid_distance(double* x, double* y, size_t l);
D_matrix_element get_best_candidate(D_matrix_element* candidates, size_t n);
void reverse_path(size_t* path, ssize_t path_len);

#endif // HELPER_H