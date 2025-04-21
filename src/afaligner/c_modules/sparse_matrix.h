// sparse_matrix.h
#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <stdlib.h>
#include <stddef.h>

// Matrix element structure
typedef struct {
    double distance;
    ssize_t prev_i;
    ssize_t prev_j;
} D_matrix_element;

// Sparse matrix element structure
typedef struct sparse_element {
    size_t i;
    size_t j;
    D_matrix_element value;
    struct sparse_element* next;
} sparse_element;

// Sparse matrix structure
typedef struct {
    size_t rows;
    size_t cols;
    sparse_element* elements;
} sparse_matrix;

// Function declarations
sparse_matrix* create_sparse_matrix(size_t rows, size_t cols);
void set_value(sparse_matrix* matrix, size_t i, size_t j, D_matrix_element value);
D_matrix_element get_value(sparse_matrix* matrix, size_t i, size_t j);
void free_sparse_matrix(sparse_matrix* matrix);

#endif // SPARSE_MATRIX_H