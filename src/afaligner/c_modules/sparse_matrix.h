#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <stdlib.h>
#include "dtwbd.h" // for D_matrix_element
#include "uthash.h"


typedef struct {
    size_t i;
    size_t j;
} Index;

typedef struct SparseMatrixEntry {
    Index key;
    D_matrix_element value;
    struct SparseMatrixEntry *hh_next;
    UT_hash_handle hh; // required by uthash
} SparseMatrixEntry;

typedef struct {
    SparseMatrixEntry *table;
} SparseMatrix;

// API
SparseMatrix *create_sparse_matrix();
void set_element(SparseMatrix *mat, size_t i, size_t j, D_matrix_element value);
D_matrix_element *get_element(SparseMatrix *mat, size_t i, size_t j);
void free_sparse_matrix(SparseMatrix *mat);

#endif
