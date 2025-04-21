#include "sparse_matrix.h"
#include "uthash.h"
#include <string.h>

SparseMatrix *create_sparse_matrix() {
    SparseMatrix *mat = malloc(sizeof(SparseMatrix));
    mat->table = NULL;
    return mat;
}

void set_element(SparseMatrix *mat, size_t i, size_t j, D_matrix_element value) {
    Index key = { i, j };
    SparseMatrixEntry *entry = NULL;

    HASH_FIND(hh, mat->table, &key, sizeof(Index), entry);
    if (entry == NULL) {
        entry = malloc(sizeof(SparseMatrixEntry));
        entry->key = key;
        HASH_ADD(hh, mat->table, key, sizeof(Index), entry);
    }
    entry->value = value;
}

D_matrix_element *get_element(SparseMatrix *mat, size_t i, size_t j) {
    Index key = { i, j };
    SparseMatrixEntry *entry = NULL;

    HASH_FIND(hh, mat->table, &key, sizeof(Index), entry);
    if (entry) {
        return &entry->value;
    }
    return NULL;
}

void free_sparse_matrix(SparseMatrix *mat) {
    SparseMatrixEntry *entry, *tmp;
    HASH_ITER(hh, mat->table, entry, tmp) {
        HASH_DEL(mat->table, entry);
        free(entry);
    }
    free(mat);
}
