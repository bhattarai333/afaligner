#include "sparse_matrix.h"
#include "uthash.h"
#include <string.h>
#include <stdbool.h>

SparseMatrix *create_sparse_matrix() {
    SparseMatrix *mat = malloc(sizeof(SparseMatrix));
    if (!mat) {
        return NULL;  // Memory allocation failed
    }
    mat->table = NULL;
    return mat;
}

bool set_element(SparseMatrix *mat, size_t i, size_t j, D_matrix_element value) {
    if (!mat) {
        return false;  // Invalid matrix pointer
    }

    Index key = { i, j };
    SparseMatrixEntry *entry = NULL;

    HASH_FIND(hh, mat->table, &key, sizeof(Index), entry);
    if (entry == NULL) {
        entry = malloc(sizeof(SparseMatrixEntry));
        if (!entry) {
            return false;  // Memory allocation failed
        }
        entry->key = key;
        HASH_ADD(hh, mat->table, key, sizeof(Index), entry);
        if (!mat->table) {
            free(entry);
            return false;  // Hash table insertion failed
        }
    }
    entry->value = value;
    return true;
}

D_matrix_element *get_element(SparseMatrix *mat, size_t i, size_t j) {
    if (!mat) {
        return NULL;  // Invalid matrix pointer
    }

    Index key = { i, j };
    SparseMatrixEntry *entry = NULL;

    HASH_FIND(hh, mat->table, &key, sizeof(Index), entry);
    if (entry) {
        return &entry->value;
    }
    return NULL;
}

void free_sparse_matrix(SparseMatrix *mat) {
    if (!mat) {
        return;  // Invalid matrix pointer
    }

    SparseMatrixEntry *entry, *tmp;
    HASH_ITER(hh, mat->table, entry, tmp) {
        if (entry) {
            HASH_DEL(mat->table, entry);
            free(entry);
        }
    }
    free(mat);
}