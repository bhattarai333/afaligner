#include "sparse_matrix.h"
#include "uthash.h"
#include <string.h>
#include <stdbool.h>
#include "logger.h"
#include <stdlib.h>

SparseMatrix *create_sparse_matrix() {
    SparseMatrix *mat = malloc(sizeof(SparseMatrix));
    if (!mat) {
        log_info("Failed to allocate memory for SparseMatrix.");
        return NULL;
    }

    LOG_ALLOC(mat, sizeof(SparseMatrix));
    log_debug("SparseMatrix created at %p", (void *)mat);

    mat->table = NULL;
    return mat;
}

bool set_element(SparseMatrix *mat, size_t i, size_t j, D_matrix_element value) {
    if (!mat) {
        log_warn("set_element called with NULL matrix pointer.");
        return false;
    }

    Index key = { i, j };
    SparseMatrixEntry *entry = NULL;

    HASH_FIND(hh, mat->table, &key, sizeof(Index), entry);
    if (entry == NULL) {
        entry = malloc(sizeof(SparseMatrixEntry));
        if (!entry) {
            log_error("Failed to allocate memory for SparseMatrixEntry at (%zu, %zu)", i, j);
            return false;
        }

        LOG_ALLOC(entry, sizeof(SparseMatrixEntry));
        entry->key = key;
        HASH_ADD(hh, mat->table, key, sizeof(Index), entry);

        if (!mat->table) {
            LOG_FREE(entry);
            free(entry);
            log_error("Failed to add entry to hash table at (%zu, %zu)", i, j);
            return false;
        }

        log_debug("Created new entry at (%zu, %zu) with value = %f", i, j, (double)value);
    } else {
        log_debug("Updated existing entry at (%zu, %zu) from %f to %f", i, j, (double)entry->value, (double)value);
    }

    entry->value = value;
    return true;
}

D_matrix_element *get_element(SparseMatrix *mat, size_t i, size_t j) {
    if (!mat) {
        log_warn("get_element called with NULL matrix pointer.");
        return NULL;
    }

    Index key = { i, j };
    SparseMatrixEntry *entry = NULL;

    HASH_FIND(hh, mat->table, &key, sizeof(Index), entry);
    if (entry) {
        log_debug("Accessed element at (%zu, %zu) with value = %f", i, j, (double)entry->value);
        return &entry->value;
    } else {
        log_debug("No entry found at (%zu, %zu)", i, j);
    }

    return NULL;
}

void free_sparse_matrix(SparseMatrix *mat) {
    if (!mat) {
        log_warn("free_sparse_matrix called with NULL matrix pointer.");
        return;
    }

    size_t count = 0;
    SparseMatrixEntry *entry, *tmp;
    HASH_ITER(hh, mat->table, entry, tmp) {
        if (entry) {
            HASH_DEL(mat->table, entry);
            free(entry);
            count++;
        }
    }

    free(mat);
    log_info("Freed sparse matrix and %zu entries.", count);
}