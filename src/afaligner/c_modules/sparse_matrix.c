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

        // Log with details from D_matrix_element (distance, prev_i, prev_j)
        log_debug("Created new entry at (%zu, %zu) with distance = %f, prev_i = %zd, prev_j = %zd",
                  i, j, value.distance, value.prev_i, value.prev_j);
    } else {
        // Log with previous entry's details and updated value
        log_debug("Updated existing entry at (%zu, %zu) from distance = %f, prev_i = %zd, prev_j = %zd to distance = %f, prev_i = %zd, prev_j = %zd",
                  i, j, entry->value.distance, entry->value.prev_i, entry->value.prev_j,
                  value.distance, value.prev_i, value.prev_j);
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
        // Log the access with details of the element
        log_debug("Accessed element at (%zu, %zu) with distance = %f, prev_i = %zd, prev_j = %zd",
                  i, j, entry->value.distance, entry->value.prev_i, entry->value.prev_j);
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