#include "sparse_matrix.h"
#include "uthash.h"
#include <string.h>
#include <stdbool.h>
#include "logger.h"
#include <stdlib.h>

SparseMatrix *create_sparse_matrix() {
    SparseMatrix *mat = malloc(sizeof(SparseMatrix));
    if (!mat) {
        log_info("[SparseMatrix] Failed to allocate memory for SparseMatrix.");
        return NULL;
    }

    LOG_ALLOC(mat, sizeof(SparseMatrix));
    log_debug("[SparseMatrix] SparseMatrix created at %p", (void *)mat);

    mat->table = NULL;
    return mat;
}

bool set_element(SparseMatrix *mat, size_t i, size_t j, D_matrix_element value) {
    if (!mat) {
        log_warn("[set_element] set_element called with NULL matrix pointer.");
        return false;
    }

    Index key = { i, j };
    SparseMatrixEntry *entry = NULL;

    // First try to find existing entry
    HASH_FIND(hh, mat->table, &key, sizeof(Index), entry);
    if (entry) {
        // Update existing entry
        log_debug("[set_element] Updating existing entry at (%zu, %zu) from distance = %f to distance = %f",
                  i, j, entry->value.distance, value.distance);
        entry->value = value;
        return true;
    }

    // Only allocate if entry doesn't exist
    entry = malloc(sizeof(SparseMatrixEntry));
    if (!entry) {
        log_error("[set_element] Failed to allocate memory for SparseMatrixEntry at (%zu, %zu)", i, j);
        return false;
    }

    LOG_ALLOC(entry, sizeof(SparseMatrixEntry));
    entry->key = key;
    entry->value = value;

    // Add to hash table
    HASH_ADD(hh, mat->table, key, sizeof(Index), entry);
    if (!mat->table) {
        LOG_FREE(entry);
        free(entry);
        log_error("[set_element] Failed to add entry to hash table at (%zu, %zu)", i, j);
        return false;
    }

    log_debug("[set_element] Created new entry at (%zu, %zu) with distance = %f",
              i, j, value.distance);
    return true;
}

D_matrix_element *get_element(SparseMatrix *mat, size_t i, size_t j) {
    if (!mat) {
        log_warn("[get_element] get_element called with NULL matrix pointer.");
        return NULL;
    }

    Index key = { i, j };
    SparseMatrixEntry *entry = NULL;

    HASH_FIND(hh, mat->table, &key, sizeof(Index), entry);
    if (entry) {
        // Log the access with details of the element
        log_debug("[get_element] Accessed element at (%zu, %zu) with distance = %f, prev_i = %zd, prev_j = %zd",
                  i, j, entry->value.distance, entry->value.prev_i, entry->value.prev_j);
        return &entry->value;
    } else {
        log_debug("[get_element] No entry found at (%zu, %zu)", i, j);
    }

    return NULL;
}

void free_sparse_matrix(SparseMatrix *mat) {
    if (!mat) {
        log_warn("[free_sparse_matrix] free_sparse_matrix called with NULL matrix pointer.");
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
    log_info("[free_sparse_matrix] Freed sparse matrix and %zu entries.", count);
}