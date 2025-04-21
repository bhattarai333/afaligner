//sparse_matrix.c
#include "sparse_matrix.h"
#include "logger.h"
#include <float.h>

sparse_matrix* create_sparse_matrix(size_t rows, size_t cols) {
    sparse_matrix* matrix = malloc(sizeof(sparse_matrix));
    if (!matrix) {
        log_error("Failed to allocate sparse matrix");
        return NULL;
    }
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->elements = NULL;
    return matrix;
}

void set_value(sparse_matrix* matrix, size_t i, size_t j, D_matrix_element value) {
    sparse_element* elem = malloc(sizeof(sparse_element));
    if (!elem) {
        log_error("Failed to allocate sparse matrix element at (%zu, %zu)", i, j);
        return;
    }
    elem->i = i;
    elem->j = j;
    elem->value = value;
    elem->next = matrix->elements;
    matrix->elements = elem;
}

D_matrix_element get_value(sparse_matrix* matrix, size_t i, size_t j) {
    sparse_element* current = matrix->elements;
    while (current) {
        if (current->i == i && current->j == j) {
            return current->value;
        }
        current = current->next;
    }
    D_matrix_element default_elem = {DBL_MAX, -1, -1};
    return default_elem;
}

void free_sparse_matrix(sparse_matrix* matrix) {
    sparse_element* current = matrix->elements;
    while (current) {
        sparse_element* next = current->next;
        free(current);
        current = next;
    }
    free(matrix);
}