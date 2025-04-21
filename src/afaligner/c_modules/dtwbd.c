#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "dtwbd.h"
#include "sparse_matrix.h"
#include "uthash.h"
#include <stdbool.h>


ssize_t dtwbd(
    double *s, size_t n,
    double *t, size_t m,
    size_t dim,
    double skip_penalty,
    size_t *window,
    size_t *path_buffer,
    double *path_distance
) {
    SparseMatrix *D_sparse = create_sparse_matrix();
    double d;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            if (window && (j < window[2*i] || j >= window[2*i+1])) {
                continue;
            }

            d = euclid_distance(&s[i*dim], &t[j*dim], dim);

            double min_prev_distance = DBL_MAX;
            ssize_t prev_i = -1, prev_j = -1;

            D_matrix_element *e;
            e = get_element(D_sparse, i - 1, j);
            if (e && e->distance < min_prev_distance) {
                min_prev_distance = e->distance;
                prev_i = i - 1;
                prev_j = j;
            }

            e = get_element(D_sparse, i, j - 1);
            if (e && e->distance < min_prev_distance) {
                min_prev_distance = e->distance;
                prev_i = i;
                prev_j = j - 1;
            }

            e = get_element(D_sparse, i - 1, j - 1);
            if (e && e->distance < min_prev_distance) {
                min_prev_distance = e->distance;
                prev_i = i - 1;
                prev_j = j - 1;
            }

            D_matrix_element cur;
            cur.distance = d + (min_prev_distance == DBL_MAX ? 0 : min_prev_distance);
            cur.prev_i = prev_i;
            cur.prev_j = prev_j;

            set_element(D_sparse, i, j, cur);
        }
    }

    double min_path_distance = DBL_MAX;
    size_t end_i = 0, end_j = 0;
    bool match = false;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            D_matrix_element *e = get_element(D_sparse, i, j);
            if (!e) continue;

            double cur_path_distance = e->distance + skip_penalty * (n - i + m - j - 2);

            if (cur_path_distance < min_path_distance) {
                min_path_distance = cur_path_distance;
                end_i = i;
                end_j = j;
                match = true;
            }
        }
    }

    ssize_t path_len = 0;

    if (match) {
        D_matrix_element *e;
        *path_distance = min_path_distance;
        for (ssize_t i = end_i, j = end_j; i != -1; ) {
            e = get_element(D_sparse, i, j);
            if (!e) break;
            path_buffer[2 * path_len] = i;
            path_buffer[2 * path_len + 1] = j;
            path_len++;
            ssize_t next_i = e->prev_i;
            ssize_t next_j = e->prev_j;
            i = next_i;
            j = next_j;
        }
        reverse_path(path_buffer, path_len);
    }

    free_sparse_matrix(D_sparse);
    return path_len;
}

double euclid_distance(double *x, double *y, size_t l) {
    double sum = 0;
    for (size_t i = 0; i < l; i++) {
        double v = x[i] - y[i];
        sum += v * v;
    }
    return sqrt(sum);
}

void reverse_path(size_t *path, ssize_t path_len) {
    for (size_t i = 0, j = path_len - 1; i < j; i++, j--) {
        size_t tmp_s = path[2 * i];
        size_t tmp_t = path[2 * i + 1];
        path[2 * i] = path[2 * j];
        path[2 * i + 1] = path[2 * j + 1];
        path[2 * j] = tmp_s;
        path[2 * j + 1] = tmp_t;
    }

    // Helper function to coarsen the sequences for FastDTWBD
double *get_coarsed_sequence(double *s, size_t n, size_t dim) {
    size_t coarse_len = n / 2;
    double *coarse = malloc(coarse_len * dim * sizeof(double));

    for (size_t i = 0; 2 * i + 1 < n; i++) {
        for (size_t j = 0; j < dim; j++) {
            coarse[i * dim + j] = (s[(2 * i) * dim + j] + s[(2 * i + 1) * dim + j]) / 2.0;
        }
    }

    return coarse;
}

// Function to generate window based on coarse path (from FastDTW)
size_t* get_window(size_t n, size_t m, size_t *path_buffer, ssize_t coarse_path_len, int radius) {
    // Adapted from your original `get_window` logic
    size_t *window = malloc(sizeof(size_t) * n * m);
    // Add logic to generate the window from coarse path data
    // ...
    return window;
}

// FastDTWBD function - recursive, multi-resolution DTWBD
ssize_t FastDTWBD(
    double *s, double *t,
    size_t n, size_t m,
    size_t dim,
    double skip_penalty,
    int radius,
    double *path_distance,
    size_t *path_buffer
) {
    const size_t min_len = 2 * (radius + 1) + 1;

    // Base case: use regular DTWBD when sequences are small enough
    if (n < min_len || m < min_len) {
        return dtwbd(s, n, t, m, dim, skip_penalty, NULL, path_buffer, path_distance);
    }

    // Coarsing the sequences for recursive calls
    double *s_coarse = get_coarsed_sequence(s, n, dim);
    double *t_coarse = get_coarsed_sequence(t, m, dim);

    // Recursive call on coarsed sequences
    ssize_t coarse_path_len = FastDTWBD(
        s_coarse, t_coarse, n / 2, m / 2, dim, skip_penalty, radius, path_distance, path_buffer
    );

    // Generate window based on coarse path data
    size_t *window = get_window(n, m, path_buffer, coarse_path_len, radius);

    // Now call regular DTWBD on finer resolution with the generated window
    ssize_t path_len = dtwbd(s, n, t, m, dim, skip_penalty, window, path_buffer, path_distance);

    // Free dynamically allocated memory
    free(s_coarse);
    free(t_coarse);
    free(window);

    return path_len;
}
}