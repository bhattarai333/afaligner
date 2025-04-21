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
}