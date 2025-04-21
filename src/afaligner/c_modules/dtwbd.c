#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "dtwbd.h"
#include "sparse_matrix.h"
#include "uthash.h"
#include <stdbool.h>
#include "fastdtwbd.h"


#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;

__declspec(dllimport) size_t FastDTWBD();
__declspec(dllimport) size_t DTWBD();
#endif


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

ssize_t FastDTWBD(
    double *s,  // first sequence of MFCC frames – n x l contiguous array
    double *t,  // second sequence of MFCC frames – m x l contiguous array
    size_t n,   // number of frames in first sequence
    size_t  m,   // number of frames in second sequence
    size_t  l,   // number of MFCCs per frame
    double skip_penalty,    // penalty for skipping one frame
    int radius,             // radius of path projection
    double *path_distance,  // place to store warping path distance
    size_t *path_buffer     // buffer to store resulting warping path – (n+m) x 2 contiguous array
) {
    ssize_t path_len;
    size_t min_sequence_len = 2 * (radius + 1) + 1;

    if (n < min_sequence_len || m < min_sequence_len) {
        return DTWBD(s, t, n, m, l, skip_penalty, NULL, path_distance, path_buffer);
    }

    double *coarsed_s = get_coarsed_sequence(s, n, l);
    double *coarsed_t = get_coarsed_sequence(t, m, l);

    path_len = FastDTWBD(coarsed_s, coarsed_t, n/2, m/2, l, skip_penalty, radius, path_distance, path_buffer);

    size_t *window = get_window(n, m, path_buffer, path_len, radius);

    path_len = DTWBD(s, t, n, m, l, skip_penalty, window, path_distance, path_buffer);

    free(coarsed_s);
    free(coarsed_t);
    free(window);

    return path_len;
}


double *get_coarsed_sequence(double *s, size_t n, size_t l) {
    size_t coarsed_sequence_len = n / 2;
    double *coarsed_sequence = malloc(coarsed_sequence_len * l * sizeof(double));

    for (size_t i = 0; 2 * i + 1 < n ; i++) {
        for (size_t j = 0; j < l; j++) {
            coarsed_sequence[l*i+j] = (s[l*(2*i)+j] + s[l*(2*i+1)+j]) / 2;
        }
    }

    return coarsed_sequence;
}


size_t *get_window(size_t n, size_t m, size_t *path_buffer, size_t path_len, int radius) {
    size_t *window = malloc(2*n*sizeof(size_t));

    for (size_t i = 0; i < n; i++) {
        window[2*i] = m;    // maximum value for lower limit
        window[2*i+1] = 0;  // minimum value for upper limit
    }

    for (size_t k = 0; k < path_len; k++) {
        size_t i = path_buffer[2*k];
        size_t j = path_buffer[2*k+1];

        for (ssize_t x = -radius; x < radius + 1; x++) {
            // update lower window limit
            update_window(window, n, m, 2*(i + x), 2*(j - radius));
            update_window(window, n, m, 2*(i + x) + 1, 2*(j - radius));

            // update upper window limit
            update_window(window, n, m, 2*(i + x), 2*(j + radius + 1) + 1);
            update_window(window, n, m, 2*(i + x) + 1, 2*(j + radius + 1) + 1);
        }
    }

    return window;
}


void update_window(size_t *window, size_t n, size_t m, ssize_t i, ssize_t j) {
    if (i < 0 || i >= n) return;

    if (j < 0) {
        j = 0;
    }
    if (j > m - 1) {
        j = m - 1;
    }
    if (j < window[2*i]) {
        window[2*i] = j;
    }
    if (j >= window[2*i+1]) {
        window[2*i+1] = j + 1;
    }
}

}