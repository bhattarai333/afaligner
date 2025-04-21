#include "helper.h"

ssize_t DTWBD(
    double *s,
    double *t,
    size_t n,
    size_t m,
    size_t l,
    double skip_penalty,
    size_t *window,
    double *path_distance,
    size_t *path_buffer
) {
    D_matrix_element *D_matrix = malloc(sizeof(D_matrix_element) * n * m);

    if (D_matrix == NULL) {
        fprintf(stderr, "ERROR: malloc() failed when allocating D_matrix\n");
        return -1;
    }

    double min_path_distance;
    double cur_path_distance;
    size_t end_i, end_j;
    bool match = false;

    min_path_distance = skip_penalty * (n + m);

    for (size_t i = 0; i < n; i++) {
        size_t from = window == NULL ? 0 : window[2*i];
        size_t to = window == NULL ? m : window[2*i+1];
        for (size_t j = from; j < to; j++) {
            double d = euclid_distance(s+i*l, t+j*l, l);

            D_matrix_element candidates[] = {
                { skip_penalty * (i + j) + d, -1, -1 },
                { get_distance(D_matrix, n, m, window, i-1, j-1) + d, i-1, j-1 },
                { get_distance(D_matrix, n, m, window, i, j-1) + d, i, j-1 },
                { get_distance(D_matrix, n, m, window, i-1, j) + d, i-1, j },
            };

            D_matrix[i*m+j] = get_best_candidate(candidates, sizeof(candidates)/sizeof(D_matrix_element));

            cur_path_distance = D_matrix[i*m+j].distance + skip_penalty * (n - i + m - j - 2);

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
        for (ssize_t i = end_i, j = end_j; i != -1; i = e->prev_i, j = e->prev_j) {
            e = &D_matrix[i*m+j];
            path_buffer[2*path_len] = i;
            path_buffer[2*path_len+1] = j;
            path_len++;
        }
        reverse_path(path_buffer, path_len);
    }

    free(D_matrix);

    return path_len;
}