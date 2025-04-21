#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <stdio.h>

#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;

__declspec(dllimport) size_t FastDTWBD();
__declspec(dllimport) size_t DTWBD();
#endif

typedef struct {
    double distance;
    ssize_t prev_i;
    ssize_t prev_j;
} D_matrix_element;

// Function declarations
ssize_t FastDTWBD(double *s, double *t, size_t n, size_t m, size_t l,
                  double skip_penalty, int radius, double *path_distance,
                  size_t *path_buffer);

double *get_coarsed_sequence(double *s, size_t n, size_t l);
size_t *get_window(size_t n, size_t m, size_t *path_buffer, size_t path_len, int radius);
void update_window(size_t *window, size_t n, size_t m, ssize_t i, ssize_t j);
double euclid_distance(double *x, double *y, size_t l);
double get_distance(D_matrix_element *D_matrix, size_t *index_map, size_t n, size_t m, 
                   size_t *window, size_t i, size_t j);
D_matrix_element get_best_candidate(D_matrix_element *candidates, size_t n);
void reverse_path(size_t *path, ssize_t path_len);

ssize_t DTWBD(
    double *s,  // first sequence of MFCC frames – n x l contiguous array
    double *t,  // second sequence of MFCC frames – m x l contiguous array
    size_t n,   // number of frames in first sequence
    size_t m,   // number of frames in second sequence
    size_t l,   // number of MFCCs per frame
    double skip_penalty,    // penalty for skipping one frame
    size_t *window,        // n x 2 contiguous array
    double *path_distance, // place to store warping path distance
    size_t *path_buffer   // buffer to store resulting warping path
) {
    // Calculate total elements needed within window
    size_t total_elements = 0;
    for (size_t i = 0; i < n; i++) {
        size_t from = window == NULL ? 0 : window[2*i];
        size_t to = window == NULL ? m : window[2*i+1];
        total_elements += (to - from);
    }

    // Allocate only needed elements
    D_matrix_element *D_matrix = malloc(sizeof(D_matrix_element) * total_elements);
    if (D_matrix == NULL) {
        fprintf(stderr, "ERROR: malloc() failed when allocating D_matrix\n");
        return -1;
    }

    // Create index mapping array
    size_t *index_map = malloc(sizeof(size_t) * n * m);
    if (index_map == NULL) {
        free(D_matrix);
        fprintf(stderr, "ERROR: malloc() failed when allocating index_map\n");
        return -1;
    }

    // Initialize index mapping
    size_t current_index = 0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            index_map[i*m + j] = (size_t)-1;  // Invalid index
        }
        size_t from = window == NULL ? 0 : window[2*i];
        size_t to = window == NULL ? m : window[2*i+1];
        for (size_t j = from; j < to; j++) {
            index_map[i*m + j] = current_index++;
        }
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
                { get_distance(D_matrix, index_map, n, m, window, i-1, j-1) + d, i-1, j-1 },
                { get_distance(D_matrix, index_map, n, m, window, i, j-1) + d, i, j-1 },
                { get_distance(D_matrix, index_map, n, m, window, i-1, j) + d, i-1, j },
            };

            size_t idx = index_map[i*m+j];
            D_matrix[idx] = get_best_candidate(candidates, sizeof(candidates)/sizeof(D_matrix_element));

            cur_path_distance = D_matrix[idx].distance + skip_penalty * (n - i + m - j - 2);

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
        ssize_t i = end_i;
        ssize_t j = end_j;
        
        while (i != -1) {
            size_t idx = index_map[i*m+j];
            e = &D_matrix[idx];
            path_buffer[2*path_len] = i;
            path_buffer[2*path_len+1] = j;
            path_len++;
            i = e->prev_i;
            j = e->prev_j;
        }
        reverse_path(path_buffer, path_len);
    }

    free(index_map);
    free(D_matrix);

    return path_len;
}

double get_distance(D_matrix_element *D_matrix, size_t *index_map, size_t n, size_t m, 
                   size_t *window, size_t i, size_t j) {
    if (i < 0 || i >= n || j < 0 || 

j >= m) {
        return DBL_MAX;
    }

    if (window == NULL || 

(j >= window[2*i] && j < window[2*i+1])) {
        size_t idx = index_map[i*m + j];
        if (idx == (size_t)-1) return DBL_MAX;
        return D_matrix[idx].distance;
    }

    return DBL_MAX;
}

// Rest of the functions remain unchanged
ssize_t FastDTWBD(double *s, double *t, size_t n, size_t m, size_t l,
                  double skip_penalty, int radius, double *path_distance,
                  size_t *path_buffer) {
    ssize_t path_len;
    size_t min_sequence_len = 2 * (radius + 1) + 1;

    if (n < min_sequence_len || 

m < min_sequence_len) {
        return DTWBD(s, t, n, m, l, skip_penalty, NULL, path_distance, path_buffer);
    }

    double *coarsed_s = get_coarsed_sequence(s, n, l);
    double *coarsed_t = get_coarsed_sequence(t, m, l);

    path_len = FastDTWBD(coarsed_s, coarsed_t, n/2, m/2, l, skip_penalty, radius, 
                        path_distance, path_buffer);

    size_t *window = get_window(n, m, path_buffer, path_len, radius);

    path_len = DTWBD(s, t, n, m, l, skip_penalty, window, path_distance, path_buffer);

    free(coarsed_s);
    free(coarsed_t);
    free(window);

    return path_len;
}

// Remaining helper functions stay the same
double euclid_distance(double *x, double *y, size_t l) {
    double sum = 0;
    for (size_t i = 0; i < l; i++) {
        double v = x[i] - y[i];
        sum += v * v;
    }
    return sqrt(sum);
}

D_matrix_element get_best_candidate(D_matrix_element *candidates, size_t n) {
    double min_distance = DBL_MAX;
    D_matrix_element best_candidate;

    for (size_t i = 0; i < n; i++) {
        if (candidates[i].distance < min_distance) {
            min_distance = candidates[i].distance;
            best_candidate = candidates[i];
        }
    }

    return best_candidate;
}

void reverse_path(size_t *path, ssize_t path_len) {
    for (size_t i = 0, j = path_len - 1; i < j; i++, j--) {
        size_t tmp_s = path[2*i];
        size_t tmp_t = path[2*i+1];
        path[2*i] = path[2*j];
        path[2*i+1] = path[2*j+1];
        path[2*j] = tmp_s;
        path[2*j+1] = tmp_t;
    }
}
