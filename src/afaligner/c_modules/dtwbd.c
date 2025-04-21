#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <stdio.h>

#ifdef _WIN32
    #define EXPORT __declspec(dllexport)
#else
    #define EXPORT __attribute__((visibility("default")))
#endif

typedef struct {
    double distance;
    ssize_t prev_i;
    ssize_t prev_j;
} D_matrix_element;

// Function declarations with EXPORT
EXPORT ssize_t FastDTWBD(double *s, double *t, size_t n, size_t m, size_t l,
                        double skip_penalty, int radius, double *path_distance,
                        size_t *path_buffer);

EXPORT double *get_coarsed_sequence(double *s, size_t n, size_t l);

EXPORT size_t *get_window(size_t n, size_t m, size_t *path_buffer, size_t path_len, int radius);

EXPORT double euclid_distance(double *x, double *y, size_t l);

EXPORT ssize_t DTWBD(
    double *s, double *t, size_t n, size_t m, size_t l,
    double skip_penalty, size_t *window,
    double *path_distance, size_t *path_buffer
) {
    // Only store three rows at a time: previous, current, and next
    const size_t ROWS_TO_STORE = 3;
    size_t max_row_elements = m;  // Maximum width of a row
    
    // Allocate memory for just 3 rows
    D_matrix_element *D_matrix = malloc(sizeof(D_matrix_element) * ROWS_TO_STORE * max_row_elements);
    if (D_matrix == NULL) {
        fprintf(stderr, "ERROR: malloc() failed when allocating D_matrix\n");
        return -1;
    }

    // Store previous results needed for backtracking
    size_t *prev_indices = malloc(sizeof(size_t) * n * 2);  // Store (i,j) pairs
    double *prev_distances = malloc(sizeof(double) * n);
    if (prev_indices == NULL || 

prev_distances == NULL) {
        free(D_matrix);
        free(prev_indices);
        free(prev_distances);
        fprintf(stderr, "ERROR: malloc() failed\n");
        return -1;
    }

    double min_path_distance = skip_penalty * (n + m);
    double cur_path_distance;
    size_t end_i = 0, end_j = 0;
    bool match = false;

    // Process matrix row by row
    for (size_t i = 0; i < n; i++) {
        size_t from = window == NULL ? 0 : window[2*i];
        size_t to = window == NULL ? m : window[2*i+1];
        size_t current_row = i % ROWS_TO_STORE;
        
        for (size_t j = from; j < to; j++) {
            double d = euclid_distance(s+i*l, t+j*l, l);
            
            // Get values from previous rows (stored in circular buffer)
            double prev_diag = (i > 0 && j > 0) ? 
                D_matrix[((i-1) % ROWS_TO_STORE)*m + (j-1)].distance : DBL_MAX;
            double prev_up = (i > 0) ? 
                D_matrix[((i-1) % ROWS_TO_STORE)*m + j].distance : DBL_MAX;
            double prev_left = (j > 0) ? 
                D_matrix[current_row*m + (j-1)].distance : DBL_MAX;
            
            // Calculate current cell
            D_matrix_element *current = &D_matrix[current_row*m + j];
            double start_cost = skip_penalty * (i + j) + d;
            
            if (start_cost <= prev_diag && start_cost <= prev_up && start_cost <= prev_left) {
                current->distance = start_cost;
                current->prev_i = -1;
                current->prev_j = -1;
            } else if (prev_diag <= prev_up && prev_diag <= prev_left) {
                current->distance = prev_diag + d;
                current->prev_i = i-1;
                current->prev_j = j-1;
            } else if (prev_up <= prev_left) {
                current->distance = prev_up + d;
                current->prev_i = i-1;
                current->prev_j = j;
            } else {
                current->distance = prev_left + d;
                current->prev_i = i;
                current->prev_j = j-1;
            }

            cur_path_distance = current->distance + skip_penalty * (n - i + m - j - 2);
            if (cur_path_distance < min_path_distance) {
                min_path_distance = cur_path_distance;
                end_i = i;
                end_j = j;
                match = true;
                
                prev_indices[i*2] = current->prev_i;
                prev_indices[i*2+1] = current->prev_j;
                prev_distances[i] = current->distance;
            }
        }
    }

    // Reconstruct path
    ssize_t path_len = 0;
    if (match) {
        *path_distance = min_path_distance;
        size_t current_i = end_i;
        size_t current_j = end_j;
        
        size_t *temp_path = malloc(sizeof(size_t) * (n+m) * 2);
        size_t temp_len = 0;
        
        while (current_i != (size_t)-1) {
            temp_path[temp_len*2] = current_i;
            temp_path[temp_len*2+1] = current_j;
            temp_len++;
            
            size_t next_i = prev_indices[current_i*2];
            size_t next_j = prev_indices[current_i*2+1];
            current_i = next_i;
            current_j = next_j;
        }
        
        for (size_t i = 0; i < temp_len; i++) {
            path_buffer[i*2] = temp_path[(temp_len-1-i)*2];
            path_buffer[i*2+1] = temp_path[(temp_len-1-i)*2+1];
        }
        path_len = temp_len;
        
        free(temp_path);
    }

    free(D_matrix);
    free(prev_indices);
    free(prev_distances);

    return path_len;
}

EXPORT ssize_t FastDTWBD(double *s, double *t, size_t n, size_t m, size_t l,
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

    if (coarsed_s == NULL || 

coarsed_t == NULL) {
        free(coarsed_s);
        free(coarsed_t);
        return -1;
    }

    path_len = FastDTWBD(coarsed_s, coarsed_t, n/2, m/2, l, skip_penalty, radius, 
                        path_distance, path_buffer);

    if (path_len < 0) {
        free(coarsed_s);
        free(coarsed_t);
        return -1;
    }

    size_t *window = get_window(n, m, path_buffer, path_len, radius);
    if (window == NULL) {
        free(coarsed_s);
        free(coarsed_t);
        return -1;
    }

    path_len = DTWBD(s, t, n, m, l, skip_penalty, window, path_distance, path_buffer);

    free(coarsed_s);
    free(coarsed_t);
    free(window);

    return path_len;
}

// Helper functions
EXPORT double *get_coarsed_sequence(double *s, size_t n, size_t l) {
    double *coarsed = malloc(sizeof(double) * ((n+1)/2) * l);
    if (coarsed == NULL) return NULL;
    
    for (size_t i = 0; i < n/2; i++) {
        for (size_t j = 0; j < l; j++) {
            coarsed[i*l + j] = (s[2*i*l + j] + s[(2*i+1)*l + j]) / 2;
        }
    }
    if (n % 2) {
        for (size_t j = 0; j < l; j++) {
            coarsed[(n/2)*l + j] = s[(n-1)*l + j];
        }
    }
    return coarsed;
}

EXPORT size_t *get_window(size_t n, size_t m, size_t *path_buffer, size_t path_len, int radius) {
    size_t *window = malloc(sizeof(size_t) * n * 2);
    if (window == NULL) return NULL;

    // Initialize window bounds
    for (size_t i = 0; i < n; i++) {
        window[2*i] = m;     // from (initialize to maximum)
        window[2*i+1] = 0;   // to (initialize to minimum)
    }

    // Expand path points to create window
    for (size_t i = 0; i < path_len; i++) {
        size_t path_i = path_buffer[2*i] * 2;
        size_t path_j = path_buffer[2*i+1];

        // Update window for nearby rows within radius
        for (ssize_t r = -radius; r <= radius; r++) {
            ssize_t curr_i = path_i + r;
            if (curr_i < 0 || 

curr_i >= (ssize_t)n) continue;

            ssize_t min_j = (ssize_t)path_j - radius;
            ssize_t max_j = (ssize_t)path_j + radius + 1;
            
            if (min_j < 0) min_j = 0;
            if (max_j > (ssize_t)m) max_j = m;

            if (min_j < (ssize_t)window[curr_i]) window[curr_i] = min_j;
            if (max_j > (ssize_t)window[curr_i+1]) window[curr_i+1] = max_j;
        }
    }

    return window;
}

EXPORT double euclid_distance(double *x, double *y, size_t l) {
    double sum = 0;
    for (size_t i = 0; i < l; i++) {
        double diff = x[i] - y[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}
