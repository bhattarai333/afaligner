ssize_t DTWBD(double *s, double *t, size_t n, size_t m, size_t l,
              double skip_penalty, size_t *window, double *path_distance,
              size_t *path_buffer) {
    log_function_entry("DTWBD");
    log_info("Parameters: n=%zu, m=%zu, l=%zu, skip_penalty=%f",
             n, m, l, skip_penalty);

    D_matrix_element *D_matrix = malloc(sizeof(D_matrix_element) * n * m);

    if (D_matrix == NULL) {
        log_error("malloc() failed when allocating D_matrix");
        return -1;
    }

    // ... rest of the function ...

    if (match) {
        log_info("Match found: end_i=%zu, end_j=%zu, min_path_distance=%f",
                end_i, end_j, min_path_distance);
    } else {
        log_warn("No match found");
    }

    free(D_matrix);
    log_function_exit("DTWBD", path_len);
    return path_len;
}