#ifndef FASTDTWBD_H
#define FASTDTWBD_H

#include <stddef.h>

ssize_t FastDTWBD(
    const double* s,
    const double* t,
    size_t n,
    size_t m,
    size_t l,
    double skip_penalty,
    int radius,
    double* path_distance,
    size_t* path_buffer
);

#endif
