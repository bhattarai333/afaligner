#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>
#include <time.h>

// Log levels
typedef enum {
    LOG_DEBUG,
    LOG_INFO,
    LOG_WARN,
    LOG_ERROR
} LogLevel;

// Initialize logger
int init_logger(const char* log_path);

// Close logger
void close_logger(void);

// Log messages with different levels
void log_debug(const char* format, ...);
void log_info(const char* format, ...);
void log_warn(const char* format, ...);
void log_error(const char* format, ...);


// Function entry/exit logging
void log_function_entry(const char* function_name);
void log_function_exit(const char* function_name, ssize_t result);

#endif