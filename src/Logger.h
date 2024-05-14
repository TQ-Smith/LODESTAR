
// File: Logger.h
// Date: 3 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: A simple thread safe logger.
// Reference: https://tuttlem.github.io/2012/12/08/simple-logging-in-c.html

#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

// Our logger structure.
typedef struct {
    char date_and_time[51];
    va_list args;
    time_t current_time;
    // Output stream of the logger.
    FILE* file;
} Logger_t;

// Our global logger.
extern Logger_t* logger;

// Our mutex to access logger.
extern pthread_mutex_t logLock;

// Initialize logger.
// Accepts:
//  const char* file -> The name of the output file for log to print to.
//                          If NULL, use stderr.
// Returns: void.
void _init_log(const char* file_name);

// Log a message to the logger.
// Accepts:
//  const char* prefix -> The prefix to the message.
//  const char* fmt, ... -> The message to log.
// Returns: void.
void _log(const char* prefix, const char* fmt, ...);

// Close the log file and free memory.
// Accepts: void.
// Returns: void.
void _close_log();

// Wrappers logging behavior.
#define INIT_LOG(file) _init_log(file)
#define LOG_INFO(...) _log("[INFO]", __VA_ARGS__)
#define LOG_WARNING(...) _log("[WARNING]", __VA_ARGS__)
#define LOG_ERROR(...) _log("[ERROR]", __VA_ARGS__)
#define LOG(...) _log("", __VA_ARGS__)
#define CLOSE_LOG() _close_log()

#endif