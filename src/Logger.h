
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

// A mutex to lock access to the global log.
static pthread_mutex_t logLock = PTHREAD_MUTEX_INITIALIZER;

// Our logger structure.
typedef struct {
    char date_and_time[51];
    va_list args;
    time_t current_time;
    // Output stream of the logger.
    FILE* file;
} Logger_t;

// Our static global logger.
static Logger_t* logger = NULL;

// Initialize logger.
// Accepts:
//  const char* file -> The name of the output file for log to print to.
//                          If NULL, use stderr.
// Returns: void.
static void _init_log(const char* file_name) {
    pthread_mutex_lock(&logLock);
    // Create the logger.
    logger = (Logger_t*) malloc(sizeof(Logger_t));
    // If no file name was given use stderr.
    if (file_name == NULL)
        logger -> file = stderr;
    else {
        // Open file to log.
        logger -> file = fopen(file_name, "w");
    }
    logger -> current_time = time(NULL);
    tzset();
    pthread_mutex_unlock(&logLock);
}

// Log a message to the logger.
// Accepts:
//  const char* prefix -> The prefix to the message.
//  const char* fmt, ... -> The message to log.
// Returns: void.
static void _log(const char* prefix, const char* fmt, ...) {
    pthread_mutex_lock(&logLock);
    if (logger == NULL) {
        pthread_mutex_unlock(&logLock);
        return;
    }
    // Print local date/time, followed by the prefix, and then, the message.
    strftime(logger -> date_and_time, sizeof(logger -> date_and_time) - 1, "%a %b %d %T %Z %Y", localtime(&logger -> current_time));
    fprintf(logger -> file, "%s: %s ", logger -> date_and_time, prefix);
    va_start(logger -> args, fmt);
    vfprintf(logger -> file, fmt, logger -> args);
    va_end(logger -> args);
    pthread_mutex_unlock(&logLock);
}

// Close the log file and free memory.
// Accepts: void.
// Returns: void.
static void _close_log() {
    pthread_mutex_lock(&logLock);
    if (logger -> file == NULL) {
        pthread_mutex_unlock(&logLock);
        return;
    }
    // If stderr is not used, close the file.
    if (logger -> file != stderr)
        fclose(logger -> file);
    // Free associated memory
    free(logger);
    // Set log back to null.
    logger = NULL;
    pthread_mutex_unlock(&logLock);
}

// Wrappers logging behavior.
#define INIT_LOG(file) _init_log(file)
#define LOG_INFO(...) _log("[INFO]", __VA_ARGS__)
#define LOG_WARNING(...) _log("[WARNING]", __VA_ARGS__)
#define LOG_ERROR(...) _log("[ERROR]", __VA_ARGS__)
#define LOG(...) _log("", __VA_ARGS__)
#define CLOSE_LOG() _close_log()

/*
int main() {
    INIT_LOG("test.log");
    LOG_INFO("Test %d: Info message.\n", 1);
    LOG_WARNING("Test %d: Warning message.\n", 2);
    LOG_ERROR("Test %d: Error message.\n", 3);
    CLOSE_LOG();
}
*/

#endif