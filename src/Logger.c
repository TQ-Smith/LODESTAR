
// File: Logger.c
// Date: 3 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: A simple thread safe logger.
// Reference: https://tuttlem.github.io/2012/12/08/simple-logging-in-c.html

#include "Logger.h"

Logger_t* logger = NULL;

pthread_mutex_t logLock = PTHREAD_MUTEX_INITIALIZER;

void _init_log(const char* file_name) {
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

void _log(const char* prefix, const char* fmt, ...) {
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

void _close_log() {
    pthread_mutex_lock(&logLock);
    if (logger == NULL) {
        pthread_mutex_unlock(&logLock);
        return;
    }
    // If stderr is not used, close the file.
    if (logger -> file != stderr)
        fclose(logger -> file);
    free(logger);
    logger = NULL;
    pthread_mutex_unlock(&logLock);
}

/*
int main() {
    INIT_LOG("test.log");
    LOG_INFO("Test %d: Info message.\n", 1);
    LOG_WARNING("Test %d: Warning message.\n", 2);
    LOG_ERROR("Test %d: Error message.\n", 3);
    CLOSE_LOG();
}
*/