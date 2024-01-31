
// File: ThreadPool.h
// Date: 31 January 2024
// Author: TQ Smith
// Purpose: A minimalistic thread pool adapted from Pthreads Programming by Nichols, Buttlar, Farrell.
//          Other resource used: https://nachtimwald.com/2019/04/12/thread-pool-in-c/

#ifndef _THREAD_POOL_
#define _THREAD_POOL_

#include <pthread.h>

#include <stdbool.h>

// Our structure that represents a work request in the queue.
typedef struct ThreadPoolWork {

    // The function for the thread to execute.
    void (*routine) ();
    // The argument to the function.
    void* arg;
    // Pointer to the next work request in the queue.
    struct ThreadPoolWork* next;

} ThreadPoolWork_t;

// Our structure that represents a thread pool.
typedef struct ThreadPool {

    // The number of thread in the pool.
    int numThreads;

    // The maximum number of requests that can be in the pool at a time.
    int maxQueueSize;

    // A flag that indicates the behavior of the main thread when queue is full.
    bool doNotBlockWhenFull;

    // Our array of threads in the pool.
    pthread_t *threads;

    // Current number of request in the queue.
    int curQueueSize;

    // Our front and back pointers to the queue.
    ThreadPoolWork_t* queueHead;
    ThreadPoolWork_t* queueTail;

    // A mutex to access the queue.
    pthread_mutex_t queueLock;

    // Condition variables used to synchronize threads for adding/removing from 
    //  the queue and indicating if the queue is empty.
    pthread_cond_t queueNotEmpty;
    pthread_cond_t queueNotFull;
    pthread_cond_t queueEmpty;

    // Flag to indicate that work in queue and processing threads should be 
    //  finished before the pool is destroyed.
    bool queueClosed;

    // Flag to indicate that the pool is being destoryed.
    bool shutdown;

} ThreadPool_t;

// Our function used to create a thread pool.
// Accepts:
//  int numWorkerThreads -> The number of threads to put into the pool.
//  int maxQueueSize -> The number of maximum requests to put into the queue.
//  bool doNotBlockWhenFull -> Flag to indicate behavior of main thread when queue is full.
// Returns:
//  ThreadPool_t*, A pointer to our created thread pool.
ThreadPool_t* init_thread_pool(int numWorkerThreads, int maxQueueSize, bool doNotBlockWhenFull);

// A function to add a request to the thread pool.
// Accepts:
//  ThreadPool_t* pool -> The thread pool to add the request to.
//  void* rountine -> The routine we want executed.
//  void* arg -> The argument to the routine.
// Returns:
//  int, Returns -1 if doNotBlockWhenFull is set and queue is full. 
//       Otherwise, return 0 if successfully added request.
int thread_pool_add_work(ThreadPool_t* pool, void* routine, void* arg);

// A function used in the main thread to wait for all processing threads
//  and all request in queue to finish.
// Accepts:
//  ThreadPool_t* pool -> The thread pool we are waiting on.
// Returns:
//  void.
void thread_pool_wait(ThreadPool_t* pool);

// A function used to deallocate all of the memory in a thread pool.
// Accepts:
//  ThreadPool_t* -> The thread pool to deallocate.
//  bool finish -> A flag used to indicate if the current threads
//                  should finish processing before being destroyed.
// Returns:
//  void.
void destroy_thread_pool(ThreadPool_t* pool, bool finish);

#endif