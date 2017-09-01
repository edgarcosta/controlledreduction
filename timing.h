// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// timing.h: header file for timing and profiling routines

#ifndef TIMING_H
#define TIMING_H

#include <cstdio>
#include <sys/resource.h>

static double timevaldiff_in_seconds(timeval start,
timeval end)
{
  /* Perform the carry for the later subtraction by updating start. */
  if (end.tv_usec < start.tv_usec) {
    int nsec = (start.tv_usec - end.tv_usec) / 1000000 + 1;
    start.tv_usec -= 1000000 * nsec;
    start.tv_sec += nsec;
  }
  if (end.tv_usec - start.tv_usec > 1000000) {
    int nsec = (end.tv_usec - start.tv_usec) / 1000000;
    start.tv_usec += 1000000 * nsec;
    start.tv_sec -= nsec;
  }

  return end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6;
}



#ifdef __APPLE__

#include <sys/time.h>

typedef struct timeval timestamp_type;

static void get_timestamp(timestamp_type *t)
{
  gettimeofday(t, NULL);
}

static double timestamp_diff_in_seconds(timestamp_type start,
timestamp_type end) {
    return timevaldiff_in_seconds(start, end);
}

#else

#include <time.h>

typedef struct timespec timestamp_type;

static void get_timestamp(timestamp_type *t)
{
  clock_gettime(CLOCK_REALTIME, t);
}

static double timestamp_diff_in_seconds(timestamp_type start, timestamp_type end)
{
  struct timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp.tv_sec + 1e-9*temp.tv_nsec;
}

#endif


inline double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}


struct timestamp_pair {
    struct rusage usage;
    timestamp_type wall_time;
} ;

inline void timestamp_mark(timestamp_pair &pair) {
    getrusage(RUSAGE_SELF, &pair.usage);
    get_timestamp( &pair.wall_time );
}
inline void timestamp_report(const timestamp_pair &pair) {
    timestamp_type wtime2;
    get_timestamp(&wtime2);
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    printf("Time: CPU %.3f s, Wall: %.3f s", timevaldiff_in_seconds(pair.usage.ru_utime, usage.ru_utime), timestamp_diff_in_seconds(pair.wall_time, wtime2) );
}
#endif
