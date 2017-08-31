// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// timing.h: header file for timing and profiling routines

#ifndef TIMING_H
#define TIMING_H

#ifdef __APPLE__

#include <sys/time.h>

typedef struct timeval timestamp_type;

static void get_timestamp(timestamp_type *t)
{
  gettimeofday(t, NULL);
}

static double timestamp_diff_in_seconds(timestamp_type start,
timestamp_type end)
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
/*static void get_cpu_timestamp(timestamp_type *t)
{
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, t);
}*/

inline double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}


#include <map>
#include <iomanip>
#pragma once
using namespace std;


inline uint64_t cycle_counter()
{
  uint32_t hi, lo;

  __asm __volatile__ ("\t"
		      "rdtsc            \n\t"
		      "movl %%edx, %0   \n\t"
		      "movl %%eax, %1   \n\t"
		      : "=r" (hi), "=r" (lo)
		      :
		      : "%edx", "%eax");

  return (((uint64_t)(hi)) << 32) + ((uint64_t) lo);
}

inline void tp(string call="Rest")
{
	extern map<string,uint64_t> tpmap;
	extern uint64_t timeaux;
	extern string lastcall;
	tpmap[lastcall]+=cycle_counter()-timeaux;
	lastcall=call;
	timeaux=cycle_counter();
}

inline void tpreport()
{
	extern map<string,uint64_t> tpmap;
	extern uint64_t timein;
	timein=cycle_counter()-timein;
	map<string,uint64_t>::const_iterator end=tpmap.end();
	cout<<"       cycles          where"<<endl;
	for(map<string,uint64_t>::iterator it=tpmap.begin();it!=end;++it)
	{
		cout << " " << setw(12) << it->second << setw(7) << setprecision(2) << fixed <<
		floor(it->second*10000/timein)*0.01 << "%  " << it->first << endl;
	}
}

#endif
