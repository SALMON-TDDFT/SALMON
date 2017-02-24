/*
 *  Copyright 2016 ARTED developers
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */
#ifdef ARTED_USE_PAPI

#include <papi.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>

int    *EventSet;
double values[2];

void papi_begin_() {
  int ret, i;

  if(PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
    fprintf(stderr, "PAPI library init error!\n");
    exit(1);
  }
  if(PAPI_thread_init((unsigned long (*)(void))(omp_get_num_threads)) != PAPI_OK) {
    fprintf(stderr, "PAPI thread init error.\n");
    exit(1);
  }
  if (PAPI_num_counters() < 2) {
    fprintf(stderr, "No hardware counters here, or PAPI not supported.\n");
    exit(1);
  }

  EventSet = (int*) malloc(sizeof(int) * omp_get_max_threads());
  for(i = 0 ; i < omp_get_max_threads() ; EventSet[i++] = PAPI_NULL);

#pragma omp parallel
  {
    int t = omp_get_thread_num();
    PAPI_create_eventset(&EventSet[t]);
    PAPI_add_event(EventSet[t], PAPI_SP_OPS);
    PAPI_add_event(EventSet[t], PAPI_DP_OPS);

    if ((ret = PAPI_start(EventSet[t])) != PAPI_OK) {
      fprintf(stderr, "PAPI failed to start counters: %s\n", PAPI_strerror(ret));
      exit(1);
    }
  }

  for(i = 0 ; i < 2 ; values[i++] = 0);
}

void papi_end_() {
  int ret, i;
  long long v[2];
  long long v0, v1;
  double vin[2];

  v0 = v1 = 0;
#pragma omp parallel shared(v) reduction(+:v0,v1)
  {
    int t = omp_get_thread_num();
    if ((ret = PAPI_stop(EventSet[t],v)) != PAPI_OK) {
      fprintf(stderr, "PAPI failed to read counters: %s\n", PAPI_strerror(ret));
      exit(1);
    }
    v0 += v[0];
    v1 += v[1];
  }
  vin[0] = v0;
  vin[1] = v1;

  MPI_Reduce(vin, values, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  PAPI_shutdown();

  free(EventSet);
}

void papi_result_(double *time) {
  printf("SP FLOP = %f\n", values[0]);
  printf("DP FLOP = %f\n", values[1]);
  printf("GFLOPS  = %.2f\n", ((values[0] + values[1]) / *time) * 1.0e-9);
}

#else

void papi_begin_()              {}
void papi_end_()                {}
void papi_result_(double* time) {}

#endif

