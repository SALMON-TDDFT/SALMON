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
#include <stdlib.h>
#include <string.h>

#define ARTED_CPU_TASK_ENV      "ARTED_CPU_TASK_RATIO"       /* 0.1 ~ 1.0 */
#define ARTED_CPU_PPN_ENV       "ARTED_CPU_PPN"              /* Process per Node */
#define ARTED_MIC_PPN_ENV       "ARTED_MIC_PPN"              /* Process per Node */
#define ARTED_LOAD_BALANCER_ENV "ARTED_ENABLE_LOAD_BALANCER" /* 1 or 0 */

void get_cpu_task_ratio_internal_(double * ret) {
  char* env = getenv(ARTED_CPU_TASK_ENV);
  if(env != NULL)
    *ret = atof(env);
  else
    *ret = 1.0;
}

void get_cpu_ppn_internal_(int * ret) {
  char* env = getenv(ARTED_CPU_PPN_ENV);
  if(env != NULL)
    *ret = atoi(env);
  else
    *ret = 1;
}

void get_mic_ppn_internal_(int * ret) {
  char* env = getenv(ARTED_MIC_PPN_ENV);
  if(env != NULL)
    *ret = atoi(env);
  else
    *ret = 1;
}

void get_load_balancer_flag_internal_(int * ret) {
  char* env = getenv(ARTED_LOAD_BALANCER_ENV);
  if(env != NULL)
    *ret = atoi(env);
  else
    *ret = 0;
}
