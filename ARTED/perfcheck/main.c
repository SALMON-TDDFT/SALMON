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
#include <stdio.h>
#include <getopt.h>
#define IS_RANGE(X, V, Y) ((X <= V) && (V <= Y))

static const int NLx = 16;
static const int NLy = 16;
static const int NLz = 16;
static const int NK  = 8*8*8;
static const int NB  = 16;
static const int Nt  = 10;

void usage(char const* name)
{
  fprintf(stderr,
          "ARTED 25-points stencil computation performance checker\n\n"
          "[Usage]\n"
          "%s [-h] [-p] [-l NXxNYxNZ] [-k NK] [-b NB] [-s Nt]\n"
          "  -l : set lattice size (default = %dx%dx%d)\n"
          "  -k : set k-point size (default = %d)\n"
          "  -b : set band size    (default %d)\n"
          "  -s : set step size    (default %d)\n"
          "  -h : show this help message\n"
          "  -p : show preprocessor definitions\n",
          name,
          NLx, NLy, NLz,
          NK,
          NB,
          Nt
  );
}

int main(int argc, char* argv[])
{
  int result;

  int nlx = NLx, nly = NLy, nlz = NLz;
  int nk = NK, nb = NB, nt = Nt;

  while((result = getopt(argc, argv, "l:k:b:s:hp")) != -1)
  {
    switch(result)
    {
      case 'l':
        {
          int x_, y_, z_;
          sscanf(optarg, "%dx%dx%d", &x_, &y_, &z_);
          if(!(IS_RANGE(4, x_, 512) && IS_RANGE(4, y_, 512) && IS_RANGE(4, z_, 512))) {
            fprintf(stderr, "lattice size is out of range. (4 <= NLx,y,z <= 512)\n");
            exit(-1);
          }
          nlx = x_;
          nly = y_;
          nlz = z_;
        }
        break;
      case 'k':
        {
          int k_;
          sscanf(optarg, "%d", &k_);
          if(!IS_RANGE(1, k_, (int)pow(2,24))) {
            fprintf(stderr, "k-point size is out of range. (1 <= NK <= 2^24)\n");
            exit(-1);
          }
          nk = k_;
        }
        break;
      case 'b':
        {
          int b_;
          sscanf(optarg, "%d", &b_);
          if(!IS_RANGE(1, b_, (int)pow(2,24))) {
            fprintf(stderr, "k-point size is out of range. (1 <= NK <= 2^24)\n");
            exit(-1);
          }
          nb = b_;
        }
        break;
      case 's':
        {
          int t_;
          sscanf(optarg, "%d", &t_);
          if(t_ < 1) {
            fprintf(stderr, "Nt is invalid value (Nt < 1).\n");
            exit(-1);
          }
          nt = t_;
        }
        break;
      case 'h':
        usage(argv[0]);
        exit(0);
      case 'p':
        print_optimize_message_();
        exit(0);
      case '?':
      default:
        usage(argv[0]);
        exit(-1);
    }
  }

  stencil_perf_check_(&nlx, &nly, &nlz, &nk, &nb, &nt);

  return 0;
}
