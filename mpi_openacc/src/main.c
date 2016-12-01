/**
 * @file main.c
 * @brief      Main program file for Parallel Distance Field Solver.
 *             "Massively Parallel Algorithm for Solving the Eikonal Equation
 *             on Multiple Accelerator Platforms"
 *             M.S. in Computer Science
 *             Boise State University
 *
 *             Required Libraries:
 *              * MPI
 *              * NetCDF-4 | HDF5 | SZIP
 *              * OpenACC
 *             Usage:
 *              <progName> -i <filename.nc> -p <outputPrefix> [--nx <val>]
 *              [--ny <val>] [--nz <val>]
 *
 * @author     Shrestha, Anup
 * @date       12 JUL 2016
 *
 * Copyright (c) 2016
 * Mechanical and Bio-medical Engineering Department
 * Boise State University
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "block_fsm.h"
#include "macrologger.h"
#include "ncwrap.h"
#include "phi3D.h"
#include "process_file.h"
#include "utilities.h"

#include <getopt.h>
#include <inttypes.h>
#include <mpi.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* MPI routine variables */
static int mpiSize;
static int mpiRank;
/*-----------------------------------------*/

/* Private method definitions */
static void pdfpStart(char *in, char *out, size_t *blockDim);
static void mprintf(const char *format, ...);
static int checkArgs(int argc, char *argv[], char **ipfile, char *opfile, size_t *dim);
/*-----------------------------------------*/

int main(int argc, char *argv[])
{
  /* Initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

  /* MPI Error Handler */
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  char * ipfile;                      /* Input filename */
  char   opfile[255];                 /* Output filename */
  size_t blockDim[NDIMS] = {1, 1, 1}; /* New Dimension of decomposed domain */

  if (checkArgs(argc, argv, &ipfile, opfile, blockDim)) {
    mprintf("\n------------------------\n");
    mprintf(" Input: %s\n", ipfile);
    mprintf("Output: %s\n", opfile);
    mprintf("------------------------\n");
    mprintf("nx: %d, ny: %d, nz: %d\n", blockDim[0], blockDim[1], blockDim[2]);
    mprintf("------------------------\n");

    pdfpStart(ipfile, opfile, blockDim);

  } else {
    mprintf("\nUsage: %s -i <filename.nc> -p <outputPrefix> [--nx <val>] [--ny "
            "<val>] [--nz <val>]\n\n",
            argv[0]);
  }

  /* Shut down MPI. */
  MPI_Finalize();

  return 0;
}

/**
 * @brief      Starts the Parallel Distance Field Preprocessor.
 *
 * @param[in]  in        Input filename.
 * @param[in]  out       Output filename.
 * @param[in]  blockDim  The dimensions of the block grid.
 */
static void pdfpStart(char *in, char *out, size_t *blockDim)
{
  int ncidIn;
  int ncidOut;

  ncOpen(in, &ncidIn);
  ncCreate(out, &ncidOut);

  LOG_DEBUG("[%d/%d] Pre Processing ....", mpiRank, mpiSize);
  preProcessFile(ncidIn, ncidOut, blockDim);

  ncClose(ncidIn);

  MPI_Barrier(MPI_COMM_WORLD);

  LOG_DEBUG("[%d/%d] Block FSM Start ...", mpiRank, mpiSize);
  startBlockFSM(ncidOut, blockDim);
  LOG_DEBUG("[%d/%d] Block FSM Complete.", mpiRank, mpiSize);

  MPI_Barrier(MPI_COMM_WORLD);

  LOG_DEBUG("[%d/%d] Post Processing ...", mpiRank, mpiSize);
  postProcessFile(ncidOut);

  LOG_DEBUG("[%d/%d] Closing ...........", mpiRank, mpiSize);
  ncClose(ncidOut);
}

/**
 * @brief      Master printf.
 *
 *             Print only when rank is the master (0) rank.
 *
 * @param[in]  format     Formatted string to print.
 * @param[in]  <unnamed>  Variable length arguments.
 */
static void mprintf(const char *format, ...)
{
  if (mpiRank == 0) {
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
  }
}

/**
 * @brief      Check the arguments.
 *
 *             Validates the arguments passed through the command line.
 *
 * @param[in]  argc    Argument count
 * @param[in]  argv    Argument values
 * @param[out] ipfile  Pointer location to store input filename.
 * @param[out] opfile  Pointer location to store output filename.
 *
 * @return     1 if successful else 0
 */
static int checkArgs(int argc, char *argv[], char **ipfile, char *opfile, size_t *dim)
{
  /* Mandatory option flags */
  int iflag = 0;
  int pflag = 0;
  int opt;

  static struct option lopts[] = {{"input", required_argument, NULL, 'i'},
                                  {"prefix", required_argument, NULL, 'p'},
                                  {"nx", required_argument, NULL, 'x'},
                                  {"ny", required_argument, NULL, 'y'},
                                  {"nz", required_argument, NULL, 'z'}};

  while ((opt = getopt_long(argc, argv, "i:p:x:y:z", lopts, NULL)) != -1) {
    switch (opt) {
      case 'i':
        *ipfile = optarg;
        iflag   = 1;
        break;
      case 'p':
        sprintf(opfile, "%s_soln.nc", optarg);
        pflag = 1;
        break;
      case 'x':
        dim[0] = strtoumax(optarg, NULL, 10);
        if (dim[0] <= 0) return 0;
        break;
      case 'y':
        dim[1] = strtoumax(optarg, NULL, 10);
        if (dim[1] <= 0) return 0;
        break;
      case 'z':
        dim[2] = strtoumax(optarg, NULL, 10);
        if (dim[2] <= 0) return 0;
        break;
      case '?':
        break;
      default:
        mprintf("Invalid option.\n");
        return 0;
    }
  }

  if (mpiSize != (dim[0] * dim[1] * dim[2])) {
    mprintf("Total processes must equal total blocks in the decmposition!\n");
    return 0;
  }

  if (iflag == 0 || pflag == 0) {
    return 0;
  }

  return 1;
}