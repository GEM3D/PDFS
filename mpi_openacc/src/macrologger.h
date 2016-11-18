/*
 * Copyright (c) 2012 David Rodrigues
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef __MACROLOGGER_H__
#define __MACROLOGGER_H__

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include "mpi.h"

// === auxiliar functions
static inline char *timenow();

#define _FILE strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__

#define NO_LOG 0x00
#define ERROR_LEVEL 0x01
#define INFO_LEVEL 0x02
#define DEBUG_LEVEL 0x03

#ifndef LOG_LEVEL
#define LOG_LEVEL DEBUG_LEVEL
#endif

#define PRINTFUNCTION(format, ...)                                             \
  {                                                                            \
    int mpiRank; char fname[20];                                               \
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);                                   \
    sprintf(fname, "GIN3D-PDFP_%d.log", mpiRank);                              \
    FILE *log = fopen(fname, "a");                                             \
    fprintf(log, format, __VA_ARGS__);                                         \
    fclose(log);                                                               \
  }

#define LOG_FMT "%s | %-5s | %s | %-10s | %s:%4d | "
#define LOG_ARGS(LOG_TAG) timenow(), LOG_TAG, hostname(), _FILE, __FUNCTION__, __LINE__

#define NEWLINE "\n"

#define ERROR_TAG "ERROR"
#define INFO_TAG "INFO"
#define DEBUG_TAG "DEBUG"

#if LOG_LEVEL >= DEBUG_LEVEL
#define LOG_DEBUG(message, args...)                                            \
  PRINTFUNCTION(LOG_FMT message NEWLINE, LOG_ARGS(DEBUG_TAG), ##args)
#else
#define LOG_DEBUG(message, args...)
#endif

#if LOG_LEVEL >= INFO_LEVEL
#define LOG_INFO(message, args...)                                             \
  PRINTFUNCTION(LOG_FMT message NEWLINE, LOG_ARGS(INFO_TAG), ##args)
#else
#define LOG_INFO(message, args...)
#endif

#if LOG_LEVEL >= ERROR_LEVEL
#define LOG_ERROR(message, args...)                                            \
  PRINTFUNCTION(LOG_FMT message NEWLINE, LOG_ARGS(ERROR_TAG), ##args)
#else
#define LOG_ERROR(message, args...)
#endif

#if LOG_LEVEL >= NO_LOGS
#define LOG_IF_ERROR(condition, message, args...)                              \
  if (condition)                                                               \
  PRINTFUNCTION(LOG_FMT message NEWLINE, LOG_ARGS(ERROR_TAG), ##args)
#else
#define LOG_IF_ERROR(condition, message, args...)
#endif

static inline char *timenow() {
  static char buffer[64];
  time_t rawtime;
  struct tm *timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);

  strftime(buffer, 64, "%Y-%m-%d %H:%M:%S", timeinfo);

  return buffer;
}

static inline char *hostname() {
  static char buffer[8];

  gethostname(buffer, sizeof(buffer));

  return buffer;
}
#endif
