/**
 * @file timer.h
 * @brief      A macro that returns the number of seconds that have elapsed
 *             since some point in the past. The timer should return times with
 *             microsecond accuracy.
 *
 * @author     Inanc Senocak
 *
 *             Example:
 *                 #include "timer.h"
 *                 . . .
 *                 double start, finish, elapsed;
 *                 . . .
 *                 GET_TIME(start);
 *                 . . .
 *                 Code to be timed
 *                 . . .
 *                 GET_TIME(finish);
 *                 elapsed = finish - start;
 *                 printf("The code to be timed took %e seconds\n", elapsed);
 *
 * @see        IPP: Section 3.6.1 (p. 121) and Section 6.1.2 (pp. 273 and ff.)
 *
 * @note       The argument now should be a double (not a pointer to a double)
 */
#ifndef _TIMER_H_
#define _TIMER_H_

#include <sys/time.h>

#define GET_TIME(now)                                                                              \
	{                                                                                              \
		struct timeval t;                                                                          \
		gettimeofday(&t, NULL);                                                                    \
		now = t.tv_sec + t.tv_usec / 1000000.0;                                                    \
	}

#endif
