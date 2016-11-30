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
