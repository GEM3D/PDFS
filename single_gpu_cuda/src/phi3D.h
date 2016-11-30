/**
 * @file phi3D.h
 * @brief      Header file for 3D Phi function.
 *
 * @author     Shrestha, Anup
 * @date       09 OCT 2015
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
#ifndef PHI3D_H_
#define PHI3D_H_

#include <math.h>

#include "file_writer.h"

#define SPEED 1 /* Speed of propagation [Equation (1) in the report] */
#define DEFAULT_BORDER_LOCATION -1
#define DEFAULT_BORDER_DISTANCE INFINITY
#define DEFAULT_INTERIOR_DISTANCE 90000

/**
 * @brief      Structure of the Phi3D datatype.
 */
typedef struct {
	int     x;        /**< Length of x-dimension*/
	int     y;        /**< Length of y-dimension*/
	int     z;        /**< Length of z-dimension*/
	double  dx;       /**< Node spacing in the x-direction */
	double  dy;       /**< Node spacing in the y-direction */
	double  dz;       /**< Node spacing in the z-direction */
	int *   location; /**< Location values array */
	double *distance; /**< Distance values array */
	double  F;        /**< Speed of propagation*/
} Phi;

/**
 * @brief      Creates a phi function.
 *
 * @param[in]  f  Pointer to the VTI file.
 * @param[out] pf Pointer to phi function.
 */
void create_phi_function(FILE *f, Phi *pf);

/**
 * @brief      Deallocates memory allocated for the phi function
 *
 * @param[out] pf    Pointer to phi function.
 */
void destroy_phi_function(Phi *pf);

/**
 * @brief      Calculates the distance field.
 *
 * @param[in,out] pf    Pointer to phi function.
 */
void calc_dist_field(Phi *pf);

/**
 * @brief      Generates the output file for the distance field.
 *
 * @param[in]  pf     Pointer to phi function.
 * @param[in]  out    The type of output file.
 * @param[in]  fname  Filename of the output file.
 */
void phi_write_file(Phi *pf, FileOut out, char *fname);

#endif /* PHI3D_H_ */
