/**
 * @file phi3D.c
 * @brief      Source file for 3D Phi function.
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

#include "phi3D.h"

#include <stdio.h>
#include <stdlib.h>
/**
 * @brief      Create the Phi3D data structure.
 *
 *             Assigns proper attributes to the Phi3D struct and allocates
 *             memory for the location and distance attribute. Returns 0 if the
 *             memory allocation failed.
 *
 * @note       This function makes the assumption that the dimension length
 *             includes the exterior nodes and therefore subtracts 2 to assign
 *             the dimension length of the interior nodes. The index 0,1,2 for
 *             dim param is the z,y,x-dimension.
 *
 * @param[in]  dim   Array of dimension length.
 * @param[in]  dlt   Array of delta values.
 * @param[out] pf    Pointer to location of Phi function.
 *
 * @return     1 - Success; 0 - Failure
 */
int makePhi3D(size_t *dim, double *dlt, Phi *pf)
{
	pf->x = dim[2] - 2;
	pf->y = dim[1] - 2;
	pf->z = dim[0] - 2;

	pf->dx = dlt[0];
	pf->dy = dlt[1];
	pf->dz = dlt[2];

	pf->F = SPEED;

	int totalNodes = dim[0] * dim[1] * dim[2];
	pf->distance   = (double *) malloc(sizeof(double) * totalNodes);
	if (pf->distance == NULL) {
		printf("Error allocating memory: %s:%d\n", __FILE__, __LINE__);
		return 0;
	}

	return 1;
}

/**
 * @brief      Destroy the Phi3D data structure.
 *
 *             Deallocates memory allocated by the makePhi3D function.
 *
 * @param[out] pf    Pointer to location of Phi function.
 *
 * @return     1 - Success; 0 - Failure
 */
int destroyPhi3D(Phi *pf)
{
	free(pf->distance);

	return 1;
}
