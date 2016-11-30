/**
 * @file process_file.h
 * @brief      Header file for pre/post processing of the NetCDF file for FSM.
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

#ifndef PROCESS_FILE_H
#define PROCESS_FILE_H

#include <stddef.h>
/**
 * @brief      Performs pre-processing of the file for FSM.
 *
 * @param[in]  ncidIn    NetCDF ID of the input file.
 * @param[in]  ncidOut   NetCDF ID of the output file.
 * @param[in]  blockDim  The dimensions of the block grid.
 */
void preProcessFile(int ncidIn, int ncidOut, size_t *blockDim);

/**
 * @brief      Performs post-processing of the file after FSM.
 *
 * @param[in]  ncid  NetCDF ID.
 */
void postProcessFile(int ncid);

#endif /* PROCESS_FILE_H */
