/**
 * @file process_file.h
 * @brief      Header file for pre/post processing of the NetCDF file for FSM.
 *
 * @author     Shrestha, Anup
 * @date       12 JUL 2016
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
