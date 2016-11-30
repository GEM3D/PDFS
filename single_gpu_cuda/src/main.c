/**
 * @file main.c
 * @brief      Main program file for Parallel Distance Field Preprocessor.
 *
 *             Calculates the distance field for an interface in a
 *             given VTI file by solving the Eikonal equation using
 *             the Fast Sweeping Method.
 *
 *             This version implements the parallel fast sweeping using
 *             Cuthill-Mckee ordering as described in the paper
 *             A parallel fast sweeping method for the Eikonal equation"
 *             by Miles Detrixhe, Federic Gibou, and Chohong Min.
 *
 *             Required Libraries:
 *              * CUDA
 *             Usage:
 *              <progName> <filename.vti> <outputPrefix>
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

#include "phi3D.h"

// Global variables
char *fileToRead;       // name of the file to read
char *prefix;           // prefix for the output file
char  outFilename[255]; // filename for the output file

int main(int argc, char *argv[])
{
  if (argc != 3) {
    printf("Usage: %s <file.vti> <prefix>\n", argv[0]);
    exit(1);
  }

  fileToRead = argv[1];
  prefix     = argv[2];

  FILE *f = fopen(fileToRead, "r");

  if (f == NULL) {
    printf("Unable to read file \"%s\"!\n", fileToRead);
    exit(1);
  }

  sprintf(outFilename, "%s_distfield", prefix);

  Phi *pf = (Phi *) malloc(sizeof(Phi));
  if (!pf) {
    printf("Error allocation memory for the phi function.\n");
    exit(1);
  }

  create_phi_function(f, pf);

  // print dimensions
  printf("Dimensions:\n");
  printf("x: %d\tdx: %f\n", pf->x, pf->dx);
  printf("y: %d\tdy: %f\n", pf->y, pf->dy);
  printf("z: %d\tdz: %f\n", pf->z, pf->dz);

  // print size
  printf("Size:\n");
  printf("x: %d\tdx: %f\n", pf->x, pf->dx);
  printf("y: %d\tdy: %f\n", pf->y, pf->dy);
  printf("z: %d\tdz: %f\n", pf->z, pf->dz);

  fclose(f);

  calc_dist_field(pf);

  phi_write_file(pf, VTI, outFilename);

  destroy_phi_function(pf);

  free(pf);

  return 0;
}
