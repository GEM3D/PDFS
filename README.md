#Parallel Distance Field Solver (PDFS)
A three dimensional Eikonal equation solver using the parallel Fast Sweeping Method (FSM) that computes the distance field of the given domain. [FSM](http://www.math.uci.edu/~zhao/homepage/research_files/FSM.pdf) is an iterative algorithm that uses upwind difference scheme for discretization and Gauss-Seidel iterations with alternating sweeping orderings to solve the discretized system.

Here, the two different versions of the parallel FSM are implemented.
* **[Shared Memory](http://www.sciencedirect.com/science/article/pii/S002199911200722X): Single GPU (CUDA)**

   ####Requirements:
   * NVIDIA Graphics Processing Unit (GPU)
   * GCC (GNU C Compiler)
   * NVCC (Nvidia C Compiler)
   * Make
   
* **[Hybrid Memory](README.md): MPI/OpenACC**

  ####Requirements:
  * Graphics Processing Unit (GPU)
  * OpenACC Compiler (PGCC)
  * MPI
  * NetCDF4 | HDF5 | SZIP
  * Make
  
###Code Compilation
Execute the make command from within the folder. Change the location of the directories of the libraries accordingly in the Makefile to correctly build the program.
    
    make
    
###Code Execution
Once you compile the code the binary executable is created within a folder called `bin` in the same directory. By default the executable is named as `PDFS`.
* **Single GPU CUDA**

        Usage: ./bin/PDFS <filename.vti> <outputPrefix>
        
        filename.vti: VTI input file
        outputPrefix: Prefix string to be added to the output file
* **MPI/OpenACC**

        Usage: ./bin/PDFS -i <filename.nc> -p <outputPrefix> [-x <val>] [-y <val>] [-z <val>]
        
        -i: input:  filename.nc:  NetCDF4 input file
        -p: prefix: outputPrefix: Prefix string to be added to the output file
        -x: Decomposition in x (Optional, default 1)
        -y: Decomposition in y (Optional, default 1)
        -z: Decomposition in z (Optional, default 1)
   
