# bem3_emf_qd2
This is the three-dimensional electromagnetic field analysis program for two-dimensional periodic arrangement objects irradiated by a plane wave. 
This is based on boundary element method, the own developed numerical solution is used.
This is full vector field three-dimensional analysis, the corner problem free. 
Intel Math Kernel Library is required. 
Gmsh is used for create a mesh data of object. 
The calculation program of quasi-periodic Green's function "d3_qpgf_d2" is used.


## Usage of example code
1. type 'make' command to compile  
   The executable d3qd2_bv_solver, example1.out, example2.out are created. 
   The d3qd2_bv_solver is the main solver of boundary integral equations. 
   The example1.out is the executable of souce code example1.c, it shows a simplest example of usage. 
   The example2.out is the executable of souce code example2.c, it shows a example of electromagnetic field intensity analysis.  

2. type './d3qd2_bv_solver' with arguments of plane wave datafile name, periodicity datafile name, medium datafile name, mesh datafile name, output datafile name,
   rotation and translation settings ( optional ).   
   For example, './d3qd2_bv_solver ipw.txt periodicity_data.txt medium_data.txt sphere_m2.msh ex.dat'. 
   The ipw.txt is the sample of incident field datafile, a plane wave is defined in it. 
   The periodicity_data.txt is the sample of periodicity datafile, periodic boundary condition and lattice vectors are defined in it. 
   The medium_data.txt is the sample of medium datafile, two mediums are defined in it. The domain number is assinged to the medium from 1 in order. 
   The sphere_m2.msh is the sample of mesh datafile, it is a two layered sphere object. It was created by Gmsh geometry file sphere_m2.geo in mesh_sample folder.
   The sphere_m2_image.pdf is the visualization result of the sphere_m2.msh.

3. type './example1.out' with a argument of datafile outputed by d3qd2_bv_solver.  
   For example, './example1.out ex.dat'. This executable calculates electromagnetic field, radiaton force and torque.  
   
4. type './example2.out' with a argument of datafile outputed by d3qd2_bv_solver.  
   For example, './example2.out ex.dat'. This executable calculates electromagnetic field intensity distributions, outputs them to text files. 
   The I_example2.pdf is the visualization result of intensity distributions, created by the Gnuplot script gscritp_example2.plt.

The folder analysis_sample1, analysis_sample2, analysis_sample3 are the analysis results using these executable. 
Please see 'd3qd2_src/bem3_emf_qd2.h' for detail of functions. The main parts of the code are parallelized by using OpenMP. The number of threads is controlled by the environment variable 'OMP_NUM_THREADS'.  

## About mesh file 
This code can use quadrangular ( bi-linear ) and triangular ( linear triangular ) elements. I recommend using quadrangular element for reduce required memory. The samples of mesh data is in the folder mesh_sample. The file with extension '.geo' is the Gmsh geometry file. The file with extension '.msh' is the mesh file created by Gmsh using the geometry file. These mesh files are created by the command 'gmsh -2 -tol 1.0e-15 xxxx.geo' in command line ( xxxx.geo is a geometry file). The domain number ( Physical Surface ) 99 is assigned to the open region in Gmsh geometry file, becase Gmsh can't use the number 0 ( assigned to open region in the code). Please refer to the manual of Gmsh for detail of geometry file.  

## Verification  
The verification results are in the folder verification. The analysis result of symmetrically arranged 19 spheres ( in the folder emf_mie_ms_result ) and the analysis result of periodic arrangement of spheres are shown. In this case, the similar results are shown.  

## Reference
1. Intel Math Kernel Library [MKL](https://software.intel.com/mkl)
2. Three-dimensional mesh generator [Gmsh](https://gmsh.info/)
3. The calculation program of quasi-periodic Green's function [d3_qpgf_d2](https://github.com/akohta/d3_qpgf_d2/)
4. The electromagnetic field analysis program [emf_mie_ms](https://github.com/akohta/emf_mie_ms/)
