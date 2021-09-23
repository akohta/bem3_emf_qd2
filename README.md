# bem3_emf_qd2
This is the three-dimensional electromagnetic field analysis program for two-dimensional periodic arrangement objects irradiated by a plane wave. 
This is based on boundary element method, the own developed numerical solution is used.
This is full vector field three-dimensional analysis, the corner problem free. 
Intel Math Kernel Library and libpng are required. 
Gmsh is used for create a mesh data of object. 
The calculation program of quasi-periodic Green's function "d3_qpgf_d2" is used.

![analysis model](model_qpbc2.png "analysis model (model_qpbc2.png)")  


## Usage of example code  

1. type 'make' command to compile.  
   The executable d3qd2_bv_solver, example1.out, example2.out, example3.out are created. 
   The d3qd2_bv_solver is the main solver of boundary integral equations. 
   The example1.out is the executable of source code example1.c, it shows a simplest example using "bem3_emf_qd2". 
   The example2.out is the executable of source code example2.c, it shows a example of electromagnetic field intensity analysis. 
   The example3.out is the executable of source code example3.c, it shows a example of outputting the instantaneous value of electromagnetic field as an image.  

2. type './d3qd2_bv_solver' with arguments of plane wave datafile name, periodicity datafile name, medium datafile name, mesh datafile name, output datafile name, rotation and translation settings (optional).   
   For example, './d3qd2_bv_solver ipw.txt periodicity_data.txt medium_data.txt sphere_m2.msh ex.dat'. 
   The ipw.txt is the sample of incident field datafile, a plane wave is defined in it. 
   The periodicity_data.txt is the sample of periodicity datafile, periodic boundary condition and lattice vectors are defined in it. 
   The medium_data.txt is the sample of medium datafile, two mediums are defined in it. The domain number is assinged to the medium from 1 in order. 
   The sphere_m2.msh is the sample of mesh datafile, it is a two-layered sphere object. 
   It was created by using Gmsh geometry file sphere_m2.geo in the mesh_sample folder.
   The sphere_m2_image.png is the visualization result of the sphere_m2.msh. 
   As a simple representation of the analysis model, the nodes used for the surface integral are output as point cloud data. 
   In this case, the file ex.particles is output and the visualization result is ex_particle.png (using ParaView).  

3. type './example1.out' with an argument of datafile name output by d3qd2_bv_solver.  
   For example, './example1.out ex.dat'. This executable calculates electromagnetic field, radiaton force and torque.  
   
4. type './example2.out' with an argument of datafile name output by d3qd2_bv_solver.  
   For example, './example2.out ex.dat'. This executable calculates electromagnetic field intensity distributions, outputs them to text files. 
   The I_example2.png is the visualization result of intensity distributions, created by using Gnuplot script gscritp_example2.plt
   (using ImageMagick to convert eps to png).  
   
5. type './example3.out' with an argument of datafile name output by d3b1_bv_solver.  
   For example, './example3.out ex.dat'. This executable calculates instantaneous value of the electromagnetic fields, outputs them to png image files. 
   The image files are output to the folder which has a name adding "images" to the datafile name specified in the argument (file-extension is excluded). 
   Each image file has a name that indicates the cross section, field component, and number of time steps (ex. xz_Ex_014.png). 
   The color bar is output as color_bar.png in the same folder. 
   The range of color bar in each cross section is output to the info.txt file (ex. xy_info.txt for z=0 plane). 
   The xz_Ex.gif, yz_Ex.gif and xy_Ex.gif are animated gifs that concatenate the png files created by using the shell script gif_animation.sh.  

Please see d3qd2_src/bem3_emf_qd2.h for detail of functions. 
The main parts of the code are parallelized by using OpenMP. 
The number of threads is controlled by the environment variable OMP_NUM_THREADS.  

![mash image 0](sphere_m2_image.png "mesh image of the object (sphere_m2_image.png)") 
![poing cloud data 0](ex_particles.png "nodes for surface integral (ex_particles.png)") 
![intensity distributions 0](I_example2.png "intensity distributions (I_example2.png)")  
![xz_Ex.gif](xz_Ex.gif "instantaneous value of the E_x on y=0 plane (xz_Ex.gif)")![yz_Ex.gif](yz_Ex.gif "instantaneous value of the E_x on x=0 plane (yz_Ex.gif)")  
![xy_Ex.gif](xy_Ex.gif "instantaneous value of the E_x on z=0 plane (xy_Ex.gif)")  


## Analysis sample 2 (in the analysis_sample2)  

This is the analysis result of plane wave scattering by cone objects. 

![mesh image 2](analysis_sample2/cone_m1_image.png "mesh image of the object (analysis_sample2/cone_m1_image.png)") 
![point cloud data 2](analysis_sample2/ex2_particles.png "nodes for surface integral (analysis_sample2/ex2_particles.png)") 
![intensity distributions 2](analysis_sample2/I_example2_logcb.png "intensity distributions (analysis_sample2/I_example2.png)") 
![xz_Ex.gif 2](analysis_sample2/xz_Ex.gif "instantaneous value of the E_x on y=0 plane (analysis_sample2/xz_Ex.gif)")![yz_Ex.gif 2](analysis_sample2/yz_Ex.gif "instantaneous value of the E_x on x=0 plane (analysis_sample2/yz_Ex.gif)")  
![xy_Ex.gif 2](analysis_sample2/xy_Ex.gif "instantaneous value of the E_x on z=0 plane (analysis_sample2/xy_Ex.gif)")  


## Verification  

The verification results are in the folder verification. 
The analysis result of symmetrically arranged 19 spheres using "emf_mie_ms" (in the folder emf_mie_ms_result) and the analysis result of periodic arrangement of spheres are shown. 

![point cloud data v](verification/emf_mie_ms_result/v1_particles.png "symmetrically arranged 19 spheres (verification/emf_mie_ms_result/v1_particles.png)") 


## About mesh file  

This code can use quadrangular ( bi-linear ) and triangular ( linear triangular ) elements. 
I recommend using quadrangular element for reduce required memory. 
The samples of mesh data are in the folder mesh_sample. 
The file with extension .geo is the Gmsh geometry file. 
The file with extension .msh is the mesh file created by Gmsh using the geometry file. 
These mesh files are created by the command 'gmsh -2 -tol 1.0e-15 xxxx.geo' in command line ( xxxx.geo is a geometry file). 
The domain number (Physical Surface) 99 is assigned to the open region in Gmsh geometry file, because Gmsh can't use the number 0 
(assigned to open region in the code). 
Please refer to the manual of Gmsh for detail of geometry file.  


## System of units  

This program use the own defined system of units (OSU), optimized for optics. 
The system of units is defined as <img src="https://latex.codecogs.com/gif.latex?c_0=1"> ( speed of light in vacuum ), 
<img src="https://latex.codecogs.com/gif.latex?\mu_0=1"> ( permeability of vacuum ). 
For the conversion from OSU to MKSA system of units, the unit of length in OSU is defined as 
<img src="https://latex.codecogs.com/gif.latex?1\times10^{-6}"> [m] in MKSA, the unit of power in OSU is defined as
<img src="https://latex.codecogs.com/gif.latex?1\times10^{-3}"> [W] in MKSA. The conversions of base unit are follows.  
<img src="https://latex.codecogs.com/gif.latex?a=1\times10^{-6}">,  
<img src="https://latex.codecogs.com/gif.latex?b=1\times10^{-3}">,  
<img src="https://latex.codecogs.com/gif.latex?a\,\mathrm{[m]}=1\,\mathrm{[L]}">,  
<img src="https://latex.codecogs.com/gif.latex?\frac{ab}{c_0^3}\,\mathrm{[kg]}=1\,\mathrm{[M]}">,  
<img src="https://latex.codecogs.com/gif.latex?\frac{a}{c_0}\,\mathrm{[s]}=1\,\mathrm{[T]}">,  
<img src="https://latex.codecogs.com/gif.latex?\sqrt{\frac{b}{c_0\mu_0}}\,\mathrm{[A]}=1\,\mathrm{[I]}">.  
Please see com_src/osu_mksa.h and com_src/osu_mksa.c for detail of conversions.  


## Reference

1. Intel Math Kernel Library [MKL](https://software.intel.com/mkl)  
2. The official PNG reference library [libpng](http://www.libpng.org/pub/png/libpng.html)  
3. Three-dimensional mesh generator [Gmsh](https://gmsh.info/)
4. The command-line driven graphing utility [gnuplot](http://www.gnuplot.info/)  
5. The utilities for manipulating images [ImageMagick](https://imagemagick.org/)  
6. The calculation program of quasi-periodic Green's function [d3_qpgf_d2](https://github.com/akohta/d3_qpgf_d2/)
7. The electromagnetic field analysis program [emf_mie_ms](https://github.com/akohta/emf_mie_ms/)  
8. The data analysis and visualization application [ParaView](https://www.paraview.org/)  
