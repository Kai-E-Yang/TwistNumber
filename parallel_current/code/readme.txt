The program main.f90 is the main program.
This program is used for calculating the squashing factor Q.
The method is proposed by Pariat & Demoulin 2012.
This version is usually used by myself.
Version v2.0
History:
	K.Y. 2016@NJU
	K.Y. Jun/2017@MSU change the compiler from ifort to gfortran.
	K.Y. Mar/2017@MSU change the reading file from the initial formatted to unformatted.
	K.Y. Jun/2018@NJU change the most of names of the subroutines to make it look better.
	K.Y. Oct/2018@NJU use Scott's method in ApJ 2017
	K.Y. Nov/2018@NJU add the namelist format for reading the parameters.
   
The steps for this program is:

1	create the magnetic field data with the program output_b.pro

2	> make -f makefile

3	copy the executable file "twist" and the parameter file "par" to the document 'data' where the magnetic field data locates.

4 	change the input parameter in "par" by hand to select the subregion, the level of grid refinement, the name of the output file, or the integral step. The default level of grid refinement is 4, the default computational region is the whole region, the default value of the integral step is 0.25.

-----------------------------------------
the format of the par file is namelist in fortran, e.g.,

&filename
	BfieldName = '*****.binary'  ! file name for the Bxyz components 
	OutFileName = '*****.binary' ! file name for the out put file
	CurrentFile = '*****.binary' ! file name for the current data, 
	                             ! if you do not give a current, the program 
	                             ! will calculate the current field automatically.
	MaskFile = '*****.binary'    ! file name for the mask data, if you do not use mask, you not have to define it
/
&cal_par
	nthreads = 1                 ! thread number used for parallel computation
	dimx = 50                    ! dimension for x
	dimy = 50                    ! dimension for y
	dimz = 50                    ! dimension for z
	x_start = 1                  ! start point in x (it starts from 1 to dim X)
	y_start = 1                  ! start point in y (it starts from 1 to dim Y)
	z_start = 1                  ! start point in z (it starts from 1 to dim Z)
	x_end = 50                   ! end point in x (it starts from 1 to dim X)
	y_end = 50                   ! end point in y (it starts from 1 to dim Y)
	z_end = 1                    ! end point in z (it starts from 1 to dim Z)
	nlevel = 10                  ! grid level
	delta_s = 0.25               ! integral step (the default value is 0.25)
	usemask = .true.             ! logical variable that define you use mask data or not, 
	                             ! default value is .false., if you do not use mask data,
	                             ! you do not have to define it.
	mask_level = 4               ! the refine level of the mask data, since the mask
	                             ! usually has a dense mesh than the B field data, 
	                             ! so that you have to give the correct value of mask_level,
	                             ! if you give the mask data.
	usecurrent = .true.          ! set this variable as .true. means that you want to use
	                             ! the current data contained in CurrentFile rather than the !
	                             ! current calculated by this program. The default is .false.
/
-----------------------------------------

5	run the program by ./twist par (par is the parameter file)

6	use the idl program read_qsl.pro to read the result.
    IDL> .com read_twist.pro
    IDL> result = read_twist(par_name); par_name is a string contains the name of the par file, e.g., 'par'

7   in the new version, read_twist is removed, you can just define an array with the same dimension 
    you find in the screen output.

8   the test dir contains a simple field that you can try, I use a simply current file that J = 0.5*B,
    so that the twist number you get should be same as field line length, except a constant 0.5.

In addition:

This parallel version is based on FORTRAN OPENMP.

If one want to use N threads for the calculation, just change the value of parameter 'nthreads' in parameter file. If this parameter is defined as 0, then the max number of threads in the computer will be used as default.

The field line integral method is Runge-Kutta 4(5).