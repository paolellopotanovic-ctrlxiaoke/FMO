# Building Netlib BLAS and LAPACK manually without CMake

These instructions assume that you are in the GAMESS directory.

1.  Ensure that the Netlib BLAS and LAPACK libraries are available:

	`tools/lapack/download-lapack.csh`

2.  Navigate to the source directory:

	`cd 3rd-party/lapack`

3.  Copy over the appropriate `make.inc` file based on your Fortran compiler choice:

	```
	# if using gfortran
	cp ../../tools/lapack/make.inc.gfortran.64 make.inc

	# if using ifort
	cp ../../tools/lapack/make.inc.ifort.64 make.inc

	# if using nvfortran
	cp ../../tools/lapack/make.inc.nvfortran.64 make.inc
	```

4.  Build BLAS:

	`make -j blaslib`

5.  Build LAPACK:

	`make -j lapacklib`

6.  Navigate one directory level up and create a `lib` folder:

	```
	cd ../
	mkdir lib
	```

7.  Copy and rename the libraries over to the `lib` folder:

	```
	cp lapack/liblapack.a lib/liblapack64.a
	cp lapack/librefblas.a lib/libblas64.a
	```

8.  Navigate to back to the GAMESS directory and complete your compilation of GAMESS:

	If you are using compall:

	```
	cd ../
	./compall
	./lked
	```

	If you are using Make:

	```
	cd ../
	make ddi
	make modules
	make -j gamess
	./lked
	```
