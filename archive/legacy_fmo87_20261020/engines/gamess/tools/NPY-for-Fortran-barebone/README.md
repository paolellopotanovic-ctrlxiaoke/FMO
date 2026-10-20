# NPY for Fortran (barebone)
Note: This is a fork that keeps only the bare minimum of functionality.
Compared to the original "NPY for Fortran" the following functionality is removed:
* data types other than integer(8) or real(8)
* dimensions other than 1 (vector) or 2 (matrix)
* ability to run on big endian machines
* generation of .npz files

The documentation given below is retained from the original project.

---

This Fortran module allows to save numerical Fortran arrays in Numpy's .npy format.
Currently supported are:
```fortran
1. integer(8)
2. real(8)
```

### *.npy files
Saving an array into a .npy-file is simply done by calling:
```fortran
call save_npy("filename.npy", array)
```
