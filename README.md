# benchio
Simple Fortran parallel IO benchmark for teaching and benchmarking purposes

Note that, before running the benchmark, you *must* set the LFS striping on the three directories `unstriped`, `striped` and `fullstriped`.

 * Set `unstriped` to have a single stripe: `lfs setstripe -c 1 unstriped`
 * Set `fullstriped` to use the maximum number of stripes: `lfs setstripe -c -1 fullstriped`
 * Set `striped` to use an intermediate number of stripes, e.g. for 4 stripes: `lfs setstripe -c 4 striped`

 The benchmark uses weak scaling, i.e. the local 3D array size is defined
 in the code and simply replicated across the paralle processes. If the local array size is n1 x n2 x n3, then the double precision
 arrays are defined with halos as: `double precision :: iodata(0:n1+1, 0:n2+1, 0:n3+1)`. A 3D cartesian topology p1 x p2 x p3 is created with dimensions suggested
 by `MPI_Dims_create()` to create a global 3D array of size l1 x l2 x l3 where l1 = p1 x n1 etc.
 
 The entries of the distributed IO array are set to a globally unique values 1, 2, ... l1xl2xl3 using the normal Fortran ordering; the
 halo values are set to -1. When writing
 to file, the halos are omitted.
 
  
 The code uses four IO methods, and for each of them writes to the three directories in turn.  The IO methods are:
 
 1. Serial IO using Fortran binary unformatted `write` with `access = stream`
 2. MPI-IO using native (i.e. binary) format and collective IO
 3. HDF5 collective IO
 4. NetCDF collective IO
 
 Note that the serial is designed to give a baseline IO figure. For simplicity, and to ensure it writes the same amount of data, rank 0 writes out its
 own local array `size` times in succession. Unlike the parallel IO formats, the contents of the file will *not* be a linearly increasing set of
 values 1, 2, 3, ..., l1xl2xl3.
