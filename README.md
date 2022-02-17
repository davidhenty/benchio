 # benchio
Simple Fortran parallel IO benchmark for teaching and benchmarking purposes.

Benchio builds on a one-dimensional parallel IO benchmark previously
developed for the EU-funded EUFORIA project. See "High Performance
I/O", Adrian Jackson, Fiona Reid, Joachim Hein, Alejandro Soba and
Xavier Saez;
[https://ieeexplore.ieee.org/document/5739034/]([https://ieeexplore.ieee.org/document/5739034/).

Note that, before running the benchmark, you *must* set the Lustre striping on the three directories `unstriped`, `striped` and `fullstriped`.

 * Set `unstriped` to have a single stripe: `lfs setstripe -c 1 unstriped`
 * Set `fullstriped` to use the maximum number of stripes: `lfs setstripe -c -1 fullstriped`
 * Set `striped` to use an intermediate number of stripes, e.g. for 4 stripes: `lfs setstripe -c 4 striped`

The program has a very basic set of command-line options. The first
three arguments must be the dimensions of the dataset; the fourth
argument specifies if these are local sizes (i.e. weak scaling), or
global sizes (strong scaling).

For example, to run using a 256 x 256 x 256 data array on every
process (i.e. weak scaling):
````
benchio 256 256 256 local
````
In this case, the total file size will scale with the number of
processes. If run on 8 processes then the total file size would be 1
GiB.

To run using a 256 x 256 x 256 global array (i.e. strong scaling):
````
benchio 256 256 256 global
````
In this case, the file size will be 128 MiB regardless of the number
of processes.

If the local array size is n1 x n2 x n3, then the double precision
arrays are defined with halos as: `double precision :: iodata(0:n1+1,
0:n2+1, 0:n3+1)`.

A 3D cartesian topology p1 x p2 x p3 is created with dimensions
suggested by `MPI_Dims_create()` to create a global 3D array of size
l1 x l2 x l3 where l1 = p1 x n1 etc.
 
 The entries of the distributed IO array are set to globally unique
 values 1, 2, ... l1xl2xl3 using the normal Fortran ordering; the halo
 values are set to -1. When writing to file, the halos are omitted.
 
  
The code can use six IO methods, and for each of them can use up to
three directories with different stripings.  The IO methods are:
 
 1. Serial IO from one controller process to a single file `serial.dat` using Fortran binary unformatted `write` with `access = stream`
 2. Multiple serial IO (file-per-process) to *P* files `rankXXXXXX.dat` using Fortran binary unformatted `write` with `access = stream`
 3. Multiple serial IO (file-per-node) to *Nnode* files `nodeXXXXXX.dat` using Fortran binary unformatted `write` with `access = stream`
 4. MPI-IO collective IO to a single file `mpiio.dat` using native (i.e. binary) format
 5. HDF5 collective IO to a single file `hdf5.dat`
 6. NetCDF collective IO single file `netcdf.dat`
 
 Note that the serial part is designed to give a baseline IO rate. For simplicity, and to ensure we write the same amount of data as for the parallel
 methods, rank 0 writes out its
 own local array `size` times in succession. Unlike the parallel IO formats, the contents of the file will therefore *not* be a linearly increasing set of
 values 1, 2, 3, ..., l1xl2xl3.

If only the first four mandatory arguments are specified then all six
IO methods and all three stripings are used. However, you can pick
subsets by setting additional optional command-line options.

The full
set of options is:
````
benchio (n1, n2, n3) (local|global)
        [serial] [proc] [node] [mpiio] [hdf5] [netcdf]
	[unstriped] [striped] [fullstriped]
````