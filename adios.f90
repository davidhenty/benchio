module adios

  use mpi
  use adios2
  implicit none


contains

subroutine adioswrite(filename, iodata, n1, n2, n3, cartcomm)

! ADIOS variables
  type(adios2_adios) :: adios2obj
  type(adios2_io) :: io
  type(adios2_engine) :: bp_writer
  type(adios2_variable) :: var_g 
  
  integer, parameter :: ndim = 3
  character*(*) :: filename
  
  integer :: n1, n2, n3
  double precision, dimension(0:n1+1,0:n2+1,0:n3+1) :: iodata
  double precision, dimension(n1,n2,n3) :: out_data

  integer*8, dimension(ndim) :: arraysize, arraystart
  integer*8, dimension(ndim) :: arraygsize, arraysubsize

  integer :: cartcomm, ierr, rank, size

  integer, dimension(ndim) :: dims, coords
  logical, dimension(ndim) :: periods


! initialise ADIOS using the MPI communicator and config file
  call adios2_init(adios2obj, "adios2.xml",cartcomm, ierr)
  call adios2_declare_io(io, adios2obj, 'Output', ierr )
        
  call MPI_Comm_size(cartcomm, size, ierr)
  call MPI_Comm_rank(cartcomm, rank, ierr)

  call MPI_Cart_get(cartcomm, ndim, dims, periods, coords, ierr)

  arraysize(:) = [n1+2, n2+2, n3+2]

! Subtract halos for array subsize

  arraysubsize(:)   = [n1, n2, n3]

! Define the global array size and the start coordinates of this ranks array
  arraygsize(:) = arraysubsize(:) * dims(:)
  arraystart(:) = arraysubsize(:) * coords(:)

! Get a copy of the array with the halos removed
  out_data = iodata(1:arraysubsize(1),1:arraysubsize(2),1:arraysubsize(3))

! Open the file
  call adios2_open (bp_writer, io, filename, adios2_mode_write, ierr)

! Define the global array
  call adios2_define_variable(var_g, io, "GlobalArray", adios2_type_dp, &
                              ndim, arraygsize, arraystart, arraysubsize, &
                              adios2_constant_dims, ierr)

! Begin ouput step
  call adios2_begin_step( bp_writer, ierr)

  call adios2_put( bp_writer, var_g, out_data, ierr)

! End the output
  call adios2_end_step(bp_writer, ierr)

! Close the file
  call adios2_close(bp_writer, ierr)

  call adios2_finalize(adios2obj, ierr)


end subroutine adioswrite

end module adios
