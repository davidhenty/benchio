!
! The "serial" IO routine takes a communicator argument. This enables
! it to be used for a variety of purposes:
!
! MPI_COMM_WORLD: standard "master" IO from a single process
! MPI_COMM_NODE:  file-per-node
! MPI_COMM_SELF:  file-per-process
!

module ioserial

  contains

subroutine serialwrite(filename, iodata, n1, n2, n3, comm)

  character*(*) :: filename
  
  integer :: n1, n2, n3
  double precision, dimension(0:n1+1,0:n2+1,0:n3+1) :: iodata

  integer :: comm, ierr, rank, size
  integer, parameter :: iounit = 10

  integer :: i

  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)

!  Write same amount of data as the parallel write but do it all from rank 0
!  This is just to get a baseline figure for serial IO performance - note
!  that the contents of the file will be different from the parallel calls

  if (rank == 0) then

     open(file=filename, unit=iounit, access='stream')

     do i = 1, size
        write(unit=iounit) iodata(1:n1, 1:n2, 1:n3)
     end do

     close(iounit)

  end if

end subroutine serialwrite

end module ioserial
