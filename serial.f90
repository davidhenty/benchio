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
!  that the contents of the file will be differnent from the parallel calls

  if (rank == 0) then

     open(file=filename, unit=iounit, access='stream')

     do i = 1, size
        write(unit=iounit) iodata(1:n1, 1:n2, 1:n3)
     end do

     close(iounit)

  end if

end subroutine serialwrite

subroutine multiwrite(filename, iodata, n1, n2, n3, comm)

  character*(*) :: filename
  
  integer :: n1, n2, n3
  double precision, dimension(0:n1+1,0:n2+1,0:n3+1) :: iodata

  integer :: comm
  integer, parameter :: iounit = 10

  open(file=filename, unit=iounit, access='stream')

  write(unit=iounit) iodata(1:n1, 1:n2, 1:n3)

  close(iounit)

end subroutine multiwrite

subroutine fdelete(filename)

  implicit none

  character *(*) :: filename
  integer, parameter :: iounit = 15
  integer :: stat

!  write(*,*) 'Deleting: ', filename

  open(unit=iounit, iostat=stat, file=filename, status='old')
  if (stat.eq.0) close(unit=iounit, status='delete')

end subroutine fdelete

end module ioserial
