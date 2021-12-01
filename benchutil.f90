module benchcomm

  ! Various MPI-related support routines

  use mpi

  implicit none

  integer :: nodecomm, nodebosscomm, nodenum

contains

  subroutine initbenchnode(comm)

    integer :: comm

    integer :: size, rank, nodesize, noderank, spansize
    integer :: ierr, tag, colour, key, irank, namelen

    integer, dimension(MPI_STATUS_SIZE) :: status
    character*(MPI_MAX_PROCESSOR_NAME) :: nodename

    call MPI_Comm_size(comm, size, ierr)
    call MPI_Comm_rank(comm, rank, ierr)

    ! Create node-local communicators

    call MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, &
                             MPI_INFO_NULL, nodecomm, ierr)

    call MPI_Comm_size(nodecomm, nodesize, ierr)
    call MPI_Comm_rank(nodecomm, noderank, ierr)

! Create spanning communicator for all the node bosses
! Put everyone else into the same (junk) comm

    colour = min(noderank, 1)
    key    = noderank

    call MPI_Comm_split(comm, colour, key, nodebosscomm, ierr)

    call MPI_Comm_size(nodebosscomm, spansize, ierr)
    call MPI_Comm_rank(nodebosscomm, nodenum,  ierr)

! Make sure all ranks on node know the node number

    call MPI_Bcast(nodenum, 1, MPI_INTEGER, 0, nodecomm, ierr)

    call MPI_Get_processor_name(nodename, namelen, ierr)

! Get the stats

    if (rank == 0) then

       if (spansize /= 1) then
          write(*,*) "Running on ", spansize, " nodes"
       else
          write(*,*) "Running on ", spansize, " node"
       end if

       do irank = 0, spansize-1

          tag = nodesize

          if (irank /= 0) then

             call MPI_Recv(nodename, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                  irank, MPI_ANY_TAG, nodebosscomm, status, ierr)

             tag = status(MPI_TAG)
             call MPI_Get_count(status, MPI_CHARACTER, namelen, ierr)

          end if

          if (tag /= 1) then

             write(*,*) "Node number", irank, " is ", nodename(1:namelen), &
                  " with ", tag, " processes"
          else

             write(*,*) "Node number", irank, " is ", nodename(1:namelen), &
                  " with ", tag, " process"
          end if

       end do

    else if (noderank == 0) then

       ! Send number of processes as tag

       tag = nodesize

       call MPI_Ssend(nodename, namelen, MPI_CHARACTER, 0, tag, nodebosscomm, ierr)

    end if

  end subroutine initbenchnode

  subroutine bossdelete(filename, comm)

    use mpi

    implicit none

    character *(*) :: filename
    integer :: comm
    
    integer, parameter :: iounit = 15
    integer :: rank, ierr, stat

    call MPI_Comm_rank(comm, rank, ierr)

    if (rank == 0) then

       open(unit=iounit, iostat=stat, file=filename, status='old')
       if (stat.eq.0) close(unit=iounit, status='delete')

    end if

  end subroutine bossdelete

end module benchcomm


module benchclock

  implicit none

  logical,          save, private :: firstcall = .true.
  double precision, save, private :: ticktime = 0.0

  integer, parameter :: int32kind = selected_int_kind( 9)
  integer, parameter :: int64kind = selected_int_kind(18)

!
!  Select high resolution clock
!

  integer, parameter :: intkind = int64kind
  integer(kind = intkind) :: clkcount, clkrate

contains

double precision function benchtime()

  double precision :: dummy

! Ensure clock is initialised  

  if (firstcall) dummy = benchtick()

  call system_clock(clkcount)

  benchtime  = dble(clkcount)*ticktime

end function benchtime


double precision function benchtick()

  if (firstcall) then

     firstcall = .false.
     call system_clock(clkcount, clkrate)
     ticktime = 1.0d0/dble(clkrate)

  end if

  benchtick = ticktime

end function benchtick

end module benchclock
