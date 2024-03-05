program benchio

  use benchclock
  use benchcomm
  use mpiio
  use ioserial
  use iohdf5
  use ionetcdf
  use adios
  use iolinear

  implicit none

  integer, parameter :: ndim = 3

  integer, parameter :: numiolayer = 8  
  integer, parameter :: numstriping = 3
  integer, parameter :: maxlen = 64

  character*(maxlen), dimension(numiolayer)  :: iostring, iolayername
  character*(maxlen), dimension(numstriping) :: stripestring
  character*(maxlen) :: argstring

  logical :: ioflag, stripeflag, globalflag
  logical, dimension(numiolayer)  :: doio
  logical, dimension(numstriping) :: dostripe

  character*(maxlen) :: filename, suffix

  integer :: iolayer, istriping, iolayermulti, iolayernode, numarg, iarg

! Set local array size - global sizes l1, l2 and l3 are scaled
! by number of processes in each dimension

  integer, parameter :: n1def = 128
  integer, parameter :: n2def = 128
  integer, parameter :: n3def = 128

  integer :: i1, i2, i3, j1, j2, j3, n1, n2, n3, l1, l2, l3, p1, p2, p3

  double precision, allocatable, dimension(:,:,:) :: iodata

  integer :: worldrank, worldsize, rank, size, ierr, comm, cartcomm, iocomm, dblesize
  integer, dimension(ndim) :: dims, coords

  integer, parameter :: iounit = 12
  integer, parameter :: kib = 1024
  integer, parameter :: mib = kib*kib
  integer, parameter :: gib = kib*mib

  logical :: reorder = .false.
  logical, dimension(ndim) :: periods = [.false., .false., .false.]

  double precision :: t0, t1, time, iorate, gibdata

  stripestring(1) = "unstriped"
  stripestring(2) = "striped"
  stripestring(3) = "fullstriped"

! These versions are special as they creates many files so need to record this

  iolayermulti = 2
  iolayernode  = 3
  
  iostring(1) = "Serial"
  iostring(2) = " Proc"
  iostring(3) = " Node "
  iostring(4) = "MPI-IO"
  iostring(5) = " HDF5 "
  iostring(6) = "NetCDF"
  iostring(7) = "Adios2"
  iostring(8) = "Linear"

  iolayername(1) = "serial"
  iolayername(2) = "proc"
  iolayername(3) = "node"
  iolayername(4) = "mpiio"
  iolayername(5) = "hdf5"
  iolayername(6) = "netcdf"
  iolayername(7) = "adios"
  iolayername(8) = "linear"

  call MPI_Init(ierr)

  comm = MPI_COMM_WORLD

  call MPI_Comm_size(comm, worldsize, ierr)
  call MPI_Comm_rank(comm, worldrank, ierr)

  ! These 3 options are block, cyclic1 and cyclic64

  
  !  call MPI_Comm_split(MPI_COMM_WORLD, (2*worldrank)/worldsize, worldrank, comm, ierr)
  !  call MPI_Comm_split(MPI_COMM_WORLD, mod(worldrank,2), worldrank, comm, ierr)
  call MPI_Comm_split(MPI_COMM_WORLD, (2*mod(worldrank,128))/128, worldrank, comm, ierr)
  

  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)

  ! These 3 options are block, cyclic1 and cyclic64

  !  if (worldrank < size) then
  !  if (mod(worldrank,2) == 0) then
  if ((2*mod(worldrank,128))/128 == 0) then

  ! Parse the arguments

  doio(:) = .false.
  dostripe(:) = .false.  
  
  numarg = command_argument_count()

  if (numarg < 4) then

     if (rank == 0) then
        write(*,*) "usage: benchio (n1, n2, n3) (local|global) [serial] [proc] [node]"
        write(*,*) "       [mpiio] [hdf5] [netcdf] [adios] [unstriped] [striped] [fullstriped]"
     end if

     call MPI_Finalize(ierr)
     stop
           
  end if

  call get_command_argument(1, argstring)
  read(argstring,*) n1
  call get_command_argument(2, argstring)
  read(argstring,*) n2
  call get_command_argument(3, argstring)
  read(argstring,*) n3

  globalflag = .true.

  call get_command_argument(4, argstring)
  if (argstring == "local") globalflag = .false.

  do iarg = 5, numarg

     ioflag = .false.
     stripeflag = .false.

     call get_command_argument(iarg, argstring)
     
     do iolayer = 1, numiolayer

        if (iolayername(iolayer) == argstring) then
           ioflag = .true.
           doio(iolayer) = .true.
        end if

     end do

     do istriping = 1, numstriping

        if (stripestring(istriping) == argstring) then
           stripeflag = .true.
           dostripe(istriping) = .true.
        end if

     end do

     if (.not.ioflag .and. .not.stripeflag) then
        
        write(*,*) "Illegal argument: ", argstring

        call MPI_Finalize(ierr)
        stop
           
     end if

  end do

  ! Check defaults

  if (count(doio(:)) == 0) then
     doio(:) = .true.
  end if

  if (count(dostripe(:)) == 0) then
     dostripe(:) = .true.
  end if

  ! Set 3D processor grid

  dims(:) = 0
  call MPI_Dims_create(size, ndim, dims, ierr)

! Reverse dimensions as MPI assumes C ordering (this is not essential)

  p1 = dims(3)
  p2 = dims(2)
  p3 = dims(1)

  ! Compute local sizes if required

  if (globalflag) then
     n1 = n1/p1
     n2 = n2/p2
     n3 = n3/p3
  end if

  ! Compute global sizes

  l1 = p1*n1
  l2 = p2*n2
  l3 = p3*n3

  call MPI_Type_size(MPI_DOUBLE_PRECISION, dblesize, ierr)

  gibdata = float(dblesize)*float(n1)*float(n2)*float(n3)*float(p1*p2*p3)
  gibdata = gibdata/float(gib)

  dims(1) = p1
  dims(2) = p2
  dims(3) = p3

  call MPI_Cart_create(comm, ndim, dims, periods, reorder, cartcomm, ierr)

  if (rank == 0) then
     write(*,*)
     write(*,*) "Simple Parallel IO benchmark"
     write(*,*) "----------------------------"
     write(*,*)
     write(*,*) "Running on ", size, " process(es)"
     write(*,*)
  end if

  ! Set up nodal stuff

  call initbenchnode(cartcomm)

  if (rank == 0) then

     write(*,*)
     write(*,*) "Process grid is (", p1, ", ", p2, ", ", p3, ")"
     write(*,*) "Array size is   (", n1, ", ", n2, ", ", n3, ")"
     write(*,*) "Global size is  (", l1, ", ", l2, ", ", l3, ")"
     write(*,*)
     write(*,*) "Total amount of data = ", gibdata, " GiB"
     write(*,*)
     write(*,*) "Clock resolution is ", benchtick()*1.0e6, ", usecs"
     write(*,*)
     write(*,*) "Using the following IO methods"
     write(*,*) "------------------------------"

     do iolayer = 1, numiolayer
        if (doio(iolayer)) write(*,*) iolayername(iolayer)
     end do

     write(*,*)
     write(*,*) "Using the following stripings"
     write(*,*) "-----------------------------"
     
     do istriping = 1, numstriping
        if (dostripe(istriping)) write(*,*) stripestring(istriping)
     end do

     write(*,*)

  end if

  ! Allocate data

  allocate(iodata(0:n1+1, 0:n2+1, 0:n3+1))

  ! Set halos to illegal values

  iodata(:,:,:) = -1
  
! Set iodata core to have unique values 1, 2, ..., p1*n1*p2*n2*p3*n3

  call MPI_Cart_coords(cartcomm, rank, ndim, coords, ierr)
  
  do i3 = 1, n3
     do i2 = 1, n2
        do i1 = 1, n1

           j1 = coords(1)*n1 + i1
           j2 = coords(2)*n2 + i2
           j3 = coords(3)*n3 + i3

           iodata(i1,i2,i3) = (j3-1)*l1*l2 + (j2-1)*l1 + j1

        end do
     end do
  end do

  do iolayer = 1, numiolayer

     if (.not. doio(iolayer)) cycle

     if (rank == 0) then
        write(*,*)
        write(*,*) "------"
        write(*,*) iostring(iolayer)
        write(*,*) "------"
        write(*,*)
     end if

!     if (iolayer == 3 .or. iolayer == 4) then
!
!        if (rank == 0) then
!           write(*,*) "WARNING: Skipping ", iostring(iolayer)
!        end if
!
!        cycle
!
!     end if

     do istriping = 1, numstriping

        if (.not. dostripe(istriping)) cycle

        filename = trim(stripestring(istriping))//"/"//trim(iolayername(iolayer))
        suffix = ""

        iocomm = cartcomm

        ! Deal with multiple files

        if (iolayer == iolayermulti) then
           iocomm = MPI_COMM_SELF
           write(suffix,fmt="(i6.6)") rank
        end if
           
        if (iolayer == iolayernode) then
           iocomm = nodecomm
           write(suffix,fmt="(i6.6)") nodenum
        end if
           
        suffix = trim(suffix)//".dat"
        filename = trim(filename)//suffix
        
        if (rank == 0) then
           write(*,*) "Writing to ", filename
        end if

        call MPI_Barrier(comm, ierr)
        t0 = benchtime()

        select case (iolayer)

        case(1:3)
           call serialwrite(filename, iodata, n1, n2, n3, iocomm)

        case(4)
           call mpiiowrite(filename, iodata, n1, n2, n3, iocomm)

        case(5)
           call hdf5write(filename, iodata, n1, n2, n3, iocomm)

        case(6)
           call netcdfwrite(filename, iodata, n1, n2, n3, iocomm)

        case(7)
           call adioswrite(filename, iodata, n1, n2, n3, iocomm)

        case(8)
           call linearwrite(filename, iodata, n1, n2, n3, iocomm)

        case default
           write(*,*) "Illegal value of iolayer = ", iolayer
           stop

        end select

        call MPI_Barrier(comm, ierr)
        t1 = benchtime()

        time = t1 - t0
        iorate = gibdata/time

        if (rank == 0) then
           write(*,*) "time = ", time, ", rate = ", iorate, " GiB/s"
        end if

        ! Rank 0 in iocomm deletes
        if (iolayer == 7) then
          ! ADIOS makes a directory so the file deletion function will not work
          ! use the shell instead

          call MPI_Barrier(comm, ierr)
          if (rank == 0) then
            call execute_command_line("rm -r "//filename)
          end if 
          call MPI_Barrier(comm, ierr)
        
        else
           call bossdelete(filename, iocomm)
        endif
        
     end do
  end do

  if (rank == 0) then
     write(*,*)
     write(*,*) "--------"
     write(*,*) "Finished"
     write(*,*) "--------"
     write(*,*)
  end if

  end if

  call MPI_Finalize(ierr)
  
end program benchio
