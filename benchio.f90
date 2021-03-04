program benchio

  use benchclock
  use benchcomm
  use mpiio
  use ioserial
  use iohdf5
  use ionetcdf

  implicit none

  integer, parameter :: numiolayer = 6
  integer, parameter :: numstriping = 3
  integer, parameter :: maxlen = 64

  character*(maxlen), dimension(numiolayer)  :: iostring, iolayername
  character*(maxlen), dimension(numstriping) :: stripestring

  character*(maxlen) :: filename, suffix

  integer :: iolayer, istriping, iolayermulti, iolayernode

! Set local array size - global sizes l1, l2 and l3 are scaled
! by number of processes in each dimension

  integer, parameter :: n1 = 128
  integer, parameter :: n2 = 128
  integer, parameter :: n3 = 128

  integer :: i1, i2, i3, j1, j2, j3, l1, l2, l3, p1, p2, p3

  double precision :: iodata(0:n1+1, 0:n2+1, 0:n3+1)

  integer :: rank, size, ierr, comm, cartcomm, iocomm, dblesize
  integer, dimension(ndim) :: dims, coords

  integer, parameter :: iounit = 12
  integer, parameter :: kib = 1024
  integer, parameter :: mib = kib*kib
  integer, parameter :: gib = kib*mib

  logical :: reorder = .false.
  logical, dimension(ndim) :: periods = [.false., .false., .false.]

  double precision :: t0, t1, time, iorate, gibdata

  stripestring(1) = 'unstriped'
  stripestring(2) = 'striped'
  stripestring(3) = 'fullstriped'

! These versions are special as they creates many files so need to record this

  iolayermulti = 2
  iolayernode  = 3

  iostring(1) = 'Serial'
  iostring(2) = ' Multi'
  iostring(3) = ' Node '
  iostring(4) = 'MPI-IO'
  iostring(5) = ' HDF5 '
  iostring(6) = 'NetCDF'

  iolayername(1) = 'serial'
  iolayername(2) = 'rank'
  iolayername(3) = 'node'
  iolayername(4) = 'mpiio'
  iolayername(5) = 'hdf5'
  iolayername(6) = 'netcdf'

  call MPI_Init(ierr)

  comm = MPI_COMM_WORLD

  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)

  dims(:) = 0

! Set 3D processor grid

  call MPI_Dims_create(size, ndim, dims, ierr)

! Reverse dimensions as MPI assumes C ordering (this is not essential)

  p1 = dims(3)
  p2 = dims(2)
  p3 = dims(1)

! Compute global sizes

  l1 = p1*n1
  l2 = p2*n2
  l3 = p3*n3

  call MPI_Type_size(MPI_DOUBLE_PRECISION, dblesize, ierr)

  gibdata = float(dblesize*n1*n2*n3)*float(p1*p2*p3)/float(gib)

  dims(1) = p1
  dims(2) = p2
  dims(3) = p3

  call MPI_Cart_create(comm, ndim, dims, periods, reorder, cartcomm, ierr)

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Simple Parallel IO benchmark'
     write(*,*) '----------------------------'
     write(*,*)
     write(*,*) 'Running on ', size, ' process(es)'
     write(*,*)
  end if

  ! Set up nodal stuff

  call initbenchnode(cartcomm)

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Process grid is (', p1, ', ', p2, ', ', p3, ')'
     write(*,*) 'Array size is   (', n1, ', ', n2, ', ', n3, ')'
     write(*,*) 'Global size is  (', l1, ', ', l2, ', ', l3, ')'
     write(*,*)
     write(*,*) 'Total amount of data = ', gibdata, ' GiB'
     write(*,*)
     write(*,*) 'Clock resolution is ', benchtick()*1.0e6, ', usecs'
  end if
  
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

     if (rank == 0) then
        write(*,*)
        write(*,*) '------'
        write(*,*) iostring(iolayer)
        write(*,*) '------'
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

        filename = trim(stripestring(istriping))//'/'//trim(iolayername(iolayer))
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
           write(*,*) 'Writing to ', filename
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

        case default
           write(*,*) 'Illegal value of iolayer = ', iolayer
           stop

        end select

        call MPI_Barrier(comm, ierr)
        t1 = benchtime()

        time = t1 - t0
        iorate = gibdata/time

        if (rank == 0) then
           write(*,*) 'time = ', time, ', rate = ', iorate, ' GiB/s'
        end if

        ! Rank 0 in iocomm deletes

        call bossdelete(filename, iocomm)
        
     end do
  end do

  if (rank == 0) then
     write(*,*)
     write(*,*) '--------'
     write(*,*) 'Finished'
     write(*,*) '--------'
     write(*,*)
  end if

  call MPI_Finalize(ierr)
  
end program benchio
