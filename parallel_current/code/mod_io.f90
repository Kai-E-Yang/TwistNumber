module mod_io
use mod_param
use mod_operator
implicit none
contains
  subroutine read_parameter()
    implicit none
    ! set default value    
    BfieldName   = 'bfield.binary'
    JfieldName   = 'none'
    OutFileName  = 'result.binary'
    MaskFile     = 'none'
    indataformat = 'binary'

    nthreads    = 1
    dimx        = 100
    dimy        = 100
    dimz        = 100
    
    x_start     = 100
    y_start     = 100
    z_start     = 100
    
    x_end       = 100
    y_end       = 100
    z_end       = 100
    
    nlevel      = 0
    mask_level  = 0
    delta_s     = 0.1d0
    usemask     = .false.
    usecurrent  = .false.
    ! read parameters   
    write(*,'(A)')'| = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ='   
    write(*,'(A)')'| Reading parameters from file: '//trim(par)
    open(Unit=1,file=par,STATUS='OLD',ACTION='READ')
    read(1,nml=filename_par,end=3)
3   rewind(1)
    read(1,nml=cal_par,end=5)
5   close(1)

  end subroutine read_parameter

  subroutine allocate_var()
    allocate(Bx(0:dimx+1,0:dimy+1,0:dimz+1))
    allocate(By(0:dimx+1,0:dimy+1,0:dimz+1))
    allocate(Bz(0:dimx+1,0:dimy+1,0:dimz+1))

    allocate(Jx(0:dimx+1,0:dimy+1,0:dimz+1))
    allocate(Jy(0:dimx+1,0:dimy+1,0:dimz+1))
    allocate(Jz(0:dimx+1,0:dimy+1,0:dimz+1))

    refine_dimx = nlevel*(x_end-x_start) + 1
    refine_dimy = nlevel*(y_end-y_start) + 1
    refine_dimz = nlevel*(z_end-z_start) + 1

    allocate(cal_data(refine_dimx,refine_dimy,refine_dimz,3))

    BoundaryEps = real(1,kind=r8)+1.0e-2
    iterMax = 20
    if(usemask) then
        allocate(mask(0: mask_level*(dimx-1) + 2, 0:mask_level*(dimy-1) + 2))
    end if
  end subroutine allocate_var

  subroutine deallocate_var()
    implicit none  
    deallocate(bx,by,bz)
    deallocate(jx,jy,jz)
    deallocate(cal_data)
    if(indataformat .eq. 'blk') then
      deallocate(x,y,z)
    end if
  end subroutine deallocate_var

  subroutine show_information()
    write(*,'(A)')'|--------------------------------------------------------'
    write(*,'(A)')'| Mission: Calculate the twist number using parallel current integral.'
    write(*,'(A)')'| 3D B field data comes from: '
    write(*,'(A)')'|                    '//trim(BfieldName)
    write(*,'(A)')'| The result is saved in file: '
    write(*,'(A)')'|                    '//trim(OutFileName)
    write(*,'(A)')'| The array of the result contains three value'
    write(*,'(A)')'| The first dimension is the sign(Bz)*log10(Q)'
    write(*,'(A)')'| The second dimension is length of the field lines'
    write(*,'(A)')'| The third dimension is mark of the boundary'
    write(*,'(A)')'| The meaning of this value is:'
    write(*,'(A)')'| '
    write(*,'(A)')'| 1: close field line with both ends root in bottom boundary'
    write(*,'(A)')'| 2: open field line, positive end roots in bottom'
    write(*,'(A)')'| 3: open field line, negative end roots in bottom'
    write(*,'(A)')'| 4: field line without end roots in the bottom boundary'
    write(*,'(A)')'| '
    write(*,'(A)')'|--------------------------------------------------------'
    write(*,'(A)')'| The NDIM of the Magnetic data is: '
    write(*,'(A,I3)')'|                  x dimension: ',dimx
    write(*,'(A,I3)')'|                  y dimension: ',dimy
    write(*,'(A,I3)')'|                  z dimension: ',dimz
    write(*,'(A)')'| '
    write(*,'(A)')'| The calculated domain is: '
    write(*,'(A,I3,A,I3)')'|                  X ',x_start,' to ',x_end
    write(*,'(A,I3,A,I3)')'|                  Y ',y_start,' to ',y_end
    write(*,'(A,I3,A,I3)')'|                  Z ',z_start,' to ',z_end
    write(*,'(A,I4)')'| The refine level is: ',nlevel
    write(*,'(A)')'| The NDIM of the result is: '
    write(*,'(A,I8)')'|                  x dimension: ',refine_dimx
    write(*,'(A,I8)')'|                  y dimension: ',refine_dimy
    write(*,'(A,I8)')'|                  z dimension: ',refine_dimz
    write(*,'(A)')'| The data dimension is: '
    write(*,'(A,I4,A,I4,A,I4,A,I2)')'|           ',refine_dimx,' X ',refine_dimy,' X ',refine_dimz,' X ',3
    write(*,'(A,I4,A)')'| This mission use ',nthreads,' threads'
    write(*,'(A)')'|--------------------------------------------------------'
  end subroutine show_information

  subroutine read_blk()
    logical :: alive
    integer(kind=i4) :: i,j,k
    integer(kind=i8) :: ntotal
    character*1024   :: header
    double precision :: ttime

    inquire(file=BfieldName,exist=alive)
    if(alive) then 
      write(*,*)'OUTPUT_B: loading data and preparing input files for QSL3D'
      OPEN(UNIT=3,STATUS='OLD',ACTION='READ',FILE=BfieldName,POSITION='REWIND',FORM='UNFORMATTED')
      read(3) header
      read(3) ntotal,dimx,dimy,dimz
      read(3) ttime
      allocate(x(0:dimx+1,0:dimy+1,0:dimz+1))
      allocate(y(0:dimx+1,0:dimy+1,0:dimz+1))
      allocate(z(0:dimx+1,0:dimy+1,0:dimz+1))
      do k=1,dimz
      do j=1,dimy
      do i=1,dimx
         read(3) x(i,j,k),y(i,j,k),z(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
      end do
      end do
      end do 
      close(3)

!$omp parallel do
      do k=1,dimz
      do j=1,dimy
      do i=1,dimx
         Jx(i,j,k) = yderiv(Bz,i,j,k) - zderiv(By,i,j,k) 
         Jy(i,j,k) = zderiv(Bx,i,j,k) - xderiv(Bz,i,j,k)
         Jz(i,j,k) = xderiv(By,i,j,k) - yderiv(Bx,i,j,k)
      enddo
      enddo
      enddo
!$omp end parallel do

    else
      write(*,'(A)')'| File '//BfieldName//' does not exist'
      stop
    end if
  end subroutine read_blk

  subroutine read_binary()
    logical          :: alive
    integer(kind=i4) :: i,j,k

    inquire(file=BfieldName,exist=alive)
    if(alive) then 
      write(*,'(A)')'| Loading magnetic field data...'
      open(3,File=BfieldName,Access="stream",Form = "unformatted" )
      read(3) Bx(1:dimx,1:dimy,1:dimz)
      read(3) By(1:dimx,1:dimy,1:dimz)
      read(3) Bz(1:dimx,1:dimy,1:dimz)
      close(3)
      write(*,'(A)')'| Data sucessfully loaded !'

      if(usecurrent) then
        open(3,File=JfieldName,Access="stream",Form = "unformatted" )
        read(3) Jx(1:dimx,1:dimy,1:dimz)
        read(3) Jy(1:dimx,1:dimy,1:dimz)
        read(3) Jz(1:dimx,1:dimy,1:dimz)
        close(3)
      else
!$omp parallel do
        do k=1,dimz
        do j=1,dimy
        do i=1,dimx
           Jx(i,j,k) = yderiv(Bz,i,j,k) - zderiv(By,i,j,k) 
           Jy(i,j,k) = zderiv(Bx,i,j,k) - xderiv(Bz,i,j,k)
           Jz(i,j,k) = xderiv(By,i,j,k) - yderiv(Bx,i,j,k)
        enddo
        enddo
        enddo
!$omp end parallel do
      end if
    else
      write(*,'(A)')'| File '//BfieldName//' does not exist'
      stop
    end if
  end subroutine read_binary

  subroutine read_data()
    if(usemask) then
      write(*,'(A)')'| Loading mask data...'
      mask = mask*0.0d0
      open(3,File=MaskFile,Access="stream",Form = "unformatted" )
      read(3) mask(1: mask_level*(dimx-1) + 1, 1:mask_level*(dimy-1) + 1)
      close(3)
      write(*,'(A)')'| Data sucessfully loaded !'
   end if

  select case(indataformat)
    case('blk')
        call read_blk()
    case('binary')
        call read_binary()
    case default
        write(*,*)'no data input'
        stop
    end select
    ! --------- get the value of B in the ghost cell ------------
    ! surface in ghost shell
    Bx(0,1:dimy,1:dimz) = Bx(1,1:dimy,1:dimz)
    By(0,1:dimy,1:dimz) = By(1,1:dimy,1:dimz)
    Bz(0,1:dimy,1:dimz) = Bz(1,1:dimy,1:dimz)
      
    Bx(dimx+1,1:dimy,1:dimz) = Bx(dimx,1:dimy,1:dimz)
    By(dimx+1,1:dimy,1:dimz) = By(dimx,1:dimy,1:dimz)
    Bz(dimx+1,1:dimy,1:dimz) = Bz(dimx,1:dimy,1:dimz)
      
    Bx(1:dimx,0,1:dimz) = Bx(1:dimx,1,1:dimz)
    By(1:dimx,0,1:dimz) = By(1:dimx,1,1:dimz)
    Bz(1:dimx,0,1:dimz) = Bz(1:dimx,1,1:dimz)
      
    Bx(1:dimx,dimy+1,1:dimz) = Bx(1:dimx,dimy,1:dimz)
    By(1:dimx,dimy+1,1:dimz) = By(1:dimx,dimy,1:dimz)
    Bz(1:dimx,dimy+1,1:dimz) = Bz(1:dimx,dimy,1:dimz)
      
    Bx(1:dimx,1:dimy,0) = Bx(1:dimx,1:dimy,1)
    By(1:dimx,1:dimy,0) = By(1:dimx,1:dimy,1)
    Bz(1:dimx,1:dimy,0) = Bz(1:dimx,1:dimy,1)
      
    Bx(1:dimx,1:dimy,dimz+1) = Bx(1:dimx,1:dimy,dimz)
    By(1:dimx,1:dimy,dimz+1) = By(1:dimx,1:dimy,dimz)
    Bz(1:dimx,1:dimy,dimz+1) = Bz(1:dimx,1:dimy,dimz)
    
    ! lines in ghost shell
    Bx(0,0,1:dimz) = Bx(1,1,1:dimz)
    By(0,0,1:dimz) = By(1,1,1:dimz)
    Bz(0,0,1:dimz) = Bz(1,1,1:dimz)
      
    Bx(dimx+1,0,1:dimz) = Bx(dimx,1,1:dimz)
    By(dimx+1,0,1:dimz) = By(dimx,1,1:dimz)
    Bz(dimx+1,0,1:dimz) = Bz(dimx,1,1:dimz)
      
    Bx(0,dimy+1,1:dimz) = Bx(1,dimy,1:dimz)
    By(0,dimy+1,1:dimz) = By(1,dimy,1:dimz)
    Bz(0,dimy+1,1:dimz) = Bz(1,dimy,1:dimz)
      
    Bx(dimx+1,dimy+1,1:dimz) = Bx(dimx,dimy,1:dimz)
    By(dimx+1,dimy+1,1:dimz) = By(dimx,dimy,1:dimz)
    Bz(dimx+1,dimy+1,1:dimz) = Bz(dimx,dimy,1:dimz)
    ! - - - - - - - - - - - - - - - - - - - - - - - 
    Bx(1:dimx,0,0) = Bx(1:dimx,1,1)
    By(1:dimx,0,0) = By(1:dimx,1,1)
    Bz(1:dimx,0,0) = Bz(1:dimx,1,1)
      
    Bx(1:dimx,dimy+1,0) = Bx(1:dimx,dimy,1)
    By(1:dimx,dimy+1,0) = By(1:dimx,dimy,1)
    Bz(1:dimx,dimy+1,0) = Bz(1:dimx,dimy,1)
      
    Bx(1:dimx,dimy+1,dimz+1) = Bx(1:dimx,dimy,dimz)
    By(1:dimx,dimy+1,dimz+1) = By(1:dimx,dimy,dimz)
    Bz(1:dimx,dimy+1,dimz+1) = Bz(1:dimx,dimy,dimz)
      
    Bx(1:dimx,0,dimz+1) = Bx(1:dimx,1,dimz)
    By(1:dimx,0,dimz+1) = By(1:dimx,1,dimz)
    Bz(1:dimx,0,dimz+1) = Bz(1:dimx,1,dimz)
    ! - - - - - - - - - - - - - - - - - - - - - - - 
    Bx(0,1:dimy,0) = Bx(1,1:dimy,1)
    By(0,1:dimy,0) = By(1,1:dimy,1)
    Bz(0,1:dimy,0) = Bz(1,1:dimy,1)
      
    Bx(dimx+1,1:dimy,0) = Bx(dimx,1:dimy,1)
    By(dimx+1,1:dimy,0) = By(dimx,1:dimy,1)
    Bz(dimx+1,1:dimy,0) = Bz(dimx,1:dimy,1)
      
    Bx(0,1:dimy,dimz+1) = Bx(1,1:dimy,dimz)
    By(0,1:dimy,dimz+1) = By(1,1:dimy,dimz)
    Bz(0,1:dimy,dimz+1) = Bz(1,1:dimy,dimz)
      
    Bx(dimx+1,1:dimy,dimz+1) = Bx(dimx,1:dimy,dimz)
    By(dimx+1,1:dimy,dimz+1) = By(dimx,1:dimy,dimz)
    Bz(dimx+1,1:dimy,dimz+1) = Bz(dimx,1:dimy,dimz)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - 
    ! points in ghost shell
    Bx(0,0,0) = Bx(1,1,1)
    By(0,0,0) = By(1,1,1)
    Bz(0,0,0) = Bz(1,1,1)
      
    Bx(dimx+1,0,0) = Bx(dimx,1,1)
    By(dimx+1,0,0) = By(dimx,1,1)
    Bz(dimx+1,0,0) = Bz(dimx,1,1)
      
    Bx(0,dimy+1,0) = Bx(1,dimy,1)
    By(0,dimy+1,0) = By(1,dimy,1)
    Bz(0,dimy+1,0) = Bz(1,dimy,1)
      
    Bx(0,0,dimz+1) = Bx(1,1,dimz)
    By(0,0,dimz+1) = By(1,1,dimz)
    Bz(0,0,dimz+1) = Bz(1,1,dimz)
      
    Bx(dimx+1,dimy+1,0) = Bx(dimx,dimy,1)
    By(dimx+1,dimy+1,0) = By(dimx,dimy,1)
    Bz(dimx+1,dimy+1,0) = Bz(dimx,dimy,1)
      
    Bx(dimx+1,0,dimz+1) = Bx(dimx,1,dimz)
    By(dimx+1,0,dimz+1) = By(dimx,1,dimz)
    Bz(dimx+1,0,dimz+1) = Bz(dimx,1,dimz)
      
    Bx(0,dimy+1,dimz+1) = Bx(1,dimy,dimz)
    By(0,dimy+1,dimz+1) = By(1,dimy,dimz)
    Bz(0,dimy+1,dimz+1) = Bz(1,dimy,dimz)
      
    Bx(dimx+1,dimy+1,dimz+1) = Bx(dimx,dimy,dimz)
    By(dimx+1,dimy+1,dimz+1) = By(dimx,dimy,dimz)
    Bz(dimx+1,dimy+1,dimz+1) = Bz(dimx,dimy,dimz)
    
    ! ----------- set the value of Jxyz in the ghost cell ---------
    ! surface in ghost shell
    Jx(0,1:dimy,1:dimz) = Jx(1,1:dimy,1:dimz)
    Jy(0,1:dimy,1:dimz) = Jy(1,1:dimy,1:dimz)
    Jz(0,1:dimy,1:dimz) = Jz(1,1:dimy,1:dimz)
      
    Jx(dimx+1,1:dimy,1:dimz) = Jx(dimx,1:dimy,1:dimz)
    Jy(dimx+1,1:dimy,1:dimz) = Jy(dimx,1:dimy,1:dimz)
    Jz(dimx+1,1:dimy,1:dimz) = Jz(dimx,1:dimy,1:dimz)
      
    Jx(1:dimx,0,1:dimz) = Jx(1:dimx,1,1:dimz)
    Jy(1:dimx,0,1:dimz) = Jy(1:dimx,1,1:dimz)
    Jz(1:dimx,0,1:dimz) = Jz(1:dimx,1,1:dimz)
      
    Jx(1:dimx,dimy+1,1:dimz) = Jx(1:dimx,dimy,1:dimz)
    Jy(1:dimx,dimy+1,1:dimz) = Jy(1:dimx,dimy,1:dimz)
    Jz(1:dimx,dimy+1,1:dimz) = Jz(1:dimx,dimy,1:dimz)
      
    Jx(1:dimx,1:dimy,0) = Jx(1:dimx,1:dimy,1)
    Jy(1:dimx,1:dimy,0) = Jy(1:dimx,1:dimy,1)
    Jz(1:dimx,1:dimy,0) = Jz(1:dimx,1:dimy,1)
      
    Jx(1:dimx,1:dimy,dimz+1) = Jx(1:dimx,1:dimy,dimz)
    Jy(1:dimx,1:dimy,dimz+1) = Jy(1:dimx,1:dimy,dimz)
    Jz(1:dimx,1:dimy,dimz+1) = Jz(1:dimx,1:dimy,dimz)
    
    ! lines in ghost shell
    Jx(0,0,1:dimz) = Jx(1,1,1:dimz)
    Jy(0,0,1:dimz) = Jy(1,1,1:dimz)
    Jz(0,0,1:dimz) = Jz(1,1,1:dimz)
      
    Jx(dimx+1,0,1:dimz) = Jx(dimx,1,1:dimz)
    Jy(dimx+1,0,1:dimz) = Jy(dimx,1,1:dimz)
    Jz(dimx+1,0,1:dimz) = Jz(dimx,1,1:dimz)
      
    Jx(0,dimy+1,1:dimz) = Jx(1,dimy,1:dimz)
    Jy(0,dimy+1,1:dimz) = Jy(1,dimy,1:dimz)
    Jz(0,dimy+1,1:dimz) = Jz(1,dimy,1:dimz)
      
    Jx(dimx+1,dimy+1,1:dimz) = Jx(dimx,dimy,1:dimz)
    Jy(dimx+1,dimy+1,1:dimz) = Jy(dimx,dimy,1:dimz)
    Jz(dimx+1,dimy+1,1:dimz) = Jz(dimx,dimy,1:dimz)
    ! - - - - - - - - - - - - - - - - - - - - - - - 
    Jx(1:dimx,0,0) = Jx(1:dimx,1,1)
    Jy(1:dimx,0,0) = Jy(1:dimx,1,1)
    Jz(1:dimx,0,0) = Jz(1:dimx,1,1)
      
    Jx(1:dimx,dimy+1,0) = Jx(1:dimx,dimy,1)
    Jy(1:dimx,dimy+1,0) = Jy(1:dimx,dimy,1)
    Jz(1:dimx,dimy+1,0) = Jz(1:dimx,dimy,1)
      
    Jx(1:dimx,dimy+1,dimz+1) = Jx(1:dimx,dimy,dimz)
    Jy(1:dimx,dimy+1,dimz+1) = Jy(1:dimx,dimy,dimz)
    Jz(1:dimx,dimy+1,dimz+1) = Jz(1:dimx,dimy,dimz)
      
    Jx(1:dimx,0,dimz+1) = Jx(1:dimx,1,dimz)
    Jy(1:dimx,0,dimz+1) = Jy(1:dimx,1,dimz)
    Jz(1:dimx,0,dimz+1) = Jz(1:dimx,1,dimz)
    ! - - - - - - - - - - - - - - - - - - - - - - - 
    Jx(0,1:dimy,0) = Jx(1,1:dimy,1)
    Jy(0,1:dimy,0) = Jy(1,1:dimy,1)
    Jz(0,1:dimy,0) = Jz(1,1:dimy,1)
      
    Jx(dimx+1,1:dimy,0) = Jx(dimx,1:dimy,1)
    Jy(dimx+1,1:dimy,0) = Jy(dimx,1:dimy,1)
    Jz(dimx+1,1:dimy,0) = Jz(dimx,1:dimy,1)
      
    Jx(0,1:dimy,dimz+1) = Jx(1,1:dimy,dimz)
    Jy(0,1:dimy,dimz+1) = Jy(1,1:dimy,dimz)
    Jz(0,1:dimy,dimz+1) = Jz(1,1:dimy,dimz)
      
    Jx(dimx+1,1:dimy,dimz+1) = Jx(dimx,1:dimy,dimz)
    Jy(dimx+1,1:dimy,dimz+1) = Jy(dimx,1:dimy,dimz)
    Jz(dimx+1,1:dimy,dimz+1) = Jz(dimx,1:dimy,dimz)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - 
    ! points in ghost shell
    Jx(0,0,0) = Jx(1,1,1)
    Jy(0,0,0) = Jy(1,1,1)
    Jz(0,0,0) = Jz(1,1,1)
      
    Jx(dimx+1,0,0) = Jx(dimx,1,1)
    Jy(dimx+1,0,0) = Jy(dimx,1,1)
    Jz(dimx+1,0,0) = Jz(dimx,1,1)
      
    Jx(0,dimy+1,0) = Jx(1,dimy,1)
    Jy(0,dimy+1,0) = Jy(1,dimy,1)
    Jz(0,dimy+1,0) = Jz(1,dimy,1)
      
    Jx(0,0,dimz+1) = Jx(1,1,dimz)
    Jy(0,0,dimz+1) = Jy(1,1,dimz)
    Jz(0,0,dimz+1) = Jz(1,1,dimz)
      
    Jx(dimx+1,dimy+1,0) = Jx(dimx,dimy,1)
    Jy(dimx+1,dimy+1,0) = Jy(dimx,dimy,1)
    Jz(dimx+1,dimy+1,0) = Jz(dimx,dimy,1)
      
    Jx(dimx+1,0,dimz+1) = Jx(dimx,1,dimz)
    Jy(dimx+1,0,dimz+1) = Jy(dimx,1,dimz)
    Jz(dimx+1,0,dimz+1) = Jz(dimx,1,dimz)
      
    Jx(0,dimy+1,dimz+1) = Jx(1,dimy,dimz)
    Jy(0,dimy+1,dimz+1) = Jy(1,dimy,dimz)
    Jz(0,dimy+1,dimz+1) = Jz(1,dimy,dimz)
      
    Jx(dimx+1,dimy+1,dimz+1) = Jx(dimx,dimy,dimz)
    Jy(dimx+1,dimy+1,dimz+1) = Jy(dimx,dimy,dimz)
    Jz(dimx+1,dimy+1,dimz+1) = Jz(dimx,dimy,dimz)

  end subroutine read_data

  subroutine write_data()
    write(*,'(A)')' Writing the result ...'
    open(4,File=OutFileName,Access="stream",STATUS="REPLACE",&
    &Form = "unformatted" )
    write(4) cal_data
    close(4)
  end subroutine write_data

end module mod_io
