module mod_param
implicit none
  save
  integer, parameter :: i4=SELECTED_INT_KIND(R=9) 
  integer, parameter :: i8=SELECTED_INT_KIND(R=18) 
  integer, parameter :: r8=SELECTED_REAL_KIND(P=15,R=300)
  character*100 :: par
  character*100 :: indataformat
  character*100 :: BfieldName, JfieldName, OutFileName, MaskFile

  integer :: dimx,dimy,dimz
  integer :: refine_dimx,refine_dimy,refine_dimz
  integer :: x_start,y_start,z_start
  integer :: x_end,y_end,z_end
  integer :: nlevel
  integer :: nthreads
  integer :: isn=10000
  integer :: iterMax
  integer :: mask_level

  real(kind=r8) :: delta_s
  real(kind=r8),dimension(:,:,:),allocatable :: Bx,By,Bz
  real(kind=r8),dimension(:,:,:),allocatable :: Jx,Jy,Jz
  real(kind=r8),dimension(:,:,:),allocatable :: x,y,z

  real(kind=r8),dimension(:,:,:,:),allocatable :: cal_data
  real(kind=r8),dimension(:,:),allocatable :: mask

  real(kind=r8) :: BoundaryEps
  real(kind=r8) :: eps_B

  logical :: usemask,usecurrent

  namelist /filename_par/ BfieldName, JfieldName, OutFileName, MaskFile, indataformat
  namelist /cal_par/ nthreads,dimx,dimy,dimz,x_start,y_start, &
    z_start,x_end,y_end,z_end,nlevel,delta_s,mask_level,&
    usemask,usecurrent
end module mod_param