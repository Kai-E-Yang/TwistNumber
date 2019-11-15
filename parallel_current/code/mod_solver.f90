module mod_solver
use mod_param
use mod_operator
use omp_lib
implicit none
contains
!-------- the main loop for calculation --------
  subroutine do_cal()
    integer(kind=selected_int_kind(8)) :: TID
    integer,dimension(:),allocatable :: mark
    integer :: sum0,flag(5),k
    integer,dimension(8) :: time_begin,time_end
    real(kind=r8) :: time_delta(8)
    
    ! Set parallel threads
    if(nthreads == 0) then
      nthreads = OMP_GET_MAX_THREADS()
    end if
    
    call OMP_SET_NUM_THREADS(nthreads)
    allocate(mark(0:(nthreads-1)))
    
    mark = 0.0d0
    sum0 = 0
    do k=1,5
       flag(k) = floor(0.2*k*refine_dimx)
    end do
    
    write(*,'(A)')'| Start the computation !'
    call date_and_time(VALUES=time_begin)
    
!$omp parallel private(TID, sum0, k)
!$omp do schedule(dynamic,1)
    
    do k=1,refine_dimx
       call cal_plane(k)
       TID = OMP_GET_THREAD_NUM()
       mark(TID) = mark(TID) + 1
       sum0 = sum(mark)
       if(sum0 .eq. flag(1) .or. sum0 .eq. flag(2) .or. sum0 &
         & .eq. flag(3) .or. sum0 .eq. flag(4) .or. sum0 .eq.  &
         & flag(5)) then
         call date_and_time(VALUES=time_end)
         time_delta = 1.0d0*(time_end - time_begin)
         write(*,'(A,f6.1,A,A,f8.2)')'| Equilibrium loop progress: '&
         & ,100.*sum(mark)/(refine_dimx*1.0),'%','; costs time (min): '&
         & ,(time_delta(7)/3600.0 + time_delta(6)/60. &
         & + time_delta(5) + time_delta(3)*24.0)*60
       end if
    end do
    
!$omp end do
!$omp end parallel
    
    call date_and_time(VALUES=time_end)
    write(*,'(A)')' Computation finished !'
  end subroutine do_cal
! ------- do parallel cal on each layer ------
  subroutine cal_plane(k)
    integer,intent(in)::k
    integer::i,j
    real(kind=r8) :: t_slice(refine_dimy,refine_dimz)
    real(kind=r8) :: end_slice(refine_dimy,refine_dimz)
    real(kind=r8) :: length_slice(refine_dimy,refine_dimz)
    real(kind=r8) :: xindex,yindex,zindex
    
    xindex = x_start+(k-1)/real(nlevel,kind=r8)
    do j=1,refine_dimz
    do i=1,refine_dimy
       yindex = y_start+(i-1)/real(nlevel,kind=r8)
       zindex = z_start+(j-1)/real(nlevel,kind=r8)
        ! print*,'x,y,z ',xindex,yindex,zindex
       call cal_point(xindex,yindex,zindex,t_slice(i,j),&
            &end_slice(i,j),length_slice(i,j))
    end do
    end do
    cal_data(k,:,:,1) = t_slice
    cal_data(k,:,:,2) = length_slice
    cal_data(k,:,:,3) = end_slice
  end subroutine cal_plane
! ------ do calculation at each point ------
  subroutine cal_point(PosiX,PosiY,PosiZ,Twist,end_mark,length)
    real(kind=r8),intent(in) :: PosiX
    real(kind=r8),intent(in) :: PosiY
    real(kind=r8),intent(in) :: PosiZ
    real(kind=r8),intent(out) :: Twist
    real(kind=r8),intent(out) :: end_mark
    real(kind=r8) :: length
    real(kind=r8) :: mask1,mask2
    ! linep is the present position
    ! linef is the forward end of the integral line
    ! linef is the backward end of the integral line
    ! linep,f,b contains the value in the order of 
    ! 1,2,3: x,y,z coordinate, 
    ! 4,5,6, 7,8,9: x,y,z components of U, V vactors,
    ! 10,11,12,13: xyz components of the tangent vector and the norm of B
    ! at the end of the line.
    real(kind=r8) :: LineP(4)
    real(kind=r8) :: LineF(4)
    real(kind=r8) :: LineB(4)
    
    LineF = real(0,kind=r8)
    LineB = real(0,kind=r8)
    LineP=(/PosiX,PosiY,PosiZ,real(0,kind=r8)/)
    call fieldline(LineP,LineF,LineB,length)
    mask1 = 1.0d0
    mask2 = 1.0d0
    if(usemask) then
      call cal_mask(LineP(1:3),mask1)
      call cal_mask(LineF(1:3),mask2)
    end if
    Twist = LineF(4) - LineB(4)
    if(mask1 .gt. 0.5 .or. mask2 .gt. 0.5) then
    ! ------------ get the mark flag for the two end ---------
    ! meaning of the value of end_mark
    ! 1: close field line
    ! 2: open field line, positive end roots in bottom
    ! 3: open field line, negative end roots in bottom
    ! 4: field line without end roots in the bottom boundary
      end_mark = real(0,kind=r8)
      if(abs(LineF(3)) .lt. BoundaryEps .and. abs(LineB(3)) &
        & .lt. BoundaryEps) then
        end_mark = real(1,kind=r8)
      end if
    
      if(abs(LineB(3)) .lt. BoundaryEps .and. abs(LineF(3)) &
        & .gt. BoundaryEps) then
        end_mark = real(2,kind=r8)
      end if
    
      if(abs(LineF(3)) .lt. BoundaryEps .and. abs(LineB(3)) &
        & .gt. BoundaryEps) then
        end_mark = real(3,kind=r8)
      end if
    
      if(abs(LineF(3)) .gt. BoundaryEps .and. abs(LineB(3)) &
        & .gt. BoundaryEps) then
        end_mark = real(4,kind=r8)
      end if
    
    else
       Twist      = real(0,kind=r8)
       length     = real(0,kind=r8)
       end_mark   = real(0,kind=r8)
    end if
  end subroutine cal_point
! -------- calculate the field line ---------
  subroutine fieldline(LineP,LineF,LineB,length)
    real(kind=r8)::LineP(4)
    real(kind=r8)::LineF(4)
    real(kind=r8)::LineB(4)
    real(kind=r8)::sig,lengthF,lengthB,length
    sig = real(1,kind=r8)
    call integralLine(LineP,LineF,sig,lengthF)
    sig = real(-1,kind=r8)
    call integralLine(LineP,LineB,sig,lengthB)
    length = abs(lengthB) + abs(lengthF)
  end subroutine fieldline
! -------- do the line integral --------
  subroutine integralLine(LineP,LineI,sig,s_end)
    integer,parameter :: neqn = 4
    integer::n_step,flag
    real(kind=r8),intent(in)::sig
    real(kind=r8)::LineP(4)
    real(kind=r8)::LineI(4)
    real(kind=r8)::tmp(4)
    real(kind=r8)::s_start,s_end,ds
    ! external rhs
    s_start = real(0,kind=r8)
    s_end   = real(0,kind=r8)
    tmp     = LineP
    ds      = real(delta_s,kind=r8)*sig
    n_step  = 0
    flag    = 1
    
    do while (n_step .le. isn .and. flag .eq. 1)
       s_end   = s_start + ds
       call rk4(rhs,neqn,s_start,s_end,tmp,flag)
       
       s_start = s_end
       n_step  = n_step+1
    end do
    LineI(1:4) = tmp
    end subroutine integralLine
end module mod_solver