module mod_operator
use mod_param
implicit none
contains
! -------- function to calculate the cross production
  function cross(a, b)
    real(kind = r8),dimension(3) :: cross
    real(kind = r8),dimension(3),intent(in) :: a, b
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross
! ------- function check whether the point outside the computational domain
  function outside_boundary(posi)
    logical:: outside_boundary
    real(kind = r8),intent(in) :: posi(3) 
    if(posi(1) .gt. 1 .and. posi(1) .lt. dimx &
    & .and. posi(2) .gt. 1 .and. posi(2) .lt. dimy &
    & .and. posi(3) .gt. 1 .and. posi(3) .lt. dimz ) then
      outside_boundary = .false.
    else
      outside_boundary = .true.
    end if
  end function outside_boundary
! -------- function give derivtive in X direction -------------
  function xderiv(array,i,j,k)
    real(kind = r8) :: xderiv
    real(kind = r8), dimension(:,:,:), intent(in) :: array
    integer,intent(in) :: i,j,k
    if ((i .ne. 1) .and. (i .ne. dimx)) then
      xderiv=0.5*(array(i+1,j,k)-array(i-1,j,k))
    else if (i .eq. 1) then
      xderiv=(-3.*array(1,j,k)+4.*array(2,j,k)-array(3,j,k))/real(2,kind=r8)
    else if (i .eq. dimx) then
      xderiv=(3.*array(dimx,j,k)-4.*array(dimx-1,j,k)+array(dimx-2,j,k))/real(2,kind=r8)
    end if
  end function xderiv
! -------- function give derivtive in Y direction -------------
  function yderiv(array,i,j,k)
    real(kind = r8) :: yderiv
    real(kind = r8), dimension(:,:,:), intent(in) :: array
    integer,intent(in) :: i,j,k
    if ((j .ne. 1) .and. (j .ne. dimy)) then
      yderiv=0.5*(array(i,j+1,k)-array(i,j-1,k))
    else if (j .eq. 1) then
      yderiv=(-3.*array(i,1,k)+4.*array(i,2,k)-array(i,3,k))/real(2,kind=r8)
    else if (j .eq. dimy) then
      yderiv=(3.*array(i,dimy,k)-4.*array(i,dimy-1,k)+array(i,dimy-2,k))/real(2,kind=r8)
    end if
  end function yderiv
! -------- function give derivtive in Z direction -------------
  function zderiv(array,i,j,k)
    real(kind = r8) :: zderiv
    real(kind = r8), dimension(:,:,:), intent(in) :: array
    integer,intent(in) :: i,j,k
    if ((k .ne. 1) .and. (k .ne. dimz)) then
      zderiv=0.5*(array(i,j,k+1)-array(i,j,k-1))
    else if (k .eq. 1) then
      zderiv=(-3.*array(i,j,1)+4.*array(i,j,2)-array(i,j,3))/real(2,kind=r8)
    else if (k .eq. dimz) then
      zderiv=(3.*array(i,j,dimz)-4.*array(i,j,dimz-1)+array(i,j,dimz-2))/real(2,kind=r8)
    end if
  end function zderiv
! ---------------------------------------------------------------
! the first dimension of xv indicates the x,y,z components of B
! the location in the second dimension of xv 
! 
!               g       g
!              /       /
!             /       /
!    g- - - -8- - - -7- - - -g
!           /|      /|   
!          / |     / |  
! g- - - -5- - - -6  |
!   g- - -|- 4 - -|- 3- - - -g
!         | /     | /
!         |/      |/     
! g- - - -1- - - -2- - - -g
!         |       |
!         |       |
!         |       |
!         g- - - -g
! 
! 
! ------ calculate the Bxyz on 8 vertex -----------
  subroutine corner (nnn,xv)
    integer,intent(in)::nnn(3)
    integer::n1,n2,n3
    integer::n1p,n2p,n3p
    real ( kind = r8 ),intent(out)::xv(3,8)
    n1 = nnn(1)
    n2 = nnn(2)
    n3 = nnn(3)
    
    n1p = nnn(1)+1
    n2p = nnn(2)+1
    n3p = nnn(3)+1
    
    xv(1,1) = Bx(n1,n2,n3)
    xv(1,2) = Bx(n1p,n2,n3)
    xv(1,3) = Bx(n1p,n2p,n3)
    xv(1,4) = Bx(n1,n2p,n3)
    xv(1,5) = Bx(n1,n2,n3p)
    xv(1,6) = Bx(n1p,n2,n3p)
    xv(1,7) = Bx(n1p,n2p,n3p)
    xv(1,8) = Bx(n1,n2p,n3p)
    
    xv(2,1) = By(n1,n2,n3)
    xv(2,2) = By(n1p,n2,n3)
    xv(2,3) = By(n1p,n2p,n3)
    xv(2,4) = By(n1,n2p,n3)
    xv(2,5) = By(n1,n2,n3p)
    xv(2,6) = By(n1p,n2,n3p)
    xv(2,7) = By(n1p,n2p,n3p)
    xv(2,8) = By(n1,n2p,n3p)
    
    xv(3,1) = Bz(n1,n2,n3)
    xv(3,2) = Bz(n1p,n2,n3)
    xv(3,3) = Bz(n1p,n2p,n3)
    xv(3,4) = Bz(n1,n2p,n3)
    xv(3,5) = Bz(n1,n2,n3p)
    xv(3,6) = Bz(n1p,n2,n3p)
    xv(3,7) = Bz(n1p,n2p,n3p)
    xv(3,8) = Bz(n1,n2p,n3p)
  end subroutine corner
! ------ calculate the Jxyz on 8 vertex -----------
  subroutine corner_J (nnn,xv)
    integer,intent(in)::nnn(3)
    integer::n1,n2,n3
    integer::n1p,n2p,n3p
    real ( kind = r8 ),intent(out)::xv(3,8)
    n1 = nnn(1)
    n2 = nnn(2)
    n3 = nnn(3)
    
    n1p = nnn(1)+1
    n2p = nnn(2)+1
    n3p = nnn(3)+1
    
    xv(1,1) = Jx(n1,n2,n3)
    xv(1,2) = Jx(n1p,n2,n3)
    xv(1,3) = Jx(n1p,n2p,n3)
    xv(1,4) = Jx(n1,n2p,n3)
    xv(1,5) = Jx(n1,n2,n3p)
    xv(1,6) = Jx(n1p,n2,n3p)
    xv(1,7) = Jx(n1p,n2p,n3p)
    xv(1,8) = Jx(n1,n2p,n3p)
    
    xv(2,1) = Jy(n1,n2,n3)
    xv(2,2) = Jy(n1p,n2,n3)
    xv(2,3) = Jy(n1p,n2p,n3)
    xv(2,4) = Jy(n1,n2p,n3)
    xv(2,5) = Jy(n1,n2,n3p)
    xv(2,6) = Jy(n1p,n2,n3p)
    xv(2,7) = Jy(n1p,n2p,n3p)
    xv(2,8) = Jy(n1,n2p,n3p)
    
    xv(3,1) = Jz(n1,n2,n3)
    xv(3,2) = Jz(n1p,n2,n3)
    xv(3,3) = Jz(n1p,n2p,n3)
    xv(3,4) = Jz(n1,n2p,n3)
    xv(3,5) = Jz(n1,n2,n3p)
    xv(3,6) = Jz(n1p,n2,n3p)
    xv(3,7) = Jz(n1p,n2p,n3p)
    xv(3,8) = Jz(n1,n2p,n3p)
  end subroutine corner_J
! --------- get the tri-linear interpolation ----------
  subroutine xitp (interp,xv,dxyz1)
    real(kind=r8),intent(out)::interp(3)
    real(kind=r8)::dxyz2(3),weight(8)
    real(kind=r8)::xv(3,8)
    real(kind=r8),intent(in)::dxyz1(3)
    
    ! the location of each index 
    !     8- - - -7
    !    /|      /|   
    !   / |     / |  
    !  5- - - -6  |
    !  |  4 - -|- 3   
    !  | /     | /
    !  |/      |/     
    !  1- - - -2
    ! 
    ! o indicates the location of dxyz in the unit square
    ! (0,1)---------(1,1)
    !     |   |     |
    !     |   |     |
    !     |---o-----|
    !     |   |     |
    ! (0,0)---------(1,0)
    ! the left bottom one is dxyz1
    ! the right up one is dxyz2=1-dxyz
    
    dxyz2 = real(1,kind=r8) - dxyz1
    
    weight(1) = dxyz2(1)*dxyz2(2)*dxyz2(3)
    weight(2) = dxyz1(1)*dxyz2(2)*dxyz2(3)
    weight(3) = dxyz1(1)*dxyz1(2)*dxyz2(3)
    weight(4) = dxyz2(1)*dxyz1(2)*dxyz2(3)
    
    weight(5) = dxyz2(1)*dxyz2(2)*dxyz1(3)
    weight(6) = dxyz1(1)*dxyz2(2)*dxyz1(3)
    weight(7) = dxyz1(1)*dxyz1(2)*dxyz1(3)
    weight(8) = dxyz2(1)*dxyz1(2)*dxyz1(3)
    
    interp(1) = dot_product(xv(1,:),weight)
    interp(2) = dot_product(xv(2,:),weight)
    interp(3) = dot_product(xv(3,:),weight)
    
  end subroutine xitp
! ------ give the equation of ODE ------
  subroutine diffLine(Posi,Tangent)
    real(kind=r8),intent(in)::Posi(3)
    real(kind=r8),intent(out)::Tangent(4)
    real(kind=r8)::dxyz(3)
    real(kind=r8)::xv(3,8)
    real(kind=r8)::bvec(3)
    real(kind=r8)::binter
    real(kind=r8)::eps
    integer :: nnn(3)
    
    eps = epsilon(real(1,kind=r8))
    dxyz = Posi - floor(Posi)
    nnn = floor(Posi)
    
    call corner(nnn,xv)
    call xitp(bvec,xv,dxyz)
    
    binter = sqrt(dot_product(bvec,bvec))
    
    if (binter .lt. eps) then
       Tangent = real(0,kind=r8)*(/1,1,1,1/)
    else
       Tangent(1:3) = bvec/binter
       Tangent(4)   = binter
    end if
  end subroutine diffLine

  subroutine diff_twist(Posi,Tangent)
    real(kind=r8),intent(in)::Posi(3)
    real(kind=r8),intent(out)::Tangent(4)
    real(kind=r8)::dxyz(3)
    real(kind=r8)::xv(3,8)
    real(kind=r8)::jvec(3)
    real(kind=r8)::jinter
    real(kind=r8)::eps
    integer :: nnn(3)
    
    eps = epsilon(real(1,kind=r8))
    dxyz = Posi - floor(Posi)
    nnn = floor(Posi)
    
    call corner_J(nnn,xv)
    call xitp(jvec,xv,dxyz)
    
    jinter = sqrt(dot_product(jvec,jvec))
    
    if (jinter .lt. eps) then
       Tangent = real(0,kind=r8)*(/1,1,1,1/)
    else
       Tangent(1:3) = jvec/jinter
       Tangent(4)   = jinter
    end if
  end subroutine diff_twist

! ---- the ODE solver runge-kutta 4 order method -----
  subroutine rk4 (f, neqn, t, t_out, y, flag)
    integer(kind = i4) neqn
    external f
    real(kind=r8) :: t,t_out,dt
    real(kind=r8) :: y(neqn)
    real(kind=r8) :: y_tmp(neqn)
    real(kind=r8) :: k1(neqn)
    real(kind=r8) :: k2(neqn)
    real(kind=r8) :: k3(neqn)
    real(kind=r8) :: k4(neqn)
    real(kind=r8) :: yp(neqn)
    real(kind=r8) :: Tangent(4)
    integer :: flag
    integer :: iterNum
    real(kind=r8) :: eps
    
    eps = epsilon(real(1,kind=r8))
    
    iterNum = 0
    flag = 1
    dt = t_out - t
    y_tmp = y
    
    call f(t,y,yp)
    k1 = dt*yp
    
    call f(t+0.5*dt,y + 0.5*k1,yp)
    k2 = dt*yp
    
    call f(t+0.5*dt,y+0.5*k2,yp)
    k3 = dt*yp
    
    call f(t+dt,y+k3,yp)
    k4 = dt*yp
    
    y_tmp =  y &
        & +( k1 &
        & +  k2*real(2,kind=r8) &
        & +  k3*real(2,kind=r8) &
        & +  k4 )/6.0d0
    
    do while(outside_boundary(y_tmp(1:3)) .and. iterNum .le. iterMax)
    
       iterNum = iterNum + 1
       dt = dt * real(0.5,kind=r8)
       
       call f(t,y,yp)
       k1 = dt*yp
    
       call f(t+0.5*dt,y + 0.5*k1,yp)
       k2 = dt*yp
    
       call f(t+0.5*dt,y+0.5*k2,yp)
       k3 = dt*yp
    
       call f(t+dt,y+k3,yp)
       k4 = dt*yp
    
       y_tmp =  y &
           & +( k1 &
           & +  k2*real(2,kind=r8) &
           & +  k3*real(2,kind=r8) &
           & +  k4 )/6.0d0
    end do
    
    y = y_tmp
    call f(t,y,yp)
    
    call diffLine(y(1:3),Tangent)
    
    if(iterNum .ge. iterMax .or. Tangent(4) .le. eps) then
      flag = 0
    end if
    
    t_out = t + dt
  end subroutine rk4
! ----------- R.H.S. of the ODE ----------
  subroutine rhs(t,y,yp)
    integer,parameter :: neqn = 4
    real(kind=r8)::t
    real(kind=r8)::y(neqn)
    real(kind=r8)::yp(neqn)
    real(kind=r8)::TangentB(4)
    real(kind=r8)::TangentJ(4)
    
    call diffLine(y(1:3),TangentB)
    call diff_twist(y(1:3),TangentJ)
    
    yp(1:3) = TangentB(1:3)
    yp(4) = dot_product(TangentJ(1:3),TangentB(1:3))*TangentJ(4)/TangentB(4)

  end subroutine rhs
! ------ if the user provide some mask value on the bottom -----
  subroutine cal_mask(posi,mask_tmp)
    real(kind=r8),intent(in) ::posi(3)
    real(kind=r8),intent(out)::mask_tmp
    real(kind=r8)::weight(4),MaskCom(4)
    real(kind=r8)::xindex,yindex
    integer::i,j
    mask_tmp = 0.0d0
    ! --------------- posi1 --------------
    xindex = posi(1)
    yindex = posi(2)
    
    i = floor( (xindex - 1)*mask_level + 1 )
    j = floor( (yindex - 1)*mask_level + 1 )
    if (xindex .lt. 1 .or. xindex .gt. dimx .or. yindex .lt. 1 .or. yindex .gt. dimy) THEN
      mask_tmp=0.0d0
    else
    ! calculate the mask component
      MaskCom(1) = mask(i,j)
      MaskCom(2) = mask(i+1,j)
      MaskCom(3) = mask(i,j+1)
      MaskCom(4) = mask(i+1,j+1)
    !  calculate the mask value
      weight(1) = (dble(i - xindex) + 1.0D0)*(dble(j-yindex) + 1.0D0)
      weight(2) = dble(xindex - i)*(dble(j-yindex) + 1.0D0)
      weight(3) = (dble(i-xindex) + 1.0D0)*dble(yindex - j)
      weight(4) = dble(xindex - i)*dble(yindex - j)
      mask_tmp  = dot_product(weight,MaskCom)
    
      if(mask_tmp .gt. 0.5) then
        mask_tmp= 1.0d0
      end if
    end if
  end subroutine cal_mask

end module mod_operator
