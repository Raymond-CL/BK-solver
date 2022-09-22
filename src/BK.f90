module BK

  use nrtype
  implicit none
  ! precision settings
  integer, parameter :: rp = sp
  ! i/o file names
  character(*), parameter :: ifname = 'input.dat'
  character(*), parameter :: ofname = 'BKtable.dat'
  ! kinematic boundary and data table
  integer(2) :: rn,yn
  real(rp) :: rmin,rmax, ymin,ymax, yh
  real(rp), allocatable, target, dimension(:,:) :: BKtable
  ! program options and flags
  integer(1) :: IniCnd,EvoMth,IntMth,EvoKer,RunCup,IntPol
  ! vegas settings
  integer :: ncall,itmax
  ! temporary pointers and parameters
  integer(2) :: vint_rind
  real(rp), dimension(:), pointer :: vint_in,vint_out
  real(rp), dimension(:), pointer :: k1,k2,k3,k4
  real(rp), dimension(:), pointer :: Nrt1,Nrt2,Nrt3

contains
  
  ! subroutine to set parameters (checked)
  subroutine setBK()
  implicit none
  integer, parameter :: ifunit = 11
  integer(1) :: err
  integer(2) :: ir
  real(rp) :: rt,r
  ! read from input.dat
  open(unit=ifunit,file=ifname,status='old')
  read(ifunit,*) rmin,rmax,rn
  read(ifunit,*) ymin,ymax,yn
  read(ifunit,*) IniCnd
  read(ifunit,*) EvoMth
  read(ifunit,*) IntMth
  read(ifunit,*) EvoKer
  read(ifunit,*) RunCup
  read(ifunit,*) IntPol
  read(ifunit,*) ncall,itmax
  close(ifunit)
  ! verify input data integraity
  if(rmin.le.0d0 .or. rmax.le.rmin .or. &
     ymin.lt.0d0 .or. ymax.le.ymin) call exit(2)
  ! initialize BK table
  yh = (ymax-ymin)/dble(yn)
  allocate(BKtable(1:rn,-1:yn),stat=err)
  if(err.ne.0) call exit(3)
  BKtable = 0d0
  ! fill initial condition data
  rt = log(rmax/rmin)/(rn-1)
  do ir = 1,rn
    r = rmin * exp(rt*(ir-1))
    BKtable(ir,-1) = r
    BKtable(ir, 0) = iniBK(r)
  enddo
  ! initialize temporary arrays for ODE
  if(EvoMth.eq.1) then
    allocate(k1(1:rn))
  elseif(EvoMth.eq.2) then
    allocate(k1(1:rn))
    allocate(k2(1:rn))
    allocate(Nrt1(1:rn))
  elseif(EvoMth.eq.3) then
    allocate(k1(1:rn))
    allocate(k2(1:rn))
    allocate(k3(1:rn))
    allocate(k4(1:rn))
    allocate(Nrt1(1:rn))
    allocate(Nrt2(1:rn))
    allocate(Nrt3(1:rn))
  endif
  end subroutine setBK



  ! BK initial condition (need update on other IC)
  function iniBK(r) result (ny0)
  implicit none
  real(rp), intent(in) :: r
  real(rp) :: ny0,Qs2,gam
  if(IniCnd.eq.1) then
    Qs2 = 0.24d0; gam = 1d0
    ny0 = 1d0 - exp( -(r*r*Qs2/4d0) )
  elseif(IniCnd.eq.2) then
    Qs2 = 0.15d0; gam = 1.13d0
    ny0 = 1d0 - exp( -(r*r*Qs2/4d0)**gam )
  elseif(IniCnd.eq.3) then
    ny0 = 1d0 - exp( -(r*r/4d0))
  else
    call exit(4)
  endif
  return
  end function iniBK



  ! prints program infomation (checked)
  subroutine BKinfo()
  write(*,*) "===========BK evolution solver==========="
  write(*,*) "given initial condition N(r,Y=0)"
  write(*,*) "generates table of dipole amplitude N(r,Y)"
  write(*,*) "by solving BK equation differential in Y"
  write(*,*) "-> reading input from: ",trim(ifname)
  write(*,*) "-> table range:"
  write(*,"(5x,es8.1,a6,es8.1,i5,a6)") rmin,"< r <",rmax,rn,"steps"
  write(*,"(5x,es8.1,a6,es8.1,i5,a6)") ymin,"< Y <",ymax,yn,"steps"
  if(IniCnd.eq.1) write(*,*) "-> initial condition: GBW"
  if(IniCnd.eq.2) write(*,*) "-> initial condition: MV"
  if(IniCnd.eq.3) write(*,*) "-> initial condition: user defined"
  if(EvoMth.eq.1) write(*,*) "-> evolution method: 1st-order Runge-Kutta"
  if(EvoMth.eq.2) write(*,*) "-> evolution method: 2nd-order Runge-Kutta"
  if(EvoMth.eq.3) write(*,*) "-> evolution method: 4th-order Runge-Kutta"
  if(IntMth.eq.1) write(*,*) "-> integrating as: dx dy"
  if(IntMth.eq.2) write(*,*) "-> integrating as: dphi r dr"
  if(IntMth.eq.3) write(*,*) "-> integrating as: dtheta r02 dr02"
  if(EvoKer.eq.1) write(*,*) "-> kernel: LO prescription"
  if(EvoKer.eq.2) write(*,*) "-> kernel: Balitsky prescription"
  if(EvoKer.eq.3) write(*,*) "-> kernel: Kovchegov-Weigert prescription"
  if(RunCup.eq.1) write(*,*) "-> alpha_s: fixed"
  if(RunCup.eq.2) write(*,*) "-> alpha_s: runs with r"
  if(RunCup.eq.3) write(*,*) "-> alpha_s: user defined"
  if(IntPol.eq.1) write(*,*) "-> Interpolation: linear in N(r) vs r"
  if(IntPol.eq.2) write(*,*) "-> Interpolation: linear in ln(N(r)) vs ln(r)"
  if(IntPol.eq.3) write(*,*) "-> Interpolation: cubic spline"
  write(*,"(a9,i12,a7,i3,a11)") "-> vegas",ncall,"calls,",itmax,"iterations"
  write(*,*) "-> output table result to: ",trim(ofname)
  write(*,*) "========================================="
  end subroutine BKinfo



  ! running evolution (Runge-Kutta)
  subroutine runBK()
  implicit none
  integer(2) :: iy
  do iy = 1,yn
  write(*,*) "solving: Y =", iy*(ymax/yn)
  if(EvoMth.eq.1) then
    vint_in => BKtable(:,iy-1)
    vint_out => k1
    call vint()
    BKtable(:,iy) = BKtable(:,iy-1) + k1(:)*yh
  elseif(EvoMth.eq.2) then
    vint_in => BKtable(:,iy-1)
    vint_out => k1
    call vint()
    Nrt1 = BKtable(:,iy-1) + k1*yh/2d0
    vint_in => Nrt1
    vint_out => k2
    call vint()
    BKtable(:,iy) = BKtable(:,iy-1) + k2(:)*yh
  elseif(EvoMth.eq.3) then
    vint_in => BKtable(:,iy-1)
    vint_out => k1
    call vint()
    Nrt1 = BKtable(:,iy-1) + k1(:)*yh/2d0
    vint_in => Nrt1
    vint_out => k2
    call vint()
    Nrt2 = BKtable(:,iy-1) + k2(:)*yh/2d0
    vint_in =>  Nrt2
    vint_out => k3
    call vint()
    Nrt3 = BKtable(:,iy-1) + k3(:)*yh
    vint_in => Nrt3
    vint_out => k4
    call vint()
    BKtable(:,iy) = BKtable(:,iy-1) + (k1(:)+2d0*k2(:)+2d0*k3(:)+k4(:))*yh/6d0
  endif
  enddo
  end subroutine runBK



  ! main integration loop
  subroutine vint()
  use nrtype
  use nr, only : vegas
  implicit none
  integer(2) :: ir
  integer(1), parameter :: ndim = 2
  integer(i4b) :: nprn,init
  real(sp) :: avgi,chi2a,sd
  real(sp), dimension(2*ndim) :: region
  nprn = -1
  if(IntMth.eq.1) then
    region(1) = -50d0;  region(3) = +50d0
    region(2) = -50d0;  region(4) = +50d0
  elseif(IntMth.eq.2) then
    region(1) = 0d0;  region(3) = 2d0*PI
    region(2) = 0d0;  region(4) = +50d0
  elseif(IntMth.eq.3) then
    region(1) = 0d0;  region(3) = 2d0*PI
    region(2) = 0d0;  region(4) = +50d0
  endif
  do ir = 1,rn
    vint_rind = ir
    !if(vint_in(vint_rind).eq.1d0) then
    !  vint_out(ir) = 0d0
    !else
      init = -1
      call vegas(region(1:2*ndim),fxn,init,ncall,itmax,nprn,avgi,sd,chi2a)
      init = +1
      call vegas(region(1:2*ndim),fxn,init,ncall,itmax,nprn,avgi,sd,chi2a)
      !if(avgi.lt.0d0) then
      !  if(avgi.lt.-0.1d0) write(*,*) 'wow, too negative there',ir,vint_in(ir),avgi
      !  avgi = 0d0
      !endif
      vint_out(ir) = avgi
    !endif
  enddo
  end subroutine vint



  ! main integrand for vegas
  function fxn(pt,wgt)
  use nrtype
  implicit none
  real(sp), dimension(:), intent(in) :: pt
  real(sp), intent(in) :: wgt
  real(sp) :: fxn
  real(rp) :: pt1,pt2
  real(rp) :: r01,r02,r12,K,Nr01,Nr02,Nr12
  fxn = 0d0
  pt1 = pt(1)
  pt2 = pt(2)
  r01 = BKtable(vint_rind,-1)
  if(IntMth.eq.1) then
    r02 = sqrt((r01/2d0+pt1)**2 + pt2**2)
    r12 = sqrt((r01/2d0-pt1)**2 + pt2**2)
  elseif(IntMth.eq.2) then
    r02 = sqrt((r01/2d0+pt2*cos(pt1))**2 + (pt2*sin(pt1))**2)
    r12 = sqrt((r01/2d0-pt2*cos(pt1))**2 + (pt2*sin(pt1))**2)
  elseif(IntMth.eq.3) then
    r02 = pt2
    r12 = sqrt( (r01-pt2*cos(pt1))**2 + (pt2*sin(pt1))**2 )
  endif
  ! divergence handling
  !if(r02.le.epsilon(r02) .or. r12.le.epsilon(r12)) return
  K = ker(r01,r02,r12)
  Nr01 = intpolr(vint_in,r01)
  Nr02 = intpolr(vint_in,r02)
  Nr12 = intpolr(vint_in,r12)
  fxn = K * (Nr02 + Nr12 - Nr01 - Nr02*Nr12)
  if(IntMth.eq.2 .or. IntMth.eq.3) fxn = fxn * pt2
  if(isnan(fxn)) then
    write(*,*) "fxn NaN:"
    write(*,*) pt1,pt2
    write(*,*) r01,r02,r12
    write(*,*) Nr01,Nr02,Nr12
    write(*,*) fxn
  endif
  return
  end function fxn
  
  ! evolution kernel
  function ker(r01,r02,r12) result(res)
  implicit none
  real(rp), intent(in) :: r01,r02,r12
  real(rp) :: res
  integer(1), parameter :: Nc = 3
  if(r02.eq.0d0 .or. r12.eq.0d0) then
    res = 0d0; return
  endif
  if(EvoKer.eq.1) then
    res = r01*r01/r02/r02/r12/r12
    res = res * Nc*as(r01)/2d0/PI/PI
  elseif(EvoKer.eq.2) then  ! Balitsky prescription
    res = r01*r01/r02/r02/r12/r12
    res = res + (as(r02)/as(r12)-1d0)/r02/r02
    res = res + (as(r12)/as(r02)-1d0)/r12/r12
    res = res * Nc*as(r01)/2d0/PI/PI
  elseif(EvoKer.eq.3) then  ! KW prescription (not yet)
    res = r01*r01/r02/r02/r12/r12
    res = res * Nc*as(r01)/2d0/PI/PI
  endif
  return
  end function ker

  ! coupling in coordinate space
  function as(r) result(res)
  implicit none
  real(rp), intent(in) :: r
  real(rp) :: res
  if(RunCup.eq.1) then
    res = 0.2d0 ! fixed
  elseif(RunCup.eq.2) then
    res = 0.7d0 ! to be implemented
  elseif(RunCup.eq.3) then
    res = 0d0   ! to be implemented
  endif
  return
  end function as

  ! simple interpolation function
  function intpolr(list,r) result(res)
  implicit none
  real(rp), dimension(1:rn), intent(in) :: list
  real(rp), intent(in) :: r
  real(rp) :: res,rtemp,r0,r1,nr0,nr1
  integer(2) :: ri0,ri1
  ! verify r within range
  if(r.lt.rmin) then
    res = 0d0;  return
  elseif(r.gt.rmax) then
    res = 1d0;  return
  endif
  ! get corresponding lower upper index
  rtemp = log(r/rmin)/log(rmax/rmin)*(rn-1) + 1
  ri0 = floor(rtemp)
  ri1 = ceiling(rtemp)
  ! get exact result if U=L index
  if(ri0.eq.ri1) then
    res = list(ri0);  return
  endif
  ! get grid points
  r0 = BKtable(ri0,-1)
  r1 = BKtable(ri1,-1)
  nr0 = list(ri0)
  nr1 = list(ri1)
  ! interpolate with selected methods
  if(IntPol.eq.1) then
    res = (r-r0)/(r1-r0)*(nr1-nr0)+nr0
  elseif(IntPol.eq.2) then
    if(nr0.eq.0d0) then ! use linear when lower bound is zero
      res = (r-r0)/(r1-r0)*nr1; return
    endif
    res = nr0 * exp( log(r/r0) / log(r1/r0) * log(nr1/nr0) )
  elseif(IntPol.eq.3) then
    ! intpolr = spline(list,r)
    res = 0d0
  endif
  ! check bound
  !if(res.lt.0d0) then
    !write(*,*) "interpolation result error",res
    !res = 0d0
  !endif
  !if(res.gt.1d0) then
    !write(*,*) "interpolation result error",res
    !res = 1d0
  !endif
  ! check NaN (develop mode)
  if(isnan(res)) then
    write(*,*) "interpolation NaN error"
    write(*,*) rtemp,ri0,ri1
    write(*,*) r0,r,r1
    write(*,*) nr0,res,nr1
    call sleep(5)
  endif
  return
  end function intpolr

  ! subroutine to print data table
  subroutine printBK()
  implicit none
  integer(2) :: ir,iy
  integer(1), parameter :: ofunit = 12
  write(*,*) "========================================="
  write(*,*) "writing BK table to ",trim(ofname)
  open(unit=ofunit,file=trim(ofname),status='replace')
  do ir = 1,rn
    do iy =-1,yn
      write(ofunit,"(es15.6)",advance="no") BKtable(ir,iy)
    enddo
    write(ofunit,*)
  enddo
  close(ofunit)
  write(*,*) "========================================="
  end subroutine printBK

end module BK
