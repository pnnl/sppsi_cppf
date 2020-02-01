! We assume pointers are type 8 integer size
#define FFTW_ADDRESS_KIND selected_int_kind(16)
! We assume C int is type 4 integers
#define C_INT_KIND selected_int_kind(8)
! We assume C communicators are type 4 integer
# define MPI_C_WORLD_KIND selected_int_kind(8)
! We assume C and Fortran Communicators are the same
#define MPI_Comm_f2c(X) X
! Maximum number of MPI ranks
#define MAXCPU 8096
! Default FFTW plan
#define FFTW_PATIENT 32
#define MY_FFTW_OPTS FFTW_PATIENT
! We only call the C interface to FFTW, since there is no MKL MPI FFTW interface

! ref(1): Zhao et al., IJP 80 (2016) 38-55
! ref(2): Ma et al., Acta Materialia 54 (2006) 2169-2179
! Alt+x:print-to-pdf
! compile with ./compile_fcc_dislocation1.pbs
! chmod 777 compile_fcc_dislocation1.pbs

program dislocation
! based on the fcc_dislocation1.f from Dr. Li

!use mpi

implicit none

integer nx,ny,nz,nx2,nx21,ny21,nz21,nxyz,nss0,n_orderp,n_iteration
integer n_load,n_rate,n_ratestep,npower,n_slip
integer iteration,iteration_plas,input_euler,nstep,npt1,n_ppt
integer mp,mp0,mpc,mc,mc0,mc1,mc2,n_grain,mp1,mpr
integer i,j,k,ii,m,nm,nm1,m0,m1,m2,m3,nprint,kstep,n_random,kt,kt1
integer interf,lamda,mpin,kstep1,n_restart,kstart,n,n1
integer ns,max_grain,kprnt,k1,j1,i1,n_elas,n_plas,n_elas1,n_plas1
integer i2,j2,k2

! real*8 = real(kind = 8), extend to 15 digits of precision, -2.2e-308 to 1.8e308
real*8  rr,rp,R0,f0

integer xp0,yp0,zp0,NL

parameter (nx=16,ny=16,nz=16,nxyz=nx*ny*nz)
parameter (nx2=nx+2,nx21=nx/2+1,ny21=ny/2+1,nz21=nz/2+1)
parameter (mp0=5) ! # of grains
parameter (nss0=12) ! # of dislocation slip systems 
parameter (n_iteration=90000)
real*8 a(5,5),b(5),inva(5,5) !linear equations

real*8 pi
real*8 dt_real,tau_pass

! normalized plane normal as a function of grainID (mp0) and slip system (nss0)
real*8 nsp(mp0,3,nss0)
! normalized slip direction as a function of grainID (mp0) and slip system (nss0)
real*8 bsp(mp0,3,nss0) !each grain
! normalized sense vector as a function of grainID (mp0) and slip system (nss0)
! t = b x n
real*8 tsp(mp0,3,nss0),burg0
! Normal of slip planes and Burgers vector
integer num_grain(nx,ny,nz)

! resolved shear stress as a function of slip system, x, y, z
real*8 msig(nss0,nx,ny,nz) !stresses associated slip system
! shear rate as a function of slip system, x, y, z
real*8 gammas_dot(nss0,nx,ny,nz) ! shear rate


real*8 temp(nx,ny,nz)
real*8 temp1(nx,ny,nz),temp2(nx,ny,nz),temp3(nx,ny,nz)
complex*16 tempk(nx21,ny,nz),temp1k(nx21,ny,nz)
complex*16 temp2k(nx21,ny,nz),temp3k(nx21,ny,nz)

real*8 xk(nx+1),yk(ny+1),zk(nz+1)  ! k-vector 
real*8 dx,dy,dz 

real*8 density_GNDs(nss0,nx,ny,nz),density_GNDet(nss0,nx,ny,nz)
real*8 density_SSD(nss0,nx,ny,nz),density_M(nss0,nx,ny,nz)

real*8 density_SSD_dot,density_SSD0
real*8 density_F,density_P

real*8 density_GND1(nx,ny,nz)
real*8 density_SSD1(nx,ny,nz),density_M1(nx,ny,nz)

real*8 density_SSD2,density_M2,density_GND2
real*8 density_SSD3(n_iteration),density_M3(n_iteration)
real*8 density_GND3(n_iteration)

real*8 Ecoe,soft
real*8 sin_nt,cos_nt,sin_nb,cos_nb

real*8 coe1,coe2, coe3,coe4, coe5,coe6, coe7,coe8
real*8 kB,Tem,d_dipole,Qbulk,Gmu,poisson

character*8  fileindex  !iteration No.
character*15 filename   !stress.xxxxxxxx 
character*19 filenamepgm !stress.xxxxxxxx.dat

! To calculate dislocation densities, the following inputs are needed
! read shear rate and resolved shear stress from DAMASK
open(unit=11,file='gammas_dot.dat')
do n=1,nss0
do k=1,nz
do j=1,ny
do i=1,nx
! read(11,110) gammas_dot(n,i,j,k),msig(n,i,j,k)
read(11,*) gammas_dot(n,i,j,k),msig(n,i,j,k)
enddo
enddo
enddo
enddo !do n=1,nss0
!110   format(2(1x,d12.5))
close(11)

open(unit=12,file='grain_num.dat')  !read output from DAMASK
do k=1,nz
do j=1,ny
do i=1,nx
! read(12,120) num_grain(i,j,k) !integer
read(12,*) num_grain(i,j,k) !integer
enddo
enddo
enddo
120   format(i5)
close(12)

open(unit=13,file='slip_systems.dat') !read output from DAMASK
do n_grain=1,mp0
 do i=1,3 !vector components
   do n_slip=1,nss0
!    read(13,130) nsp(n_grain,i,n_slip)
!&                  ,bsp(n_grain,i,n_slip)
!&                  ,tsp(n_grain,i,n_slip)
    read(13, *) nsp(n_grain,i,n_slip), bsp(n_grain,i,n_slip), tsp(n_grain,i,n_slip)
enddo
enddo
enddo

130   format(3(1x,d12.5))
close(13)

! print first few items for verification
!rite(*, *) "shear rate and shear stress"
!o i = 1, 16
!   write(*, *) gammas_dot(1, i, 1, 1), msig(1, i, 1, 1)
!nd do
!rite(*, *) "grain ID"
!o i = 1, 16
!   write(*, *) num_grain(i, 1, 1)
!nd do
!rite(*, *) "slip system"
!o i = 1, 12
!   write(*, *) nsp(1, 1, i), bsp(1, 1, i), tsp(1, 1, i)
!nd do




! codes below are commented together in a time

! The COMMON statement defines a block of main memory storage so that different program units can
! share the same data without using arguments.
! common with gradk, kkxyz and dxyz are the common block names
common/kxyz/xk,yk,zk
common/dxyz/dx,dy,dz


! interface begin
integer, allocatable :: seed(:)
intrinsic :: random_seed, random_number

integer(kind=FFTW_ADDRESS_KIND) ::  alloc_local,local_nz,
&                           local_nz_offset,i8i,i8j,i8k
integer(kind=MPI_C_WORLD_KIND) :: c_world
integer myproc,numproc,ierr
integer :: mpi_pos1(MAXCPU),mpi_pos2(MAXCPU),mpi_size1(MAXCPU),
&                                  mpi_size2(MAXCPU)
logical msg
common /parallel/ alloc_local,local_nz,local_nz_offset,myproc,
&              numproc,mpi_pos1,mpi_pos2,mpi_size1,mpi_size2,msg

interface
 integer(kind=FFTW_ADDRESS_KIND) function fftw_mpi_local_size_3d
&        (i1,i2,i3,w,i4,i5)  bind(C, name='fftw_mpi_local_size_3d')
 integer(kind=FFTW_ADDRESS_KIND), value :: i1,i2,i3
 integer(kind=MPI_C_WORLD_KIND), value  :: w
 integer(kind=FFTW_ADDRESS_KIND) :: i4,i5
 end function
 subroutine fftw_mpi_init() bind(C, name='fftw_mpi_init')
 end subroutine fftw_mpi_init
end interface

call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierr)
if(numproc .gt. MAXCPU) then
 write(*,*) 'Too many CPU'
 stop
endif
call MPI_COMM_RANK(MPI_COMM_WORLD, myproc , ierr)
if(myproc .eq. 0) then
 msg=.true.
else
 msg=.false.
endif
call fftw_mpi_init()

! We have to call C API for fftw MPI Stuff.
! Crazy hoops to jump throught, but it works.
i8i = nz
i8j = ny
i8k = nx21
c_world = MPI_Comm_f2c(MPI_COMM_WORLD)
alloc_local=fftw_mpi_local_size_3d(i8i,i8j,i8k,
&         c_world,local_nz,local_nz_offset)
ii = local_nz*nx*ny
call mpi_gather(ii,1,MPI_INTEGER,mpi_size1   ,1,MPI_INTEGER,0,
&                         MPI_COMM_WORLD,ierr)
ii = local_nz*nx21*ny
call mpi_gather(ii,1,MPI_INTEGER,mpi_size2   ,1,MPI_INTEGER,0,
&                         MPI_COMM_WORLD,ierr)
if(msg) then
 mpi_pos1(1) = 0
 mpi_pos2(1) = 0
 do i=2,numproc
   mpi_pos1(i) = mpi_pos1(i-1) + mpi_size1(i-1)
   mpi_pos2(i) = mpi_pos2(i-1) + mpi_size2(i-1)
 enddo
endif
!      call second(time1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
pi=asin(1.0d0)*2.0d0
kB=1.38d-23 !(JK^-1)

!open(unit=1,file='input.dat')      !input
!
!read(1,*) dt_real
!read(1,*) dx,dy,dz       !e_rate0=0.0008 works
!read(1,*) nstep,kprnt
!read(1,*) burg0,density_SSD0 !burg0=1.0
!read(1,*) coe1,coe2,coe3,coe4
!read(1,*) coe5,coe6,coe7,coe8
!read(1,*) Ecoe,soft !for nucleation & grain growth
!read(1,*) Qbulk,Tem
!close(1)

dt_real = 1.0
dx = 1.0
dy = 1.0
dz = 1.0
burg0 = 1.0
density_DDS0 = 0.0
coe1 = 0.18
coe2 = 5.0
coe3 = 5.0
coe4 = 8.0e6
coe5 = 15.0
coe6 = 5.0e12
coe7 = 1.0e-29
coe8 = 0.33
qbulk = 2.4e-19 ! in J



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc start from here

call gradk  !output xk(nx),yk(ny),zk(nz)

do k=1+local_nz_offset,local_nz+local_nz_offset
do j=1,ny
do i=1,nx

 do n=1,nss0  !nss0=9
   density_SSD(n,i,j,k)=density_SSD0
   density_GNDs(n,i,j,k)=0.0
   density_GNDet(n,i,j,k)=0.0
   density_M(n,i,j,k)=0.0
 enddo

enddo
enddo
enddo

c   calculate gradient of gammas_dot

do n=1,nss0
do k=1+local_nz_offset,local_nz+local_nz_offset
do j=1,ny
do i=1,nx
 temp(i,j,k)=gammas_dot(n,i,j,k)  !read output from DAMASK
enddo
enddo
enddo
call forward(temp,tempk,nx,ny,nz)

do k=1+local_nz_offset,local_nz+local_nz_offset
do j=1,ny
do i=1,nx21
 temp1k(i,j,k)=(0.0,1.0)*xk(i)*tempk(i,j,k) !gradient
 temp2k(i,j,k)=(0.0,1.0)*yk(j)*tempk(i,j,k)
 temp3k(i,j,k)=(0.0,1.0)*zk(k)*tempk(i,j,k)
enddo
enddo
enddo
call backward(temp1k,temp1,nx,ny,nz)
call backward(temp2k,temp2,nx,ny,nz)
call backward(temp3k,temp3,nx,ny,nz)

do k=1+local_nz_offset,local_nz+local_nz_offset
do j=1,ny
do i=1,nx
 n_grain=num_grain(i,j,k) !for orientation
 density_GNDs(n,i,j,k)= density_GNDs(n,i,j,k) !eq(7) of Ref(1)
&            -(temp1(i,j,k)*tsp(n_grain,1,n)    !n--slip system
&             +temp2(i,j,k)*tsp(n_grain,2,n)    !not the rate
&             +temp3(i,j,k)*tsp(n_grain,3,n))/burg0*dt_real

 density_GNDet(n,i,j,k)= density_GNDet(n,i,j,k)
&            +(temp1(i,j,k)*bsp(n_grain,1,n)
&             +temp2(i,j,k)*bsp(n_grain,2,n)
&             +temp3(i,j,k)*bsp(n_grain,3,n))/burg0*dt_real
enddo
enddo
enddo
enddo !do n=1,nss0

c        goto 1330 !not going to calculate density_M & density_SSD

do n=1,nss0  !see eqs.(4) and (5) of ref(1)
do k=1+local_nz_offset,local_nz+local_nz_offset
do j=1,ny
do i=1,nx

 n_grain=num_grain(i,j,k)
 density_F=0.0
 density_P=0.0  ! tsp() is sense vector of slip system n1 on n_grain
 do n1=1,nss0  !summation
   rr=sqrt(tsp(n_grain,1,n1)**2+tsp(n_grain,2,n1)**2
&                                 +tsp(n_grain,3,n1)**2) !rr should be =1
   cos_nt=(nsp(n_grain,1,n)*tsp(n_grain,1,n1)
&            +nsp(n_grain,2,n)*tsp(n_grain,2,n1)
&            +nsp(n_grain,3,n)*tsp(n_grain,3,n1))/rr

   sin_nt=sqrt(abs(1.0d0-(cos_nt)**2))
   cos_nb=nsp(n_grain,1,n)*bsp(n_grain,1,n1)
! bsp() is Burgers vector (unit vector=slip direction)
&           +nsp(n_grain,2,n)*bsp(n_grain,2,n1)
&           +nsp(n_grain,3,n)*bsp(n_grain,3,n1)

   sin_nb=sqrt(abs(1.0d0-(cos_nb)**2))
   density_F=density_F            !eq(4) of Ref(2)
&                        +abs(density_SSD(n1,i,j,k)*cos_nt)
&                        +abs(density_GNDs(n1,i,j,k)*cos_nb)
&                        +abs(density_GNDet(n1,i,j,k)*cos_nt)
   density_P=density_P            !eq(5) of Ref(2)
&                        +abs(density_SSD(n1,i,j,k)*sin_nt)
&                        +abs(density_GNDs(n1,i,j,k)*sin_nb)
&                        +abs(density_GNDet(n1,i,j,k)*sin_nt)

 enddo !n1=1,nss0 !summation

 density_M(n,i,j,k)= 2.0*kB/(coe1*coe2*coe3*Gmu*burg0**3)
&                  *Tem*sqrt(density_P)*sqrt(density_F) !eq(13) of ref(2)

 if(abs(msig(n,i,j,k)).gt.tau_pass) then
                   ! critical shear stress on slip system n
  d_dipole=sqrt(3.0)*Gmu*burg0
&                 /(16.0*pi*(1.0-poisson)*msig(n,i,j,k))
      !see p6 of ref(1) between eqs(6) & (7)
 else
  d_dipole=0.0
 endif

      !see p2173 equ(18) Ma et al. Acta Mater. 54(2006)2169, Ref(2)
 density_SSD_dot= (coe4*sqrt(density_F)
&                   +coe6*d_dipole*density_M(n,i,j,k)
&                   -coe5*density_SSD(n,i,j,k))
&                                        *abs(gammas_dot(n,i,j,k))
                    ! why need abs()??  HU
&     -coe7*exp(-Qbulk/(kB*Tem))*abs(msig(n,i,j,k))/(kB*Tem)
&          *density_SSD(n,i,j,k)**2
&           *abs(gammas_dot(n,i,j,k))**coe8
      ! should be the von Mises equivalent shear rate   HU

c         density_SSD00(n,i,j,k)= density_SSD(n,i,j,k)
c     &                        +density_SSD_dot*dt_real
     !no need of density_SSD00(n,i,j,k)
 density_SSD(n,i,j,k)= density_SSD(n,i,j,k)
&                        +density_SSD_dot*dt_real

enddo
enddo
enddo

enddo !n=1,nss0

1330    continue !without calculating density_M & density_SSD

c      updating tau value
do k=1+local_nz_offset,local_nz+local_nz_offset
do j=1,ny
do i=1,nx
 density_M1(i,j,k)=0.0 !at each grid, sum of all slip system
 density_SSD1(i,j,k)=0.0
 density_GND1(i,j,k)=0.0

 do n=1,nss0

   density_M1(i,j,k)=density_M1(i,j,k)+density_M(n,i,j,k)
   density_SSD1(i,j,k)=density_SSD1(i,j,k)+density_SSD(n,i,j,k)
   density_GND1(i,j,k)=density_GND1(i,j,k)
&                        +sqrt(density_GNDs(n,i,j,k)**2
&                             +density_GNDet(n,i,j,k)**2)
 enddo

enddo
enddo
enddo

c  ----------------dislocation density end-----------------------------
c    output results
   call gather_real(density_M1)
   call gather_real(density_SSD1)
   call gather_real(density_GND1)

if(msg) then
c       write(fileindex,'(I8)') kt

c       do ispace=1,len(fileindex)
c          if (fileindex(ispace:ispace).eq." ") then
c             fileindex(ispace:ispace)="0"
c          endif
c       end do

c       filename='orderp.'//fileindex
c       filenamepgm=filename//'.dat'
c            open(unit=11,file=filenamepgm)

   open(unit=10,file='dis_M.dat')
    do i=1,nx
    do j=1,ny
    do k=1,nz
     write(10,1108) density_M1(i,j,k)
    enddo
    enddo
    enddo
    close(10)

1108       format(d12.5)

   open(unit=10,file='dis_SSD.dat')
    do i=1,nx
    do j=1,ny
    do k=1,nz
     write(10,1108) density_SSD1(i,j,k)
    enddo
    enddo
    enddo
    close(10)

   open(unit=10,file='dis_GND.dat')
    do i=1,nx
    do j=1,ny
    do k=1,nz
     write(10,1108) density_GND1(i,j,k)
    enddo
    enddo
    enddo
    close(10)


 density_M2=0.0  !average
 density_SSD2=0.0
 density_GND2=0.0
do k=1,nz
do j=1,ny
do i=1,nx

  density_M2=density_M2+density_M1(i,j,k)/nxyz
  density_SSD2=density_SSD2+density_SSD1(i,j,k)/nxyz
  density_GND2=density_GND2+density_GND1(i,j,k)/nxyz

enddo
enddo
enddo

 density_M3(kt)=density_M2  !stored for each kt step
 density_SSD3(kt)=density_SSD2
 density_GND3(kt)=density_GND2
   open(unit=10,file='dis_total.dat')
     write(10,1110) density_M2,density_SSD2,density_GND2
1110        format(3(1x,d12.5))
   close(10)

 endif !if(msg)

call MPI_FINALIZE(ierr)
stop

contains

C for FFTW
!     ========================================================
subroutine gradk
implicit none
integer nx,ny,nz,nx2,nx21,ny21,nz21,i,j,k,ny2,nz2,nxyz

parameter (nx=64,ny=64,nz=64,nxyz=nx*ny*nz)
parameter (nx2=nx+2,nx21=nx/2+1,ny21=ny/2+1,nz21=nz/2+1)

!Real Variables
real*8 fnx,fny,fnz
real*8 fk1,fk2,fk3
real*8 rr,d0
real*8 dx,dy,dz
real*8 pi,twopi
!Arrays
real*8 xk(nx+1),yk(ny+1),zk(nz+1)

!Common variables

common /dxyz/dx,dy,dz
common /kxyz/xk,yk,zk


pi=4.0*atan(1.0)
twopi=2.0*pi

!     ***
ny2=ny+2
nz2=nz+2

fnx=twopi/(nx*dx)
fny=twopi/(ny*dy)
fnz=twopi/(nz*dz)

do i=1,nx21
 fk1=float(i-1)*fnx
  xk(i)=fk1
  xk(nx2-i)=-fk1
enddo

do j=1,ny21
 fk2=float(j-1)*fny
  yk(j)=fk2
  yk(ny2-j)=-fk2
enddo

do k=1,nz21
 fk3=float(k-1)*fnz
  zk(k)=fk3
  zk(nz2-k)=-fk3
enddo

return
end subroutine gradk

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2@cc
!c     A code to do FFT using FFTW 3.0.1
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!      program  newFFT
! fft2Df
subroutine forward(r,c,n1,n2,n3)

!      parameter (n1=256,n2=1,n3=256)  !,n1=nx,n2=ny,n3=nz)
use mpi
    implicit none
!            include 'fftw3.f'
integer n1,n2,n3
logical ifirst
save ifirst,planfw
real*8 r(n1,n2,n3)
complex*16 c(n1/2+1,n2,n3)
real*8 scale
integer(kind=FFTW_ADDRESS_KIND) ::   planfw
data ifirst/.true./
integer loc
integer i,j,k

integer(kind=FFTW_ADDRESS_KIND) ::  alloc_local,local_nz,
&                                local_nz_offset
integer myproc,numproc,ierr
integer :: mpi_pos1(MAXCPU),mpi_pos2(MAXCPU),mpi_size1(MAXCPU),
&                                       mpi_size2(MAXCPU)
logical msg
common /parallel/ alloc_local,local_nz,local_nz_offset,myproc,
&             numproc,mpi_pos1,mpi_pos2,mpi_size1,mpi_size2,msg
real*8, allocatable, save :: r_loc(:)
complex*16, allocatable, save :: c_loc(:)
integer(kind=FFTW_ADDRESS_KIND) ::  i8i,i8j,i8k
integer(kind=MPI_C_WORLD_KIND) :: c_world

interface
integer(kind=FFTW_ADDRESS_KIND) function fftw_mpi_plan_dft_r2c_3d
&                            (n0,n1,n2,in,out,comm,flags)
&           bind(C, name='fftw_mpi_plan_dft_r2c_3d')
integer(kind=FFTW_ADDRESS_KIND), value :: n0
integer(kind=FFTW_ADDRESS_KIND), value :: n1
integer(kind=FFTW_ADDRESS_KIND), value :: n2
real*8, dimension(*), intent(inout) :: in
complex*16, dimension(*), intent(inout) :: out
integer(kind=MPI_C_WORLD_KIND), value :: comm
integer(kind=C_INT_KIND), value :: flags
end function fftw_mpi_plan_dft_r2c_3d
subroutine fftw_execute(plan) bind(C,name='fftw_execute')
integer(kind=FFTW_ADDRESS_KIND), value:: plan
end subroutine
end interface

if (ifirst) then
allocate(r_loc(2*alloc_local),c_loc(alloc_local))
i8i = n3
i8j = n2
i8k = n1
c_world = MPI_Comm_f2c(MPI_COMM_WORLD)
planfw = fftw_mpi_plan_dft_r2c_3d(i8i,i8j,i8k,
&                 r_loc,c_loc,c_world,MY_FFTW_OPTS)
ifirst=.false.
end if

loc = 1
do k=1+local_nz_offset,local_nz+local_nz_offset
do j=1,n2
do i=1,n1
 r_loc(loc)=r(i,j,k)  ! Copy to local copy
 loc=loc+1
enddo
loc=loc+2 ! MPI padding
enddo
enddo

call fftw_execute(planfw)

loc = 1
do k=1+local_nz_offset,local_nz+local_nz_offset
do j=1,n2
do i=1,n1/2 + 1
 c(i,j,k)=c_loc(loc)  ! Copy from local copy
 loc=loc+1
enddo
enddo
enddo

return
end subroutine forward


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      program  newFFTb
!      subroutine fft2Db(c,r,n1,n2,n3)
subroutine backward(c,r,n1,n2,n3)
!      parameter (n1=256,n2=1,n3=256)  !,n1=nx,n2=ny,n3=nz)
use mpi
    implicit none
!            include 'fftw3.f'
integer n1,n2,n3,i,j,k,loc
logical ifirst
save ifirst,planbw
real*8 r(n1,n2,n3),scale
complex*16 c(n1/2+1,n2,n3),tempk(n1/2+1,n2,n3)
integer(kind=FFTW_ADDRESS_KIND) ::   planbw
data ifirst/.true./

integer(kind=FFTW_ADDRESS_KIND) ::  alloc_local,local_nz,
&                             local_nz_offset
integer myproc,numproc,ierr
integer :: mpi_pos1(MAXCPU),mpi_pos2(MAXCPU),mpi_size1(MAXCPU),
&                       mpi_size2(MAXCPU)
logical msg
common /parallel/ alloc_local,local_nz,local_nz_offset,myproc,
&               numproc,mpi_pos1,mpi_pos2,mpi_size1,mpi_size2,msg
real*8, allocatable, save :: r_loc(:)
complex*16, allocatable, save :: c_loc(:)
integer(kind=FFTW_ADDRESS_KIND) ::  i8i,i8j,i8k
integer(kind=MPI_C_WORLD_KIND) c_world

interface
integer(kind=FFTW_ADDRESS_KIND) function fftw_mpi_plan_dft_c2r_3d
&                                     (n0,n1,n2,in,out,comm,flags)
&            bind(C, name='fftw_mpi_plan_dft_c2r_3d')
integer(kind=FFTW_ADDRESS_KIND), value :: n0
integer(kind=FFTW_ADDRESS_KIND), value :: n1
integer(kind=FFTW_ADDRESS_KIND), value :: n2
complex*16, dimension(*), intent(inout) :: in
real*8, dimension(*), intent(inout) :: out
integer(kind=MPI_C_WORLD_KIND), value :: comm
integer(kind=C_INT_KIND), value :: flags
end function fftw_mpi_plan_dft_c2r_3d
subroutine fftw_execute(plan) bind(C,name='fftw_execute')
integer(kind=FFTW_ADDRESS_KIND), value:: plan
end subroutine
end interface

if (ifirst) then
allocate(r_loc(2*alloc_local),c_loc(alloc_local))
i8i = n3
i8j = n2
i8k = n1
c_world = MPI_Comm_f2c(MPI_COMM_WORLD)
planbw = fftw_mpi_plan_dft_c2r_3d(i8i,i8j,i8k,
&                      c_loc,r_loc,c_world,MY_FFTW_OPTS)
ifirst=.false.
end if

loc = 1
do k=1+local_nz_offset,local_nz+local_nz_offset
do j=1,n2
do i=1,n1/2 + 1
 c_loc(loc) = c(i,j,k) ! Copy to local copy
 loc=loc+1
enddo
enddo
enddo

call fftw_execute(planbw)

scale = 1.d0/(n1*n2*n3)

loc = 1
do k=1+local_nz_offset,local_nz+local_nz_offset
do j=1,n2
do i=1,n1
 r(i,j,k) = r_loc(loc)*scale ! Copy from local copy
 loc=loc+1
enddo
loc=loc+2 ! MPI padding
enddo
enddo

return
end subroutine backward


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine gather_real(stuff)
use mpi
implicit none
integer nx,ny,nz
parameter (nx=64,ny=16,nz=32)
real*8, intent(inout) :: stuff(nx,ny,nz)
integer length
integer(kind=FFTW_ADDRESS_KIND) ::  alloc_local,local_nz,
&                               local_nz_offset
integer myproc,numproc,ierr
integer :: mpi_pos1(MAXCPU),mpi_pos2(MAXCPU),mpi_size1(MAXCPU),
&                         mpi_size2(MAXCPU)
logical msg
common /parallel/ alloc_local,local_nz,local_nz_offset,myproc,
&               numproc,mpi_pos1,mpi_pos2,mpi_size1,mpi_size2,msg

length = nx*ny*local_nz

IF(MSG) THEN
   call mpi_gatherv(MPI_IN_PLACE,
&          length,MPI_DOUBLE_PRECISION,stuff,
&          mpi_size1,mpi_pos1,
&          MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
ELSE
   call mpi_gatherv(stuff(1,1,1+local_nz_offset),
&          length,MPI_DOUBLE_PRECISION,stuff,
&          mpi_size1,mpi_pos1,
&          MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
ENDIF
end subroutine gather_real

end program dislocation
