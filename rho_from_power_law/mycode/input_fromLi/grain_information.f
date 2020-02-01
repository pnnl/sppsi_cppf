! modified from pa_VP_solution3D64_newb_Li6b.f, subroutine read_order1
! 2020/01/16/
! gfortran -o grain_information.out grain_information.f

      program grains
      implicit none
      integer nx,ny,nz,nxyz
      integer ns,nss0,mp0
      integer i,j,k,m1,m2,n,n1,nn,number

      parameter (nx=64,ny=64,nz=64,nxyz=nx*ny*nz) 
      parameter (nss0=12) ! # of dislocation slip systems 
      parameter (mp0=40)

      real*8 nsp0(3,nss0),bsp0(3,nss0),tsp0(3,nss0) !cross product
      real*8 p1,rr
      integer num(mp0),num_grain(nx,ny,nz)

      real*8 phi0(mp0),theta0(mp0),shi0(mp0)
      real*8 trans(mp0,3,3),trans0(3,3)
      
      real*8 nsp(mp0,3,nss0)
      real*8 bsp(mp0,3,nss0)
      real*8 tsp(mp0,3,nss0)

      open(unit=1,file='slipsystem.dat')      !input, slip system

c   slip systems
      do n=1,nss0  !nss0=12
       read(1,*) nsp0(1,n),nsp0(2,n),nsp0(3,n)
     &          ,bsp0(1,n),bsp0(2,n),bsp0(3,n)
      enddo !nss0=12

      close(1)

      open(unit=2,file='record.dat') !output
      write(2,*)' slip systems:'

!to change into unit vectors
      do n=1,nss0  !nss0=12
       rr=sqrt(nsp0(1,n)**2+nsp0(2,n)**2+nsp0(3,n)**2)
       nsp0(1,n)=nsp0(1,n)/rr  !unit normalization
       nsp0(2,n)=nsp0(2,n)/rr  !given normal vector
       nsp0(3,n)=nsp0(3,n)/rr

       rr=sqrt(bsp0(1,n)**2+bsp0(2,n)**2+bsp0(3,n)**2)
       bsp0(1,n)=bsp0(1,n)/rr   !given burger-vector
       bsp0(2,n)=bsp0(2,n)/rr
       bsp0(3,n)=bsp0(3,n)/rr
c    cross product
       tsp0(1,n)= bsp0(3,n)*nsp0(2,n) - bsp0(2,n)*nsp0(3,n)
       tsp0(2,n)=-bsp0(3,n)*nsp0(1,n) + bsp0(1,n)*nsp0(3,n)
       tsp0(3,n)= bsp0(2,n)*nsp0(1,n) - bsp0(1,n)*nsp0(2,n)
       write(2,203) n,nsp0(1,n),nsp0(2,n),nsp0(3,n)
     &               ,bsp0(1,n),bsp0(2,n),bsp0(3,n)
     &               ,tsp0(1,n),tsp0(2,n),tsp0(3,n)
      enddo !nss0=12

 203  format(i2,9(1x,e11.4))

     
      open (unit=101,file='c_ordp.dat') !input
      open (unit=102,file='grain_num.dat') !output

 1015  format(1x,e11.4,I4)
 1017  format(I5)

       do k=1,nz
       do j=1,ny  !see PolyCrytFrom64newCP1Li.F,unit=11
       do i=1,nx  !they are the ouput order: k,j,i
c  spatial each point
        read(101,1015) p1,n  !n-order parameter/grain No.
         num_grain(i,j,k)=n !store the order parameter number at the given point

        write(102,1017)  num_grain(i,j,k)  

       enddo
       enddo
       enddo 

c  output & check if all is right
      close(101)
      close(102)
c had grain number of each point

       open(unit=55,file='euler_angle.dat')
        do n=1,mp0
         read(55,550) nn,phi0(nn),theta0(nn),shi0(nn)
        enddo  !n=1,mp0
 550   format(i3,3(1x,f7.2))
       close(55)


       do n1=1,mp0  !grain orientation number

        call rotation(phi0(n1),theta0(n1),shi0(n1)
     &                 ,trans0) !angles in degree

        do m1=1,3
        do m2=1,3
         trans(n1,m1,m2)=trans0(m1,m2)
        enddo  !m2=1,3
        enddo  !m1=1,3

c  output transformed matrix
        do ns=1,nss0 
         do m1=1,3
          bsp(n1,m1,ns)=0.0
          nsp(n1,m1,ns)=0.0
          tsp(n1,m1,ns)=0.0
         do m2=1,3
           bsp(n1,m1,ns)=bsp(n1,m1,ns)+trans(n1,m1,m2)*bsp0(m2,ns)
           nsp(n1,m1,ns)=nsp(n1,m1,ns)+trans(n1,m1,m2)*nsp0(m2,ns)
           tsp(n1,m1,ns)=tsp(n1,m1,ns)+trans(n1,m1,m2)*tsp0(m2,ns)
         enddo !m2=1,3
         enddo !m1=1,3
        enddo  !ns=1,nss0

       enddo !do n=1,mp0  !grain orientation number

!   output for checking
      do n=1,mp0
       write(2,*) 
       write(2,108) n,phi0(n),theta0(n),shi0(n)
 108   format(' n_orientation=',i3,' phi0=',f7.2,' theta0=',f7.2
     &                           ,' shi0=',f7.2)
       write(2,*) 
       write(2,*) ' transform matrix'
       do m1=1,3
         write(2,1012) (trans(n,m1,m2),m2=1,3)
       enddo

       write(2,*) 
       write(2,*) ' slip plane normal and Burger vector'
       do ns=1,nss0
         write(2,1010) (nsp(n,m1,ns),m1=1,3),(bsp(n,m1,ns),m1=1,3)
     &                ,(tsp(n,m1,ns),m1=1,3)
       enddo  !do ns=1,nss0
      enddo  !n=1,mp0
 1010   format(9(1x,e11.4))
 1012   format(3(1x,e11.4))
       write(2,*) 

c real stuff of the sub
       do n=1,mp0
        num(n)=0 !points of orientation n
       enddo

      do k=1,nz
      do j=1,ny 
      do i=1,nx  !they are the ouput order: k,j,i
c  spatial each point
        n=num_grain(i,j,k)
        num(n)=num(n)+1

      enddo !i=1,nx
      enddo !j=1,ny
      enddo !k=1,nz
  
       write(2,*)
       number=0
       do n=1,mp0
        write(2,128) n,num(n) !orientations & sizes
         number=number+num(n)
       enddo

       write(2,*) '  total points=',number,nx*ny*nz

 128   format(' grain No.=',i3,' points=',i8)

 
          close(2)

       end 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rotation(phi0,theta0,shi0,trans) !angles in deg.
      implicit none

      real*8 trans(3,3)
      real*8 pi,phi0,theta0,shi0,phi,theta,shi

      real*8 tra11,tra12,tra13,tra21,tra22,tra23,tra31,tra32,tra33
      real*8 tr11,tr12,tr13
      real*8 tr21,tr22,tr23
      real*8 tr31,tr32,tr33

       pi=dacos(-1.d0)

          phi=phi0*pi/180.
          theta=theta0*pi/180.
          shi=shi0*pi/180.

c   Creating the tranformation matrix from local to global and global to local
c      P(i)(local)=tra(i,j)*P(j)(global)
c      P(i)(local)=tra(i,j)*P(j)(global) !P=vector
c   'tra' for global to local  
c   'tr' is local to global which is just the transpose matrix of 'tra'
         tra11=cos(phi)*cos(shi)
     %         -cos(theta)*sin(phi)*sin(shi) 

         tra12=cos(shi)*sin(phi)
     %         +cos(theta)*cos(phi)*sin(shi)
         tra13=sin(theta)*sin(shi)
         tra21=-sin(phi)*cos(shi)*cos(theta)
     %     -cos(phi)*sin(shi)
         tra22=cos(shi)*cos(phi)*cos(theta)
     %     -sin(phi)*sin(shi)
         tra23=sin(theta)*cos(shi)
         tra31=sin(phi)*sin(theta)
         tra32=-cos(phi)*sin(theta)
         tra33=cos(theta) ! end of global to local transformation matrix
c they are the components of TR of Eq.(2.4) in Samrat's thesis

c       write(33,1051) tra11,tra12,tra13
c       write(33,1051) tra21,tra22,tra23
c       write(33,1051) tra31,tra32,tra33
c       write(33,*)

       tr11=tra11  !transpose of tra
       tr12=tra21  !so,P(i)(global)=tr(i,j)*P(j)(local)
       tr13=tra31
       tr21=tra12
       tr22=tra22
       tr23=tra32
       tr31=tra13
       tr32=tra23
       tr33=tra33! end of local to global tranformation matrix

! input into a matrix
       trans(1,1)=tr11  !so,P(i)(global)=tr(i,j)*P(j)(local)
       trans(1,2)=tr12
       trans(1,3)=tr13

       trans(2,1)=tr21
       trans(2,2)=tr22
       trans(2,3)=tr23

       trans(3,1)=tr31
       trans(3,2)=tr32
       trans(3,3)=tr33 ! end of local to global tranformation matri

       return 

       end !subroutine stiffT_sub
cccccccccccccccccccccccccccccccccccccccccccccccccc
