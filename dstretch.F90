      subroutine MULTIGRID(P1,RHS1,xind)

use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mg_vari 
implicit none

!      include 'header_duct'
!c-----------------------------------------------------------------------
!c     this is an example of a main program using mgd9v
!c-----------------------------------------------------------------------
!c
      integer xind,i,j,k

!c
!c     parameter( nm=   774, nxf= 65, nyf= 33 )
!c     parameter( nm=  2919, nxf=129, nyf= 65 )
!c      parameter( nm= 11304, nxf=257, nyf=129 )
!c
!c      double precision
!c     +      v(nxf*nyf), vc(nm), vb(nxf*nyf), vbc(nm),val(nxf,nyf),
!c     +      rhs(nxf*nyf), rhsc(nm), a(nxf*nyf*9), ac(9*nm),
!c     +      ldu(3*nxf*nyf), lduc(3*nm), work(nxf*12),
!c     +      wa(nxf*nyf), wac(nm), wb(nxf*nyf), wbc(nm), tol, resno,
!c     +      gyf(0:nyf+1),gzf(0:nxf+1),hx,hy,x,y,pi

!      double precision    hx,hy,x,y     
      real(r8)          ::  P1(NY*NZ), RHS1(NY*NZ),bv(4)

!c     user data statements,
!c
!c     data nxc,nyc,levels/5,3,5/
!c     data nxc,nyc,levels/5,3,6/
!c      data nxc,nyc,levels/5,3,7/
!c
!c      data maxit,istart,iprep/30,1,0/
!c      data iout/1,0,0,2,1,0/
!      data iout/0,0,0,0,0,0/
!c     
!      data iprep/0/

      iprep = 0 
      iout(:)  = 0
!c      data nout/6/
!c
!c     open(unit=nout,file='output')
!c
!c     problem set up
!c
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!c     boundary condition      
      bc(1) = 1
      bc(2) = 2
      bc(3) = 2
      bc(4) = 2

      bv(:) = 0.d0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc

      v(:,xind)   = P1
      rhs(:,xind) = RHS1


      IF ( INIT_FLAG ) THEN

!c      write(nout,1)
!c    1 format('1')
!c
      call siam( NZ,NY,bc,KX2P(xind),bv,a(:,xind),dz,dzf,dy,dyf)
!c      call bc_con(NZ,NY,bc,bv,rhs(:,xind))
!c
!c     solution of the linear system
!c
!c
      iprep = 0
      tol=1.0d-9
      ifail=111
      istart=1
      maxit=0
!c
!      open(40,file='flow_field.txt',form='formatted',status='unknown')

      
      call mgd9v(levels, nzc, nyc, NZ, NY, nm,              &
                iout, istart, iprep, maxit, tol,            &
                rhs(:,xind), rhsc(:,xind), a(:,xind),       &
                ac(:,xind), v(:,xind), vc(:,xind),          &
                vb(:,xind), vbc(:,xind), work(:,xind),      &
                ldu(:,xind), lduc(:,xind), wa(:,xind),      & 
                wac(:,xind), wb(:,xind), wbc(:,xind),       &
                resno, ndid, ifail)

    
!c      write(6,*) "*MG grids have been initialized FOR*kx=",XIND,
!c     +    KX2(xind)

      ELSE
      tol=1.0d-8
      ifail=1
      iprep=1
      istart=1
      maxit=200


      call mgd9v(levels, nzc, nyc, NZ, NY, nm,              &
                iout, istart, iprep, maxit, tol,           &
                rhs(:,xind), rhsc(:,xind), a(:,xind),      &
                ac(:,xind), v(:,xind), vc(:,xind),         &  
                vb(:,xind), vbc(:,xind), work(:,xind),     &
                ldu(:,xind), lduc(:,xind), wa(:,xind),     &
                wac(:,xind), wb(:,xind), wbc(:,xind),      &
                resno, ndid, ifail)

      P1 = v(:,xind)


      ENDIF


      return
      end

      subroutine siam( nx,ny,bc,KX2,bv,l,dz,dzf,dy,dyf)
!c-----------------------------------------------------------------------
use ntypes
use mpi_var, only:rank 
      implicit none
!c
!c     testproblem to be found in
!c     p.m. de zeeuw and e.j. van asselt
!c     the convergence rate of multi-level algorithms applied to the
!c     convection-diffusion equation
!c     siam j. sci. stat. comput. 6(2) april 1985, p499.
!c
!c     ( u  + u  ) + a u = f(x,y)
!c        xx   yy
!c
!c     version dd911009
!c-----------------------------------------------------------------------
      integer nx, ny, nout
      integer bc(4)
      real(r8) :: eps, l(nx,ny,9), f(nx,ny),bv(4),dy(0:ny+1), &
          dyf(0:ny+1),dz(0:nx+1), dzf(0:nx+1),pi 
      integer i, j, k,nxm
      real(r8)   :: x, y, hx, hy, mu, mux, muy, a, axy, ah, b, bxy, bh, dir, r, KX2
      parameter ( pi=3.141592654) 

  

      call zeros(nx*ny*9,l)
      call zeros(nx*ny  ,f)

      do 20 j = 1, ny-1
        do 10 i = 1, nx-1
          l(i+1,j  ,4)= -1.d0/(dz(i+1)*dzf(i+1))
          l(i+1,j+1,4)= -1.d0/(dz(i+1)*dzf(i+1))
          l(i  ,j  ,6)= -1.d0/(dz(i+1)*dzf(i))
          l(i  ,j+1,6)= -1.d0/(dz(i+1)*dzf(i))

        
          l(i  ,j  ,8)= -1.d0/(dy(j+1)*dyf(j))
          l(i+1,j  ,8)= -1.d0/(dy(j+1)*dyf(j))
          l(i  ,j+1,2)= -1.d0/(dy(j+1)*dyf(j+1))
          l(i+1,j+1,2)= -1.d0/(dy(j+1)*dyf(j+1))

   10   continue
   20 continue

      do 40 j = 1, ny
        do 30 i = 1, nx
          l(i,j,5)= -(l(i,j,2)+l(i,j,4)+l(i,j,6)+l(i,j,8))   
!          l(i,j,5)=2.d0/hx**2.0+2.d0/hy**2.0                    
          l(i,j,5) = 1.d0/(dz(i+1)*dzf(i))+1.d0/(dz(i)*dzf(i))  &
               + 1.d0/(dy(j)*dyf(j)) +  1.d0/(dy(j+1)*dyf(j))  &
               + KX2
   30   continue
   40 continue


      IF ( (KX2 .LE. 0.d0 ) .and. (rank .eq. 0) ) THEN
      write(6,*)'#################################################'
      write(6,*)'boundary condition'
      write(6,*) 'Bottom wall', bc(1), 'Right wall', bc(2),'Top wall', &
                bc(3),'Left wall', bc(4)
      write(6,*)'#################################################'    
      write(6,*)'boundary value', 'Bottom wall', bv(1), 'Right wall',  &
                bv(2),'Top wall',  bv(3),'Left wall', bv(4)
      ENDIF


!c     Boundary conditions
       do j = 1,ny
!c      Dirichlet/periodic left face
        if( bc(4).eq.1 ) then
         l(1,j,:)   =  0.0d0
         l(1,j,5)   =  1.d0
         f(1,j)     =  bv(4)
        else
         l(1,j,:)   =  0.0d0
         l(1,j,5)   = 1/dz(2)**2.0
         l(1,j,6)   = -1/dz(2)**2.0
         f(1,j)     = -bv(4)/dz(2)
        endif

!c      Dirichlet/periodic right face
        if( bc(2).eq.1 ) then
         l(nx,j,:) =  0.0d0
         l(nx,j,5) =  1.0d0
         f(nx,j)   =  bv(2)
        else
         l(nx,j,:) =  0.0d0
         l(nx,j,5)   = 1/dz(nx)**2.0
         l(nx,j,4)   = -1/dz(nx)**2.0
         f(nx,j)     = bv(2)/dz(nx)
        endif
       enddo

       do i = 1,nx
!C      Dirichlet/periodic top face
        if( bc(3).eq.1 ) then
         l(i,ny,:) =  0.0d0
         l(i,ny,5) =  1.0d0
         f(i,ny)   =  bv(3)
        else
         l(i,ny,:) =  0.0d0
         l(i,ny,5) =  1/dy(ny)**2.0
         l(i,ny,2) =  -1/dy(ny)**2.0 
         f(i,ny)   =  bv(3)/dy(ny)  
        endif   
       enddo
       do i = 1,nx
!c      Dirichlet/periodic bottom face
        if( bc(1).eq.1 ) then
         l(i,1,:)   =  0.0d0
         l(i,1,5)   =  1.0d0
         f(i,1)     =  bv(1)
        else
         l(i,1,:)   =  0.0d0
         l(i,1,5)   =  1/dy(2)**2.0
         l(i,1,8)   =  -1/dy(2)**2.0
         f(i,1)     =  -bv(1)/dy(2)
        endif
       enddo


      return
      end

      subroutine bc_con(nx,ny,bc,bv,f)
      implicit none

      integer i,j,k,nx, ny,bc(4)
      real*8 f(nx,ny),bv(4)


     
!c     Boundary conditions
       do j = 1,ny
!c      Dirichlet/periodic left face
        if( bc(4).eq.1 ) then
         f(1,j)     =  bv(4)
        else
        endif

!c      Dirichlet/periodic right face
        if( bc(2).eq.1 ) then
         f(nx,j)   =  bv(2)
        else
        endif
       enddo

       do i = 1,nx
!c      Dirichlet/periodic bottom face
        if( bc(1).eq.1 ) then
         f(i,1)     =  bv(1)
        else
        endif

!c      Dirichlet/periodic top face
        if( bc(3).eq.1 ) then
         f(i,ny)   =  bv(3)
        else
        endif
       enddo 


      return
      end
 
      double precision function mu(eps,ah)
      double precision eps, ah
!c
           if( ah.gt.eps )    then
             mu=    eps/(2.0d0*ah)
      else if( ah.lt.(-eps) ) then
             mu=1.0d0+eps/(2.0d0*ah)
      else
             mu=0.5d0
      end if
      return
      end
