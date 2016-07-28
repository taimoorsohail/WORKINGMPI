
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS_DUCT(FINAL)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!      INCLUDE 'header_duct'
 
use ntypes
use Domain
use Grid
use Fft_var, only :  pi, CIKX
use TIME_STEP_VAR
use run_variable
use variable_stat
use mpi_var

implicit none

      LOGICAL FINAL, TKE_BUDGET, MEAN_KE_BUDGET,MEAN_BACKGROUND,SAVE_3D, &
              ZX_PLANE, XY_PLANE, ZX_PLANE_2,XY_PLANE_comb, ZX_PLANE_comb, &
                ZX_PLANE_comb_2, PLANE_3D, XY_PLANE_comb_2    
      CHARACTER*35 FNAME
      integer i,j,k,n
      real(r8) :: uc, ubulk

       OPEN (11,file='signal_plane_data',form='formatted',status='old')
        READ(11,*)
        READ(11,*) XY_PLANE_comb,XY_PLANE_comb_2, ZX_PLANE_comb,ZX_PLANE_comb_2, PLANE_3D
       close(11)

    
       XY_PLANE        = .false.
       ZX_PLANE        = .false.
       ZX_PLANE_2      = .false.
!      PLANE_3D        = .true.
!      ZX_PLANE_comb_2 = .true.
!      XY_PLANE_comb   = .true.
!      ZX_PLANE_comb   = .true.
       SAVE_3D         = .false.
       TKE_BUDGET      = .false.
       MEAN_KE_BUDGET  = .false.
       MEAN_BACKGROUND = .false.

 IF ( rank .eq. 0) THEN
  WRITE(6,*) 'Saving flow statistics.'
!  write(6,*) 'Allocate all the tmp arrays'
  IF (SAVE_3D) THEN
   write(6,*) 'Writing 3D data field'
  ENDIF
ENDIF
!  allocating stat variables      
!      call allocate_temps
!C Apply Boundary conditions to velocity field
!       IF (USE_MPI) THEN
!         CALL APPLY_BC_VEL_MPI
!       ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
        CALL APPLY_BC_VEL_LEFT
        CALL APPLY_BC_VEL_RIGHT
!      END IF


      IF (RANK .eq. 0) THEN
      
! Compute and write out bulk velocity
! Integrat the instantaneous mean profile numerically at GY points
      UBULK= 0.0d0
      area = 0.0d0
      DO J=2,NY
       DO K=2,NZ
        area = area + DZ(k)*DY(j)
        UBULK=UBULK+0.25*(dble(CU3X(0,K,J))+dble(CU3X(0,K-1,J)) +  &
                dble(CU3X(0,K-1,J-1)) + dble(CU3X(0,K,J-1)) )*DY(j)*DZ(K)
       ENDDO
      ENDDO
      UBULK=UBULK/area      
! Write out UBULK
      write(6,*) 'UBULK: ',UBULK
     ENDIF

! Save CUi
      do k=0,NZ+1
        do i=0,NX2P
          do j=0,NY+1
            CR1X(i,k,j)=CU1X(i,k,j)
            CR2X(i,k,j)=CU2X(i,k,j)
            CR3X(i,k,j)=CU3X(i,k,j)
!   THIS STEPs ARE REQURIED WHEN DERIVATIVES w.r.t X IS REQURIED LATER
            CF1X(i,k,j)=CU1X(i,k,j)
            CF2X(i,k,j)=CU2X(i,k,j)
            CF3X(i,k,j)=CU3X(i,k,j)
          end do
        end do
      end do

    !  write(6,*)'I am here flow_stst_0', rank

!    TRANSFERRING MEAN TO ALL NODES 

      do k=0,NZ+1
          do j=0,NY+1
            p_mean(k,j)=CPX(0,k,j)
            CALL MPI_BCAST_COMPLEX(p_mean(k,j),1,1)       
          enddo
      enddo

    

      do k=0,NZV-1
          do j=0,NY+1
            CALL MPI_BCAST_COMPLEX(CR1X(0,k,j),1,1)
            CALL MPI_BCAST_COMPLEX(CR2X(0,k,j),1,1)
            CALL MPI_BCAST_COMPLEX(CR3X(0,k,j),1,1)       
          enddo
      enddo

! Convert to physical space
      
      

      CALL REAL_FOURIER_TRANS_U1 (.false.)
      CALL REAL_FOURIER_TRANS_U2 (.false.)
      CALL REAL_FOURIER_TRANS_U3 (.false.)
      CALL REAL_FOURIER_TRANS_P (.false.)


! Get the turbulent kinetic energy at each level
      do k=0,NZ+1 
        do j=0,NY+1
          urms(k,j)=0.d0
          vrms(k,j)=0.d0
          wrms(k,j)=0.d0
      do i=0,min(NXP,NXP_L) 
        urms(k,j)=urms(k,j)+(U1X(i,k,j)-dble(CR1X(0,k,j)))**2
        vrms(k,j)=vrms(k,j)+(U2X(i,k,j)-dble(CR2X(0,k,j)))**2
        wrms(k,j)=wrms(k,j)+(U3X(i,k,j)-dble(CR3X(0,k,j)))**2
      end do
!        urms(k,j)=dsqrt(urms(k,j)/(dble(NX)))
!        vrms(k,j)=dsqrt(vrms(k,j)/(dble(NX)))
!        wrms(k,j)=dsqrt(wrms(k,j)/(dble(NX)))
      end do 
      end do



        uv(:,:)=0. 
        uw(:,:)=0.
        wv(:,:)=0.
        pv(:,:)=0.
        pu(:,:)=0.

! Compute the Reynolds stress and turbulent flux
      do k=0,NZ
      do j=0,NY
        uv(k,j)=0. 
        uw(k,j)=0.
        wv(k,j)=0.
        pv(k,j)=0.
        pu(k,j)=0.
      do i=0,min(NXP,NXP_L)
        uv(k,j)=uv(k,j)+(U1X(i,k,j)-dble(CR1X(0,k,j)))             &
                       *0.5*((U2X(i,k,j+1)-dble(CR2X(0,k,j+1))) +  &
                             (U2X(i,k,j)-dble(CR2X(0,k,j))) )

        wv(k,j)=wv(k,j)+ 0.5*((U3X(i,k+1,j)-dble(CR3X(0,k+1,j))) +  &
                              (U3X(i,k,j)-dble(CR3X(0,k,j))) )*     &
                         0.5*((U2X(i,k,j+1)-dble(CR2X(0,k,j+1))) +  &
                              (U2X(i,k,j)-dble(CR2X(0,k,j))) )
 
        uw(k,j)=uw(k,j)+(U1X(i,k,j)-dble(CR1X(0,k,j)))             &
                       *0.5*((U3X(i,k+1,j)-dble(CR3X(0,k+1,j))) +  &
                             (U3X(i,k,j)-dble(CR3X(0,k,j))) ) 

        pu(k,j)=pv(k,j)+0.5*((U3X(i,k+1,j)-dble(CR3X(0,k+1,j))) +  &
                             (U3X(i,k,j)-dble(CR3X(0,k,j))) )      &
                           *(PX(i,k,j)-dble(p_mean(k,j)))

        pv(k,j)=pv(k,j)+ 0.5*((U2X(i,k,j+1)-dble(CR2X(0,k,j+1))) +     &
                              (U2X(i,k,j)-dble(CR2X(0,k,j))) )         &                   
                           *(PX(i,k,j)-dble(p_mean(k,j)))
      end do
        uv(k,j)=uv(k,j)/(float(NX))
        uw(k,j)=uw(k,j)/(float(NX))
        wv(k,j)=wv(k,j)/(float(NX))
        pu(k,j)=pv(k,j)/(float(NX))
        pv(k,j)=pv(k,j)/(float(NX))
      end do
      end do

!Compute the  mean velocity gradient
      do k=1,NZ 
      do j=1,NY
        dudy(k,j)=dble(CR2X(0,k,j)-CR2X(0,k,j-1))/(GYF(j)-GYF(j-1))
        dwdy(k,j)=dble(CR3X(0,k,j)-CR3X(0,k,j-1))/(GYF(j)-GYF(j-1))
      end do
      end do

      do k=1,NZ
       dudy(k,1)=dudy(k,2)
      enddo
      

      do k=1,NZ              
      do j=1,NY
        dudz(k,j)=dble(CR2X(0,k,j)-CR2X(0,k-1,j))/(GZF(k)-GZF(k-1))
        dwdz(k,j)=dble(CR3X(0,k,j)-CR3X(0,k-1,j))/(GZF(k)-GZF(k-1))
      end do
      enddo
      
      do j=1,NY
       dwdz(1,j)=dwdz(2,j)
       dudz(1,j)=dudz(2,j)
      end do
        
 
 

! Calculate the mean square shear
      do k=1,NZ
      do j=1,NY
        shear(k,j)=0.d0
          do i=0,min(NXP,NXP_L)
           shear(k,j)=shear(k,j)                                   &
                 +((U1X(i,k,j+1)-U1X(i,k,j-1))/(2.d0*DYF(j)))**2.d0  &
                 +((U3X(i,k,j+1)-U3X(i,k,j-1))/(2.d0*DYF(j)))**2.d0  &
                 +((U1X(i,k+1,j)-U1X(i,k-1,j))/(2.d0*DZF(k)))**2.d0  &
                 +((U3X(i,k+1,j)-U3X(i,k-1,j))/(2.d0*DZF(k)))**2.d0
          end do
        end do
        shear(k,j)=shear(k,j)/dble(NX)
      end do





 
!c      call tkebudget_chan_1   

!C      Call netcdf
!c       call   NETCDF_WRITE
      

! Do over the number of passive scalars
      do n=1,N_TH

! Save CTHX
      do k=0,NZ+1
        do i=0,NX2P
          do j=0,NY+1
            CRTHX(i,k,j,n)=CTHX(i,k,j,n)
            CFTHX(i,k,j,n)=CTHX(i,k,j,n)
          end do
        end do
      end do

       do k=0,NZV-1
          do j=0,NY+1
           CALL MPI_BCAST_COMPLEX(CRTHX(0,k,j,N),1,1)
          enddo
       enddo

      end do  ! done with passive scalar do loop

! Convert to physical space
 
      if (N_TH .gt. 0) then
       CALL REAL_FOURIER_TRANS_TH (.false.)
      endif


! Do over the number of passive scalars
      do n=1,N_TH

      do j=0,NY+1
      do k=0,NZ+1
        thrms(k,j,n)=0.
      do i=0,min(NXP,NXP_L)
        thrms(k,j,n)=thrms(k,j,n) + (abs(THX(i,k,j,n)  &
                -dble(CRTHX(0,k,j,n))))**2.
      end do
      end do
      end do

! Compute the Reynolds stress and mean velocity gradient
      do j=1,NY
      do k=1,NZ
        thv(k,j,n)=0.
        thw(k,j,n)=0.
      do i=0,min(NXP,NXP_L)
        thv(k,j,n)=thv(k,j,n)+(THX(i,k,j,n)-dble(CRTHX(0,k,j,n)))  &
           *(0.5*(U2X(i,k,j)+U2X(i,k,j+1) )                        &
           - 0.5*dble(CR2X(0,k,j)+CR2X(0,k,j+1)) )

       thw(k,j,n)=thw(k,j,n)+(THX(i,k,j,n)-dble(CRTHX(0,k,j,n)))   &
         *(0.5*(U3X(i,k+1,j)+U3X(i,k,j) )                          &
           - 0.5*dble(CR3X(0,k+1,j)+CR3X(0,k,j)) )
      end do
      thv(k,j,n)=thv(k,j,n)/(float(NX))
      thw(k,j,n)=thw(k,j,n)/(float(NX))
      end do
      end do

! Get the y-derivative of the mean scalar at GYF points
      do k=1,NZ
      do j=1,NY
        dthdy(k,j,n)=dble(CRTHX(0,k,j+1,n)-CRTHX(0,k,j-1,n))/(2.*DYF(j))
      end do
      end do
      
      do k=1,NZ
      do j=1,NY
        dthdz(k,j,n)=dble(CRTHX(0,k+1,j,n)-CRTHX(0,k-1,j,n))/(2.*DZF(K))
      end do
      end do     

      UBULK= 0.0d0
      area = 0.0d0

      DO J=2,NY
       DO K=2,NZ
        area  = area + DY(j)*DZ(k)
        UBULK=UBULK+0.25*(dble(CRTHX(0,K,J,n))+dble(CRTHX(0,K-1,J,n)) +  &
          dble(CRTHX(0,K-1,J-1,n)) + dble(CRTHX(0,K,J-1,n)) )*DY(j)*DZ(K)
       ENDDO
      ENDDO
      UBULK=UBULK/area

! Write out UBULK
      If (rank .eq. 0) then
! Write out UBULK
      write(*,*) 'THBULK: ',UBULK, 'Ri', RI_TAU(N)
      write(99,*) 'THBULK: ',UBULK, 'Ri', RI_TAU(N)
      endif
      
! End do over number of passive scalars, n
      end do

!     IF (XY_PLANE_comb)THEN
!       CALL  MPI_COMBINE_XY (U1X,U1_XY,NXP+1,NY+2)
!       CALL  MPI_COMBINE_XY (U2X,U2_XY,NXP+1,NY+2) 
!       CALL  MPI_COMBINE_XY (U3X,U3_XY,NXP+1,NY+2)
!       CALL  MPI_COMBINE_XY (THX,TH_XY,NXP+1,NY+2) 
!      ENDDO

       if (rank .eq. 0) then
       open(66,file='plane_data_yz/time_bulk.txt',form='formatted', &
       status='unknown', position='append' )
       write(6,*) 'TIME STEP = ', TIME_STEP, ' TIME= ', TIME/60.0,' dt = ', DELTA_T
       write(66,565) TIME/(60.0d0*60.0d0),DELTA_T,UBULK
       close(66)
565    format(f12.5,2f13.8)
       
      k = time_step/SAVE_STATS_INT
      call plane_parav_plane_yz(k)

      end if
  

      IF (XY_PLANE_comb)THEN
        k = time_step/SAVE_STATS_INT
        call plane_parav_plane_xy(k)
      ENDIF

      IF (XY_PLANE_comb_2)THEN
        k = time_step/SAVE_STATS_INT
        call plane_parav_plane_xy_2(k)
      ENDIF

      IF (ZX_PLANE_comb) THEN
        k=time_step/SAVE_STATS_INT
        call plane_parav_plane_zx(k)
      ENDIF

      IF (ZX_PLANE_comb_2) THEN
        k=time_step/SAVE_STATS_INT
        call plane_parav_plane_zx_2(k)
      ENDIF

      IF (PLANE_3D) THEN
        k=time_step/SAVE_STATS_INT
        call plane_parav_3D(k)
      ENDIF




!        IF (XY_PLANE) THEN
!!!!!!!!!!!!SAVING SPAN-WISE PLANE DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      i = int((NZ+1)/2) !(NZ+1)/2
!      k = time_step/SAVE_STATS_INT
!      call plane_XY_binary(k,i)
!      ENDIF

!      IF (ZX_PLANE) THEN
!!!!!!!!!!!!SAVING SPAN-WISE PLANE DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      j = NY !int((NY+1)/2) !(NZ+1)/2
!      k = time_step/SAVE_STATS_INT
!      call plane_ZX_binary(k,j)
!      ENDIF

!      IF (ZX_PLANE_2) THEN
!      j = NY-16 !int((NY+1)/2) !(NZ+1)/2
!      k = time_step/SAVE_STATS_INT
!      call plane_ZX_binary_2(k,j)
!      ENDIF


      IF (SAVE_3D) THEN
       k = time_step/SAVE_STATS_INT
       call plane_3D_binary(k)
      ENDIF

      if ( MEAN_BACKGROUND ) then 
        call  mean_background_potential_duct
     endif

! Convert back to Fourier space
     if (N_TH .gt. 0) then 
       CALL REAL_FOURIER_TRANS_TH (.true.)
     endif


!      write(6,*) 'I am here 1', rank

!      CALL MPI_COMBINE_STATS(urms,NZ+2,NY+2)

!      write(6,*) 'I am here 2', rank
!      CALL MPI_COMBINE_STATS(vrms,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(wrms,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(uv,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(uw,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(wv,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(pu,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(pv,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(shear,NZ+2,NY+2) 

!      CALL MPI_COMBINE_STATS(dudz,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(dudy,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(dwdz,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(dwdy,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(omega_x,NZ+2,NY+2)


! Get the bulk rms value

      if(rank .eq. 0) then

      do k=0,NZ+1
      do j=0,NY+1
        urms(k,j)=dsqrt(urms(k,j)/(dble(NX)))
        vrms(k,j)=dsqrt(vrms(k,j)/(dble(NX)))
        wrms(k,j)=dsqrt(wrms(k,j)/(dble(NX)))
      end do
      end do      
  
      urms_b=0.d0
      vrms_b=0.d0
      wrms_b=0.d0
      area  = 0.d0
      do k=1,NZ
      do j=1,NY
        area  = area + DY(j)*DZ(k)
        urms_b=urms_b + 0.25d0*(urms(k,j)+urms(k,j-1) &
            +  urms(k-1,j)+urms(k-1,j-1))*DY(j)*DZ(k)
        vrms_b=vrms_b + 0.25d0*(vrms(k,j)+vrms(k,j-1) &
            +  vrms(k-1,j)+vrms(k-1,j-1))*DY(j)*DZ(k)
        wrms_b=wrms_b + 0.25d0*(wrms(k,j)+wrms(k,j-1) &
            +  wrms(k-1,j)+wrms(k-1,j-1))*DY(j)*DZ(k)
      end do
      enddo
      urms_b=urms_b/(area)
      vrms_b=vrms_b/(area)
      wrms_b=wrms_b/(area)
 

! Write out the bulk rms velocity
      write(6,*) '<U_rms>: ',urms_b
      write(6,*) '<V_rms>: ',vrms_b, area
      write(6,*) '<W_rms>: ',wrms_b


      endif

      Do n =1,N_th
!      CALL MPI_COMBINE_STATS(thrms(0,0,n),NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(thv(0,0,n),NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(thw(0,0,n),NZ+2,NY+2)

!      CALL MPI_COMBINE_STATS(dthdz(0,0,n),NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(dthdy(0,0,n),NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(Rig,NZ+2,NY+2)

      if (rank .eq. 0) then
      do j=0,NY+1
      do k=0,NZ+1
        thrms(k,j,n)=sqrt(thrms(k,j,n)/float(NX))
      end do
      end do
      endif

      enddo

      IF (MEAN_KE_BUDGET) THEN
        call mean_energy_budget_duct
      ENDIF 
      !!!!!!!!!!!!!TKE BUDGET!!!!!!!!!!!!!!!!!!!!!
      IF (TKE_BUDGET) THEN
       call tke_budget_duct
      ENDIF

      
!C Convert velocity back to Fourier space
! C Convert velocity back to Fourier space
      CALL REAL_FOURIER_TRANS_U1 (.true.) 
      CALL REAL_FOURIER_TRANS_U2 (.true.)
      CALL REAL_FOURIER_TRANS_U3 (.true.)
      CALL REAL_FOURIER_TRANS_P (.true.)

      CF1X = (0.d0,0.d0)
      CF2X = (0.d0,0.d0)
      CF3X = (0.d0,0.d0) 

       
  

!     Saving files for paraview   
      if (rank .eq. 0 ) then
       write(6,*) 'done save_stats chan'
       k = time_step/SAVE_STATS_INT
       call plane_parav_vel(k)
      endif

! Deallocating stat variables
!     call deallocate_temps

!     CCCCCCCCCCCCCCCCCCCCCCCCCC
      if (rank .eq. 0 ) then
      write(*,*) 'Deallocate all the tmp arrays'
      endif

      RETURN
      END

subroutine tke_budget_duct
! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps
! This subroutine should be called in SAVE_STATS_CHAN after computing
! plane averaged statistics, with the velocity in physical space, and 
! CRi containing the velocity in Fourier space

use ntypes
use Domain
use Grid
use Fft_var, only : CIKXP
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use variable_stat
use mpi_var      
use les_chan_var

implicit none

      integer i,j,k,n

      

!      CHARACTER*28 file_tke

      if (rank .eq. 0) then    
      write(6,*) 'Calculating  tke budget', rank
      endif

!      tke_mean(:,:) = 0.0d0
      tke_2_1(:,:)  = 0.0d0
      tke_2_2(:,:)  = 0.0d0
      tke_3(:,:)    = 0.0d0
      tke_4(:,:)    = 0.0d0
   

      IF (TIME_STEP .eq. 0) THEN        
       tke_1(:,:) = 0.0d0
      ELSE 
      do j=0,NY
       do k=0,NZ
! tke_mean defined at GY points
!        tke_mean_old(k,j)=tke_mean(k,j)
        tke_mean(k,j)=0.5d0*( urms(k,j)**2.d0         &
             + 0.25*(vrms(k,j) + vrms(k,j+1))**2.d0     &
             + 0.25*(wrms(k,j) + wrms(k+1,j))**2.d0 ) 
        tke_1(k,j)=(tke_mean(k,j)-tke_mean_old(k,j))/(TIME-TIME_old)
       end do
      end do
      ENDIF

!      time_old=TIME

!     U2*dtke/dy

      do j=1,NY-1
       do k=1,NZ-1
        tke_2_1(k,j)=-dble(CR2X(0,K,J))*(tke_mean(k,j+1)-tke_mean(k,j-1))/(2.0d0*DYF(j))                  
       end do
      end do

!     U3*dtke/dz

      do j=1,NY-1
      do k=1,NZ-1
       tke_2_2(k,j)=-dble(CR3X(0,K,J))*(tke_mean(k+1,j)-tke_mean(k-1,j))/(2.0d0*DZF(k))  
       tke_2(k,j) = tke_2_1(k,j) + tke_2_2(k,j)
      end do
      end do
       


! Get the production at GYF and GZF points (cell center)
      do j=1,NY
      do k=1,NZ
!        tke_3(k,j)=       -uv(k,j)*dUdy(k,j)-uw(k,j)*dUdz(k,j)  
!                          -vv(k,j)*dVdy(k,j)-vw(k,j)*dVdz(k,j)
!                          -wv(k,j)*dWdy(k,j)-ww(k,j)*dWdy(k,j)
       tke_3(k,j)=-uv(k,j)*dble(CR1X(0,k,j+1)-CR1X(0,k,j-1))/(2.0d0*DYF(j))  &
                  -uw(k,j)*dble(CR1X(0,k+1,j)-CR1X(0,k-1,j))/(2.0d0*DZF(k))  &
                  -0.25d0*(vrms(k,j)+vrms(k,j+1))**2.0*dble(CR2X(0,k,j+1)-CR2X(0,k,j))/DYF(j)     &
                  -wv(k,j)*0.50d0*dble(CR2X(0,k+1,j+1)+CR2X(0,k+1,j)-CR2X(0,k-1,j+1)-CR2X(0,k-1,j))/(2.0d0*DZF(k))   &
                  -wv(k,j)*0.50d0*dble(CR3X(0,k+1,j+1)+CR3X(0,k,j+1)-CR3X(0,k+1,j-1)-CR3X(0,k,j-1))/(2.0d0*DYF(j))   &
                  -0.25d0*(wrms(k,j)+wrms(k+1,j))**2.0*dble(CR3X(0,k+1,j)-CR3X(0,k,j))/DZF(k)     
      end do
      end do



!    viscous diffusion NU*d2tke/dxidxi

      do j=1,NY-1
       do k=1,NZ-1
        tke_4(k,j)= NU*( ( (tke_mean(K+1,J) - tke_mean(K,J)) / DZ(K+1)              &
                          -(tke_mean(K,J)   - tke_mean(K-1,J)) / DZ(K)) /DZF(K)               &
                 +        ( (tke_mean(K,J+1) - tke_mean(K,J)) / DY(J+1)             &
                           -(tke_mean(K,J)   - tke_mean(K,J-1)) / DY(J)) /DYF(J)  )           
      end do 
      end do

      do n=1,N_TH
        do j=1,NY
         do k=1,NZ
!          tke_5(k,j,n)=-RI_TAU(n)*thv(k,j,n)
           tke_5(k,j,n)= RI_TAU(N)*thv(k,j,n)
         end do
        end do
      end do
  
! Construct the resolved turbulent transport terms 
! Convert the pressure to physical space for computation of the
! pressure transport
! Get the mean of the pressure
!      do j=0,NY+1
!       do k=0,NZ+1
!        p_mean(k,j)=dble(CF1X(0,k,j))
!       end do
!      end do


      transport(:,:) = 0.0d0;

      do j=1,NY+1      
      do k=0,NZ+1
        transport(k,j)=0.d0
        do i=0,min(NXP,NXP_L)
          transport(k,j)=transport(k,j)     &
! Vertical Pressure Transport term:
             +0.50d0*(U2X(I,K,J) -dble(CR2X(0,k,j)) )*(PX(I,K,J)+PX(I,K,J-1)-p_mean(k,J)-p_mean(k,J-1)) 
          end do        
        transport(k,j)=transport(k,j)/dble(NX)
      end do
      end do

     CALL MPI_COMBINE_STATS(transport,NZ+2,NY+2)

      do j=1,NY
       do k=1,NZ
        tke_6_1_1(k,j)=- (transport(k,j+1)-transport(k,j))/DYF(j)
       end do 
      end do

      do j=0,NY+1      
      do k=1,NZ+1
        transport(k,j)=0.d0
        do i=0,min(NXP,NXP_L)
          transport(k,j)=transport(k,j)  &
! Horizontal Pressure Transport term:
             + 0.50d0*(U3X(I,K,J)-dble(CR3X(0,k,j)))*(PX(I,K,J)+PX(I,K-1,J)-p_mean(k,J)-p_mean(k-1,J))
          end do        
        transport(k,j)=transport(k,j)/dble(NX)
      end do
      end do

      CALL MPI_COMBINE_STATS(transport,NZ+2,NY+2)

      do j=1,NY
       do k=1,NZ
        tke_6_1_2(k,j)=-(transport(k+1,j)-transport(k,j))/DZF(k)
        tke_6_1(k,j)=tke_6_1_1(k,j) + tke_6_1_2(k,j)
       end do 
      end do

! Turbulent transport terms:
!   d(0.5*u_i^2.0*v)dy


      do j=1,NY+1
      do k=1,NZ
      transport(k,j)=0.d0
          do i=0,min(NXP,NXP_L)
            transport(k,j)=transport(k,j)             &
! u1^2*u2
          + 0.5d0*0.250d0*(U1X(I,K,J)+U1X(I,K,J-1)-dble(CR1X(0,k,J)+CR1X(0,k,J-1)))**2.d0   &
                 *(U2X(I,K,J)- dble(CR2X(0,k,J)) )       &
! U2^3
          + 0.5d0*(U2X(I,K,J)- dble(CR2X(0,k,J)))**3.d0  &
! U3^2*U2
          + 0.5d0*0.06250d0*(U3X(I,K+1,J)+U3X(I,K+1,J-1)+U3X(I,K,J)+U3X(I,K,J-1) &
                 -dble(CR3X(0,k+1,J)+CR3X(0,k+1,J-1)+CR3X(0,k,J)+CR3X(0,k,J-1)))**2.d0    &
                 *(U2X(I,K,J)- dble(CR2X(0,k,J)))
         end do
        transport(k,j)=transport(k,j)/dble(NX)
      end do
      end do

      CALL MPI_COMBINE_STATS(transport,NZ+2,NY+2)
! Now, the vertical derivative of the transport term:
      do j=1,NY
      do k=1,NZ
        tke_6_2_1(k,j)=-(transport(k,j+1)-transport(k,j))/DYF(j)
      end do
      end do


!   d(0.5*u_i^2.0*w)dz

      do j=1,NY
      do k=1,NZ+1
      transport(k,j)=0.d0
          do i=0,min(NXP,NXP_L)
            transport(k,j)=transport(k,j)   &
! U1^2*U3
          + 0.5d0*0.250d0*(U1X(I,K,J)+U1X(I,K-1,J)-dble(CR1X(0,K,J)+CR1X(0,K-1,J)))**2.d0 &
                 *(U3X(I,K,J)- dble(CR3X(0,K,J)) ) &
! U3^3
          + 0.5d0*(U3X(I,K,J)- dble(CR3X(0,K,J)))**3.d0 &
! U2^2.0*U3
          + 0.5d0*0.06250d0*(U2X(I,K,J+1)+U2X(I,K-1,J+1)+U2X(I,K,J)+U2X(I,K-1,J) &
                 -dble(CR2X(0,K,J+1)+CR2X(0,K-1,J+1)+CR2X(0,K,J)+CR2X(0,K-1,J)))**2.d0 &
                 *(U3X(I,K,J)- dble(CR3X(0,K,J)))
         end do
        transport(k,j)=transport(k,j)/dble(NX)
        end do
      end do

      CALL MPI_COMBINE_STATS(transport,NZ+2,NY+2)
! Now, the horizontal derivative of the transport term:

      do j=1,NY
       do k=1,NZ
        tke_6_2_2(k,j)=-(transport(k+1,j)-transport(k,j))/DZF(k)   
                 
        tke_6_2(k,j) = tke_6_2_1(k,j) + tke_6_2_2(k,j)
       end do
      end do
       
 
! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
      do j=0,NY+1
       do k=0,NZ+1
        epsilon(k,j)=0.d0
       end do
      end do

! Store du/dx in CS1
      do j=1,NY
      do k=1,NZ
      do i=0,NX2P
        CS1X(i,k,j)=CIKXP(i)*CF1X(i,k,j) ! WE HAVE ALSO STORED CU1 IN CF1X
!        CS1(i,k,j)=CIKXP(i)*CR1X(i,k,j) 
      end do
      end do
      end do
! Convert to physical space

      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
      ENDDO
      S1Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)

!C      call fft_x_to_physical(CS1,S1X,0,NY+1,0,NZ+1)


      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
!        epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
         epsilon(k,j)=epsilon(k,j) + S1X(i,k,j)**2.0
      end do
      end do
      end do
! Store dv/dx in CS1
      do j=1,NY
      do k=1,NZ
      do i=0,NX2P
        CS1X(i,k,j)=CIKXP(i)*CF2X(i,k,j)   ! WE HAVE ALSO STORED CU2 IN CF2X
!        CS1(i,k,j)=CIKXP(i)*CR2X(i,k,j)
      end do
      end do
      end do

! Convert to physical space
      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
      ENDDO
      S1Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)

!C       call fft_x_to_physical(CS1,S1X,0,NY+1,0,NZ+1)
      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
      end do
      end do
      end do
      
 ! Store dw/dx in CS1
      do j=1,NY
      do k=1,NZ
      do i=0,NX2P
        CS1X(i,k,j)=CIKXP(i)*CF3X(i,k,j)  ! WE HAVE ALSO STORED CU3 IN CF3X
      end do
      end do
      end do

! Convert to physical space
      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
      ENDDO
      S1Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)

!C      call fft_x_to_physical(CS1,S1X,0,NY+1,0,NZ+1)
      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
         epsilon(k,j)=epsilon(k,j)+ (S1X(i,k,j)**2.0)
      end do
      end do
      end do
     
      
! Compute du/dy note remove mean
      do j=0,NY+1
      do k=0,NZ+1
      do i=0,min(NXP,NXP_L)
         S2X(i,k,j)= U1X(i,k,j)-dble(CR1X(0,k,j))
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
       S1X(i,k,j) = (S2X(i,k,j+1)-S2X(i,k,j-1))/(2.0d0*DYF(j))
      end do
      end do
      end do


      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
      end do
      end do
      end do
      
      
      


      
! Compute du/dz 
      do j=0,NY+1
      do k=0,NZ+1
      do i=0,min(NXP,NXP_L)
         S2X(i,k,j)= U1X(i,k,j)-dble(CR1X(0,k,j))
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        S1X(i,k,j) = (S2X(i,k+1,j)-S2X(i,k-1,j))/(2.0d0*DZF(k))
      end do
      end do
      end do    

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
         epsilon(k,j)=epsilon(k,j)+ (S1X(i,k,j)**2.0)
      end do
      end do
      end do

! Compute dv/dy  note remove mean

      do j=0,NY+1
      do k=0,NZ+1
      do i=0,min(NXP,NXP_L)
       S2X(i,k,j)= U2X(i,k,j)-dble(CR2X(0,k,j))
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
       S1X(i,k,j) = (S2X(i,k,j+1)-S2X(i,k,j))/DYF(j)
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
      end do
      end do
      end do

! Compute dw/dy, note remove mean
      do j=0,NY+1
      do k=0,NZ+1
      do i=0,min(NXP,NXP_L)
       S2X(i,k,j)= U3X(i,k,j)-dble(CR3X(0,k,j))
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
       S1X(i,k,j)=0.50d0*(S2X(i,k+1,j+1)+S2X(i,k,j+1)-S2X(i,k+1,j-1)-S2X(i,k,j-1))/(2.0d0*DYF(j))
      end do
      end do
      end do
     

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
      end do
      end do
      end do

! Compute dv/dz, note remove mean
      do j=0,NY+1
      do k=0,NZ+1
      do i=0,min(NXP,NXP_L)
         S2X(i,k,j)= U2X(i,k,j)-dble(CR2X(0,k,j))
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        S1X(i,k,j) = 0.50d0*(S2X(i,k+1,j+1)+S2X(i,k+1,j)-S2X(i,k-1,j+1)-S2X(i,k-1,j))/(2.0d0*DZF(k))
      end do
      end do
      end do 

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
      end do
      end do
      end do

! Store dw/dz in CS1
      do j=0,NY+1
      do k=0,NZ+1
      do i=0,min(NXP,NXP_L)
         S2X(i,k,j)= U3X(i,k,j)-dble(CR3X(0,k,j))
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        S1X(i,k,j) =(S2X(i,k+1,j)-S2X(i,k,j))/DZF(k)                      
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
         epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
      end do
      end do
      end do

!!!!!! NEED ADD FOR ALL PROCESSORS!!!!!!!!!!!!
      CALL MPI_COMBINE_STATS(epsilon,NZ+2,NY+2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=1,NY
      do k=1,NZ
        tke_7(k,j)=-NU*epsilon(k,j)/float(NX)
      end do
      end do

      
      do k=1,NZ
        tke_7(k,0)=tke_7(k,1)
      end do


       if (rank .eq. 0) then
       open(66,file='plane_tke/tke_budget.dat',form='formatted', &
       status='unknown' , position='append' )

       open(77,file='plane_tke/tke_budget_sgs.dat',form='formatted', &
       status='unknown' , position='append' )       

      sum_Prod  = 0.0d0 ;
      sum_dissp = 0.0d0 ;
      sum_buoyF = 0.0d0 ;       
      sum_dtkdt = 0.0d0 ;
      sum_vis   = 0.0d0 ;
      sum_trans = 0.0d0 ;
      sum_advec = 0.0d0 ;
      sum_dissp_sgs = 0.0d0 
      sum_prod_sgs   = 0.0d0

      area = 0.0d0 ;

 
      do j=1,NY
      do k=1,NZ
       area  = area + DY(j)*DZ(k)
       sum_Prod  = sum_Prod  + tke_3(k,j)*DY(j)*DZ(k)
       sum_dissp = sum_dissp + tke_7(k,j)*DY(j)*DZ(k)
       sum_buoyF = sum_buoyF + tke_5(k,j,1)*DY(j)*DZ(k)
       sum_dtkdt = sum_dtkdt + tke_1(k,j)*DY(j)*DZ(k)
       sum_vis   = sum_vis   + tke_4(k,j)*DY(j)*DZ(k) 
       sum_trans = sum_trans + (tke_6_1(k,j) + tke_6_2(k,j))*DY(j)*DZ(k) 
       sum_advec = sum_advec + tke_2(k,j)*DY(j)*DZ(k) 
       IF(LES)THEN
       sum_prod_sgs = sum_prod_sgs + tke_sgs_p(k,j)*DY(j)*DZ(k)
       sum_dissp_sgs = sum_dissp_sgs + tke_sgs_diss(k,j)*DY(j)*DZ(k)
       ENDIF
      end do
      end do

      sum_Prod = sum_Prod/area ;
      sum_dissp = sum_dissp/area ;
      sum_buoyF = sum_buoyF/area ;
      sum_dtkdt = sum_dtkdt/area ;
      sum_vis   = sum_vis/area ;
      sum_trans = sum_trans/area ;
      sum_advec = sum_advec/area ;
 
      sum_Prod_sgs = sum_Prod_sgs/area ;
      sum_dissp_sgs = sum_dissp_sgs/area ;

      write(66,566) TIME/(60.0d0*60.0d0),sum_dtkdt,sum_advec, sum_Prod, sum_dissp, sum_buoyF, sum_vis, &
                    sum_Prod + sum_dissp + sum_buoyF

      write(77,566) TIME/(60.0d0*60.0d0),sum_Prod, sum_dissp,sum_Prod_sgs,-abs(sum_dissp_sgs)
      close(66)
      close(77)

566   format(f10.5,7E18.9) 
      endif


      IF (RANK .eq. 0 ) THEN
      k = time_step/SAVE_STATS_INT
      call plot_para_tke(k)        
      ENDIF

      if (rank .eq. 0) then
      write(6,*) 'Done with  tke budget', rank
      endif

      return 
       end
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!-------------------------------------C---
subroutine mean_energy_budget_duct
!-------------------------------------C--



! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps
! This subroutine should be called in SAVE_STATS_CHAN after computing
! plane averaged statistics, with the velocity in physical space, and 
! CRi containing the velocity in Fourier space
use ntypes
use Domain
use Grid
use Fft_var, only :  CIKXP
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use variable_stat
use mpi_var
implicit none

      integer i,j,k,n
      CHARACTER*31 file_energy

     if (rank .eq. 0) then
       write(6,*) 'Calculating energy_budget', rank
      endif

      tke_2_1(:,:)  = 0.0d0
      tke_2_2(:,:)  = 0.0d0
      tke_4(:,:) = 0.0d0

      IF (TIME_STEP .eq. 0) THEN        
       tke_1(:,:) = 0.0d0
      ELSE
      do j=0,NY
       do k=0,NZ
! tke_mean defined at GY points
!        energy_mean_old(k,j)=energy_mean(k,j)
        energy_mean(k,j)=0.5d0*( dble(CR1X(0,k,j))**2.d0  &
             + 0.25d0*dble(CR2X(0,k,j+1)+CR2X(0,k,j))**2.d0                 &
             + 0.25d0*dble(CR3X(0,k+1,j)+CR3X(0,k,j))**2.d0 )
        tke_1(k,j)=(energy_mean(k,j)-energy_mean_old(k,j)) &
             /(TIME-TIME_old_eng)
       end do
      end do
      ENDIF 

!      time_old_eng = TIME

!     U2*d<E>/dy

      do j=1,NY-1
       do k=1,NZ-1
        tke_2_1(k,j)=-dble(CR2X(0,K,J))*(energy_mean(k,j+1)-energy_mean(k,j-1))/(2.0d0*DYF(j))
       end do
      end do

!     U3*d<E>/dz

      do j=1,NY-1
      do k=1,NZ-1
       tke_2_2(k,j)=-dble(CR3X(0,K,J))*(energy_mean(k+1,j)-energy_mean(k-1,j))/(2.0d0*DZF(k))  
       tke_2(k,j) = tke_2_1(k,j) + tke_2_2(k,j)
      end do
      end do


    


! Get the production 

      do j=1,NY
      do k=1,NZ
!        tke_3(k,j)=       -uv(k,j)*dUdy(k,j)-uw(k,j)*dUdz(k,j)  
!                          -vv(k,j)*dVdy(k,j)-vw(k,j)*dVdz(k,j)
!                          -wv(k,j)*dWdy(k,j)-ww(k,j)*dWdy(k,j)
       tke_3(k,j)=uv(k,j)*dble(CR1X(0,k,j+1)-CR1X(0,k,j-1))/(2.0d0*DYF(j))  &
                  +uw(k,j)*dble(CR1X(0,k+1,j)-CR1X(0,k-1,j))/(2.0d0*DZF(k))  &
                  +0.25d0*(vrms(k,j)+vrms(k,j+1))**2.0*dble(CR2X(0,k,j+1)-CR2X(0,k,j))/DYF(j)     &
                  +wv(k,j)*0.50d0*dble(CR2X(0,k+1,j+1)+CR2X(0,k+1,j)-CR2X(0,k-1,j+1)-CR2X(0,k-1,j))/(2.0d0*DZF(k))   &
                  +wv(k,j)*0.50d0*dble(CR3X(0,k+1,j+1)+CR3X(0,k,j+1)-CR3X(0,k+1,j-1)-CR3X(0,k,j-1))/(2.0d0*DYF(j))   &
                  +0.25d0*(wrms(k,j)+wrms(k+1,j))**2.0*dble(CR3X(0,k+1,j)-CR3X(0,k,j))/DZF(k)   
      end do
      end do


!    viscous diffusion NU*d2<ke>/dxidxi

      do j=1,NY-1
       do k=1,NZ-1
        tke_4(k,j)= NU*( ( (energy_mean(K+1,J) - energy_mean(K,J)) / DZ(K+1)              &
                          -(energy_mean(K,J)   - energy_mean(K-1,J)) / DZ(K)) /DZF(K)               &
                 +        ( (energy_mean(K,J+1) - energy_mean(K,J)) / DY(J+1)             &
                           -(energy_mean(K,J)   - energy_mean(K,J-1)) / DY(J)) /DYF(J)  ) 
      end do
      end do

      do n=1,N_TH
        do j=1,NY
         do k=1,NZ
!          tke_5(k,j,n)=-RI_TAU(n)*dble(CTHX(0,k,j,n))*dble(CR2X(0,k,j))
           tke_5(k,j,n)=RI_TAU(N)*dble(CTHX(0,k,j,n))*dble(CR2X(0,k,j))
         end do
        end do
      end do

!     write(6,*) 'I am here ', rank

! Construct the resolved turbulent transport terms 
! Convert the pressure to physical space for computation of the
! pressure transport
! Get the mean of the pressure
!      do j=0,NY+1
!       do k=0,NZ+1
!        p_mean(k,j)=dble(CF1X(0,k,j))
!       end do
!      end do

      transport(:,:) = 0.d0
      
      do j=1,NY+1      
      do k=0,NZ+1
! Vertical Pressure Transport term:
        transport(k,j)=0.50d0*dble(CR2X(0,k,j))*(p_mean(k,J)+p_mean(k,J-1))        
      end do
      end do

      do j=1,NY
       do k=1,NZ
        tke_6_1_1(k,j)=- (transport(k,j+1)-transport(k,j))/DYF(j)
       end do 
      end do

      do j=0,NY+1      
      do k=1,NZ+1
! Horizontal (streamwise) Pressure Transport term:
          transport(k,j)=0.50d0*dble(CR3X(0,k,j))*(p_mean(k,J)+p_mean(k-1,J))
      end do
      end do

      do j=1,NY
       do k=1,NZ
        tke_6_1_2(k,j)=-(transport(k+1,j)-transport(k,j))/DZF(k)
        tke_6_1(k,j)=tke_6_1_1(k,j) + tke_6_1_2(k,j)
       end do 
      end do

! Turbulent transport terms:
!   d(0.5*u_i^2.0*v)dy

      transport(:,:) = 0.d0 ;

      do j=0,NY
       do k=0,NZ
        transport(k,j)=                 &
! <U1*u2>*<U1>
          + uv(k,j)*dble(CR1X(0,K,J))                   &
! <u3*u2><U3>
          + 0.50d0*wv(k,j)*dble(CR3X(0,K+1,J)+CR3X(0,K,J))                   &
! <u2*u2>*<U2>
          + 0.50d0*vrms(k,j)**2.0*dble(CR2X(0,K,J+1)+CR2X(0,K,J))

       end do
      end do

! Now, the vertical derivative of the transport term:

      do j=1,NY
      do k=1,NZ
        tke_6_2_1(k,j)=-(transport(k,j+1)-transport(k,j-1))/(2.0d0*DYF(j))
      end do
      end do


!   d(0.5*u_i^2.0*w)dz
      
      transport(:,:) = 0.d0 ;
      do j=0,NY
       do k=0,NZ
         transport(k,j)=  &
! <u1*u3>*<U1>
          + uw(k,j)*dble(CR1X(0,K,J))     &
! <u2*u3>*<U2>
          + 0.50d0*wv(k,j)*dble(CR2X(0,K,J+1)+CR2X(0,K,J))     &
! <u3*u3>*<U3>
          + 0.50d0*wrms(k,j)**2.0*dble(CR3X(0,K+1,J)+CR3X(0,K,J))
        end do
      end do

! Now, the horizontal derivative of the transport term:

      do j=1,NY
       do k=1,NZ
        tke_6_2_2(k,j)=-(transport(k+1,j)-transport(k-1,j))/(2.0d0*DZF(k))
        tke_6_2(k,j) = tke_6_2_1(k,j) + tke_6_2_2(k,j)
       end do
      end do
       
!      write(6,*) 'I am here 1', rank 
 
! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
      do j=0,NY+1
       do k=0,NZ+1
        epsilon(k,j)=0.d0
       end do
      end do

      
      
! Compute du/dy note remove mean

      do j=1,NY
      do k=1,NZ
       tke_7(k,j) = (dble(CR1X(0,k,j+1))-dble(CR1X(0,k,j-1)))/(2.0d0*DYF(j))

       epsilon(k,j)=epsilon(k,j)+(tke_7(k,j)**2.0)
      end do
      end do

      
      
! Compute du/dz 
      do j=1,NY
      do k=1,NZ
        tke_7(k,j) = (dble(CR1X(0,k+1,j))-dble(CR1X(0,k-1,j)))/(2.0d0*DZF(k))

        epsilon(k,j)=epsilon(k,j)+ (tke_7(k,j)**2.0)
      end do
      end do
    

! Compute d<V>/dy  

      do j=1,NY
      do k=1,NZ
       tke_7(k,j) = (dble(CR2X(0,k,j+1))-dble(CR2X(0,k,j)))/DYF(j)          
                
 
       epsilon(k,j)=epsilon(k,j)+ (tke_7(k,j)**2.0)
      end do
      end do

! Compute dw/dy

      do j=1,NY
      do k=1,NZ
       tke_7(k,j)= 0.50d0*(CR3X(0,k+1,j+1)+CR3X(0,k,j+1)-CR3X(0,k+1,j-1)-CR3X(0,k,j-1))/(2.0d0*DYF(j)) 

       epsilon(k,j)=epsilon(k,j)+ (tke_7(k,j)**2.0)
      end do
      end do
     
! Store d<V>/dz 
      

      do j=1,NY
      do k=1,NZ
        tke_7(k,j) = 0.50d0*(CR2X(0,k+1,j+1)+CR2X(0,k+1,j)-CR2X(0,k-1,j+1)-CR2X(0,k-1,j))/(2.0d0*DZF(k))
         epsilon(k,j)=epsilon(k,j)+ (tke_7(k,j)**2.0)
      end do
      end do 


! Store d<W>/dz 
      

      do j=1,NY
      do k=1,NZ
        tke_7(k,j) = (dble(CR3X(0,k+1,j))-dble(CR3X(0,k,j)))/DZF(k) 
        epsilon(k,j)=epsilon(k,j)+ (tke_7(k,j)**2.0)
      end do
      end do

      do j=1,NY
      do k=1,NZ
        tke_7(k,j)=-NU*epsilon(k,j)
      end do
      end do

!       write(6,*) 'I am here 2', rank
      
      do k=1,NZ
        tke_7(k,0)=tke_7(k,1)
      end do

      if (rank .eq. 0) then
       open(66,file='plane_energy/mean_ke_budget.dat',form='formatted', &
       status='unknown' , position='append' )
       
      sum_Prod  = 0.0d0 ;
      sum_dissp = 0.0d0 ;
      sum_buoyF = 0.0d0 ;       
      sum_dtkdt = 0.0d0 ;
      sum_vis   = 0.0d0 ;
      sum_trans = 0.0d0 ;
      sum_advec = 0.0d0 ;

      area = 0.0d0 ;

 
      do j=1,NY
      do k=1,NZ
       area  = area + DY(j)*DZ(k)
       sum_Prod  = sum_Prod  + tke_3(k,j)*DY(j)*DZ(k)
       sum_dissp = sum_dissp + tke_7(k,j)*DY(j)*DZ(k)
       sum_buoyF = sum_buoyF + tke_5(k,j,1)*DY(j)*DZ(k)
       sum_dtkdt = sum_dtkdt + tke_1(k,j)*DY(j)*DZ(k)
       sum_vis   = sum_vis   + tke_4(k,j)*DY(j)*DZ(k) 
       sum_trans = sum_trans + (tke_6_1(k,j) + tke_6_2(k,j))*DY(j)*DZ(k) 
       sum_advec = sum_advec + tke_2(k,j)*DY(j)*DZ(k) 
      end do
      end do

      sum_Prod = sum_Prod/area ;
      sum_dissp = sum_dissp/area ;
      sum_buoyF = sum_buoyF/area ;
      sum_dtkdt = sum_dtkdt/area ;
      sum_vis   = sum_vis/area ;
      sum_trans = sum_trans/area ;
      sum_advec = sum_advec/area ;

      write(66,566) TIME/(60.0d0*60.0d0),sum_dtkdt,sum_advec, sum_Prod, sum_dissp, sum_buoyF, sum_vis, &
                    sum_Prod + sum_dissp + sum_buoyF
      close(66)

566   format(f10.5,7E18.9) 
      endif

      if (rank .eq. 0) then
      write(6,*) 'Done with energy_budget', rank
      endif

      if (rank .eq. 0 ) then
       k = time_step/SAVE_STATS_INT
       call plot_para_eng(k)
      endif




      return
      end


subroutine precalculation_for_budget 

use ntypes
use Domain
use Grid
use Fft_var, only :  pi, CIKX
use TIME_STEP_VAR
use run_variable
use variable_stat
use mpi_var

implicit none

        integer i,j,k,n
!  allocating stat variables      
!      call allocate_temps
!C Apply Boundary conditions to velocity field
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
        CALL APPLY_BC_VEL_LEFT
        CALL APPLY_BC_VEL_RIGHT

        time_old_eng = TIME
        time_old     = TIME
        do j=0,NY
        do k=0,NZ
         energy_mean_old(k,j)=0.5d0*( dble(CU1X(0,k,j))**2.d0  &
             + 0.25d0*dble(CU2X(0,k,j+1)+CU2X(0,k,j))**2.d0    &
             + 0.25d0*dble(CU3X(0,k+1,j)+CU3X(0,k,j))**2.d0)
        enddo
        enddo

! Save CUi
      do k=0,NZ+1
        do i=0,NX2P
          do j=0,NY+1
            CR1X(i,k,j)=CU1X(i,k,j)
            CR2X(i,k,j)=CU2X(i,k,j)
            CR3X(i,k,j)=CU3X(i,k,j)
          end do
        end do
      end do


      do k=0,NZV-1
          do j=0,NY+1
            CALL MPI_BCAST_COMPLEX(CR1X(0,k,j),1,1)
            CALL MPI_BCAST_COMPLEX(CR2X(0,k,j),1,1)
            CALL MPI_BCAST_COMPLEX(CR3X(0,k,j),1,1)       
          enddo
      enddo

! Convert to physical space
      

      CALL REAL_FOURIER_TRANS_U1 (.false.)
      CALL REAL_FOURIER_TRANS_U2 (.false.)
      CALL REAL_FOURIER_TRANS_U3 (.false.)



! Get the turbulent kinetic energy at each level
      do k=0,NZ+1 
        do j=0,NY+1
          urms(k,j)=0.d0
          vrms(k,j)=0.d0
          wrms(k,j)=0.d0
           do i=0,min(NXP,NXP_L) 
            urms(k,j)=urms(k,j)+(U1X(i,k,j)-dble(CR1X(0,k,j)))**2
            vrms(k,j)=vrms(k,j)+(U2X(i,k,j)-dble(CR2X(0,k,j)))**2
            wrms(k,j)=wrms(k,j)+(U3X(i,k,j)-dble(CR3X(0,k,j)))**2
           end do
        end do 
      end do


      CALL MPI_COMBINE_STATS(urms,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(vrms,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(wrms,NZ+2,NY+2)


! Get the bulk rms value

      if(rank .eq. 0) then

      do k=0,NZ+1
      do j=0,NY+1
        urms(k,j)=dsqrt(urms(k,j)/(dble(NX)))
        vrms(k,j)=dsqrt(vrms(k,j)/(dble(NX)))
        wrms(k,j)=dsqrt(wrms(k,j)/(dble(NX)))
      end do
      end do 
      
      do j=0,NY
       do k=0,NZ
! tke_mean defined at GY points
        tke_mean_old(k,j)=0.5d0*( urms(k,j)**2.d0         &
             + 0.25*(vrms(k,j) + vrms(k,j+1))**2.d0     &
             + 0.25*(wrms(k,j) + wrms(k+1,j))**2.d0 ) 
       end do
      end do

      endif

! C Convert velocity back to Fourier space
      CALL REAL_FOURIER_TRANS_U1 (.true.) 
      CALL REAL_FOURIER_TRANS_U2 (.true.)
      CALL REAL_FOURIER_TRANS_U3 (.true.)


     return
     end 

!-------------------------------------C---
subroutine mean_background_potential_duct
!-------------------------------------C--



! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps
! This subroutine should be called in SAVE_STATS_CHAN after computing
! plane averaged statistics, with the velocity in physical space, and 
! CRi containing the velocity in Fourier space
use ntypes
use Domain
use Grid
use Fft_var, only :  CIKXP
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use variable_stat
use mpi_var
implicit none

      integer i,j,k,n
      CHARACTER*31 file_energy

      if (rank .eq. 0) then
       write(6,*) 'Calculating Z_star', rank
      endif


      call z_star_cal

      return
      end

!-------------------------------------C---
      subroutine z_star_cal
!-------------------------------------C--

! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps
! This subroutine should be called in SAVE_STATS_CHAN after computing
! plane averaged statistics, with the velocity in physical space, and 
! CRi containing the velocity in Fourier space
use ntypes
use Domain
use Grid
use Fft_var, only :  CIKXP,KX2P
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use variable_stat
use mpi_var
use qsort_c_module
implicit none

INTERFACE
   FUNCTION Heavi_side (r)
     REAL(8)             ::  Heavi_side
     REAL(8), INTENT(IN) :: r
   END FUNCTION Heavi_side
END INTERFACE


      integer i,j,k,n,ii,jj,kk,NP_I, tot_ind
      CHARACTER*31 file_energy
!      real(r8)  ::  rho_temp, rho_trans, rho_diff, DV
      real(r8)  ::   rho_trans_all(1:NP), H_domian, DV!, sum_1(1:NP), sum_temp
      real(r8)                          :: rho_max, rho_min, sum_vol, z_star_pre,BIG_NUM
      real(r8)  ::   temp_var(1:10,1:2)     
 
      real(r8),allocatable,dimension(:,:)   ::  rho_1D, rho_1D_tot
      real(r8),allocatable,dimension(:)     ::  drho_1D_totdz_s, drho_1D_totdz_s_check
      parameter (BIG_NUM = 10.d0**8.0)

      

      if (rank .eq. 0) then
       write(6,*) 'Calculating Potential_energy_budget', rank
      endif

      area = 0.0d0
      do k=0,NZ+1
          area = area + LX*DZ(k)
      enddo

      H_domian = GYF(NY) - GYF(1)

      Z_star_bar (:,:) = 0.0d0 ; 

      tot_ind = (NZ+2)*(NY+2)*(NXP+1)


!Started  z* calculation
      allocate (rho_1D(1:tot_ind,1:2))

      rho_1D(:,1) = 100.d0 ;
      rho_1D(:,2) = -10.d0 ;

 !    H_domian = dble(rank)

      jj = 0
      do i=0,min(NXP,NXP_L)
        do k=0,NZ+1
          do j=0,NY+1
           jj =  jj + 1
           ii = tot_ind*rank + (NY+2)*(NZ+2)*i + (NY+2)*k + j +1           ! Calculate unique grid index 
           rho_1D(jj,1) = -alpha_T*THX(i,k,j,1) + THBAR(j,N_TH+1)          ! Storing the rho in a 1D array
           rho_1D(jj,2) =  dble(ii)  ! dble(tot_ind*rank + jj)             ! Storing the index for each grid to a 1D array
          enddo
         enddo
       enddo

      if (rank .eq. 0) then
       allocate (rho_1D_tot(1:tot_ind*NP,1:2))                            
       allocate (drho_1D_totdz_s(1:tot_ind*NP))
       allocate (drho_1D_totdz_s_check(1:tot_ind*NP))
       
       drho_1D_totdz_s_check(:) = 0.0d0
       drho_1D_totdz_s(:)       = 0.0d0
       rho_1D_tot(:,:)          = 0.0d0 
      endif

      CALL MPI_GATHER(rho_1D(:,1),tot_ind,MPI_DOUBLE_PRECISION,rho_1D_tot(:,1), &
         tot_ind,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_GATHER(rho_1D(:,2),tot_ind,MPI_DOUBLE_PRECISION,rho_1D_tot(:,2),  &
         tot_ind,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

      
  !    CALL MPI_GATHER(H_domian,1,MPI_DOUBLE_PRECISION,rho_trans_all,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

      if (rank .eq. 0) then

!       sorting for rho values : keep the grid index
        call QsortC(rho_1D_tot(:,1),rho_1D_tot(:,2))

        DO ii = 1,tot_ind*NP
         drho_1D_totdz_s_check(ii) = rho_1D_tot(ii,1) !!!!!!!!!!!!!!!!!!!!!!!1

         if ( rho_1D_tot(ii,1) .gt. 90.0d0 ) then
          kk = ii-1
          goto 212
         else

! steps to calculate dz*/drho
! keeping 1/dp(ii) information
 
          sum_vol = rho_1D_tot(ii+1,1)-rho_1D_tot(ii,1)
          if ( abs(sum_vol) .lt. 1.d0/BIG_NUM ) then
            if (sum_vol .lt. 0.0d0 ) then
             drho_1D_totdz_s(ii) = - BIG_NUM 
            else
             drho_1D_totdz_s(ii) =   BIG_NUM
            endif 
          else
           drho_1D_totdz_s(ii) = 1.0d0/sum_vol
          endif
         endif

        ENDDO

212      continue

        rho_1D_tot(:,1) = 0.0d0 
        sum_vol = 0.0d0 ;
        z_star_pre = 0.0d0 ;
          
        open(66,file='plane_energy/z_star_prof.dat',form='formatted', &
       status='unknown' )

         DO ii = kk, 1,-1 
           NP_I = floor(rho_1D_tot(ii,2)/tot_ind)
              i = floor((rho_1D_tot(ii,2)- (NP_I)*tot_ind)/((NZ+2)*(NY+2)))
              k = floor((rho_1D_tot(ii,2)- (NP_I)*tot_ind - i*((NZ+2)*(NY+2)))/(NY+2))
              j = floor(rho_1D_tot(ii,2)- (NP_I)*tot_ind - i*((NZ+2)*(NY+2)) - k*(NY+2) -1)

              DV  = DX(i)*DY(j)*DZ(k)
!             calculating z* at different ii level and  mapped back to physical grid(i,j,k)
              rho_1D_tot(int(rho_1D_tot(ii,2)),1) = (sum_vol + 0.50d0*DV)/area
              sum_vol =  sum_vol +  DV  
 
!             calculating dz*drho at different ii level : not yet mapped to physical grid(i,j,k) 
              drho_1D_totdz_s(ii) =  drho_1D_totdz_s(ii)*(z_star_pre-(sum_vol + 0.50d0*DV)/area)  

              z_star_pre =   (sum_vol + 0.50d0*DV)/area  !  storing the z* for (ii) level to calculate z*(ii)-z*(ii-1)
              write(66,565) real(ii),drho_1D_totdz_s_check(ii),(sum_vol + 0.50d0*DV)/area, drho_1D_totdz_s(ii)
         ENDDO
         close(66)
565      format(F14.2,3E18.9)
!        interpolating dz*drhp at last point of the array (here is kk)
         drho_1D_totdz_s(kk) = 2d0*drho_1D_totdz_s(kk-1) - drho_1D_totdz_s(kk-2)  
      endif

      call MPI_SCATTER(rho_1D_tot(:,1),tot_ind,MPI_DOUBLE_PRECISION, rho_1D(:,1),  &
                  tot_ind,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,IERROR)



      S1X(:,:,:) = 0.0d0
      S2X(:,:,:) = 0.0d0

       jj = 0
       do i=0,min(NXP,NXP_L)
        do k=0,NZ+1
          do j=0,NY+1
           jj =  jj + 1              
              S1X(I,K,J) = rho_1D(jj,1)
              S2X(I,K,J) = rho_1D(jj,1)      
          enddo
         enddo
        enddo 

! Done with z* calculation 
!      

        Z_star_bar (:,:) = 0.0d0 ; 

       do j=1,NY
       do k=1,NZ
        Z_star_bar(k,j)=0.
       do i=0,min(NXP,NXP_L)
        Z_star_bar(k,j)=Z_star_bar(k,j) + S1X(i,k,j)
       end do
        Z_star_bar(k,j)=Z_star_bar(k,j)/(float(NX))
      end do
      end do


       do k=1,NZ
        phi_b2(k)= 0.0d0 
        do i=0,min(NXP,NXP_L) 
!         phi_b2(k)= phi_b2(k) + S1X(i,k,1)*(-alpha_T*THX(i,k,2,1) + alpha_T*THX(i,k,1,1)  &
!            + THBAR(2,N_TH+1) - THBAR(1,N_TH+1))/(DY(2))
 
          phi_b2(k)= phi_b2(k) + S1X(i,k,1)*(-alpha_T)*( 2.0*(THX(i,k,2,1) - THX(i,k,1,1))/DY(2) &
                                - (THX(i,k,3,1) - THX(i,k,2,1))/DY(3))
       end do
        phi_b2(k)= -(NU/PR(1))*Gravity*phi_b2(k)/dble(NX)
      end do


!    Calculate the Eb (background potential energy)

     

      do k=0,NZ+1 
        do j=0,NY+1
          back_potential(k,j)=0.d0
           do i=0,min(NXP,NXP_L) 
            back_potential(k,j)= back_potential(k,j) + Gravity*(-alpha_T*THX(i,k,j,1))*S1X(i,k,j)
           end do
          back_potential(k,j)=back_potential(k,j)/(float(NX))
        end do 
      end do


      CALL MPI_COMBINE_STATS(back_potential,NZ+2,NY+2)

     IF (TIME_STEP .eq. 0) THEN        
       tke_1(:,:) = 0.0d0
      ELSE
      do j=0,NY+1
       do k=0,NZ+1
        tke_1(k,j)=(back_potential(k,j)-back_potential_old(k,j)) &
             /(TIME-TIME_old_pot)
       end do
      end do
      ENDIF

      back_potential_old = back_potential
      TIME_old_pot = TIME


      if (rank .eq. 0) then
        rho_1D_tot(:,1) = 0.0d0  ! reuse z* array for storing dz*drho values according to physical grid(i,j,k) 
        DO ii = 1,kk
        rho_1D_tot(int(rho_1D_tot(ii,2)),1) = drho_1D_totdz_s(ii)
        ENDDO
      endif

      call MPI_SCATTER(rho_1D_tot(:,1),tot_ind,MPI_DOUBLE_PRECISION, rho_1D(:,1),  &
                  tot_ind,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,IERROR)


! Compute the background  dissipation rate, epsilon=-kappa*g*<z*(d2rho/dx_i2)>
      do j=0,NY+1
       do k=0,NZ+1
        tke_2(k,j)=0.d0
       end do
      end do

! Store d2rho/dx2 in CS1
      do j=1,NY
      do k=1,NZ
      do i=0,NX2P
        CS1X(i,k,j)=-KX2P(I)*(-alpha_T)*CFTHX(i,k,j,1) ! WE HAVE ALSO STORED CTHX IN CFTHX
      end do
      end do
      end do
    
! Convert to physical space

      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
      ENDDO
      S1Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)


      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
         tke_2(k,j)=tke_2(k,j) + S2X(I,K,J)*S1X(i,k,j)
      end do
      end do
      end do

! Compute d2<..>/dy2
      

      do j=2,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
       S1X(i,k,j) = -alpha_T*(((THX(I,K,J+1,1) - THX(I,K,J,1)) / DY(J+1)            &
              -(THX(I,K,J,1)   - THX(I,K,J-1,1)) / DY(J)) /DYF(J))
      end do
      end do
      end do


      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        tke_2(k,j)=tke_2(k,j) + S2X(I,K,J)*S1X(i,k,j)
!         tke_2(k,j) = S1X(i,k,j) 
      end do
      end do
      end do

      
      do k=1,NZ
         tke_2(k,1) = 2.0d0*tke_2(k,2) - tke_2(k,3)
      end do
     
! Compute d2<..>/dz2 


      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        S1X(i,k,j) = -alpha_T*(((THX(I,K+1,J,1) - THX(I,K,J,1)) / DZ(K+1)             &
                -(THX(I,K,J,1)   - THX(I,K-1,J,1)) / DZ(K)) /DZF(K))
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
         tke_2(k,j)=tke_2(k,j) + S2X(I,K,J)*S1X(i,k,j)
      end do
      end do
      end do

      CALL MPI_COMBINE_STATS(tke_2,NZ+2,NY+2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=1,NY
      do k=1,NZ
        tke_2(k,j)=(NU/PR(1))*Gravity*tke_2(k,j)/float(NX)
      end do
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      S2X(:,:,:) = 0.0d0

       jj = 0
       do i=0,min(NXP,NXP_L)
        do k=0,NZ+1
          do j=0,NY+1
           jj =  jj + 1              
              S2X(I,K,J) = rho_1D(jj,1)
          enddo
         enddo
        enddo

! Compute the irreversible mixing/APE dissipation rate, epsilon=-kappa*g*<dz*drho*(drho/dx_i)^2.0>
      do j=0,NY+1
       do k=0,NZ+1
        epsilon_rho(k,j)=0.d0
       end do
      end do
   
! Store du/dx in CS1
      do j=1,NY
      do k=1,NZ
      do i=0,NX2P
        CS1X(i,k,j)=CIKXP(i)*(-alpha_T)*CFTHX(i,k,j,1) ! WE HAVE ALSO STORED CTHX IN CFTHX
      end do
      end do
      end do
! Convert to physical space

      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
      ENDDO
      S1Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)


      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
         epsilon_rho(k,j)=epsilon_rho(k,j) + S2X(I,K,J)*S1X(i,k,j)**2.0
      end do
      end do
      end do

     
      
! Compute d<..>/dy note remove mean


      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
       S1X(i,k,j) = -alpha_T*(THX(i,k,j+1,1)-THX(i,k,j-1,1))/(2.0d0*DYF(j)) &
                    + (THBAR(j+1,N_TH+1)-THBAR(j-1,N_TH+1))/(2.*DYF(j))
      end do
      end do
      end do


      do j=2,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        epsilon_rho(k,j)=epsilon_rho(k,j)+S2X(I,K,J)*S1X(i,k,j)**2.0
      end do
      end do
      end do
      
  
      
! Compute d<..>/dz 


      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        S1X(i,k,j) = -alpha_T*(THX(i,k+1,j,1)-THX(i,k-1,j,1))/(2.0d0*DZF(k))
      end do
      end do
      end do    

      do j=2,NY
      do k=1,NZ-2
      do i=0,min(NXP,NXP_L)
         epsilon_rho(k,j)=epsilon_rho(k,j)+ S2X(I,K,J)*S1X(i,k,j)**2.0
      end do
      end do
      end do

      do k=1,NZ
        epsilon_rho(k,1) = 2.0*epsilon_rho(k,2) - epsilon_rho(k,3)
      end do

      do j=1,NY
        epsilon_rho(NZ-1,j) = 2.0*epsilon_rho(NZ-2,j) - epsilon_rho(NZ-3,j)
      end do

      do j=1,NY
        epsilon_rho(NZ,j) = 2.0*epsilon_rho(NZ-1,j) - epsilon_rho(NZ-2,j)
      end do

!!!!!! NEED ADD FOR ALL PROCESSORS!!!!!!!!!!!!
      CALL MPI_COMBINE_STATS(epsilon_rho,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(phi_b2,NZ+2,1)
      CALL MPI_COMBINE_STATS(Z_star_bar,NZ+2,NY+2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=1,NY
      do k=1,NZ
        epsilon_rho(k,j)=-(NU/PR(1))*Gravity*epsilon_rho(k,j)/float(NX)
      end do
      end do

      
      do k=1,NZ
        Z_star_bar(k,j)=Z_star_bar(k,j)/(float(NX))
      end do
      

!          Calculating area integration       
      
      IF (TIME_STEP .eq. 0) THEN        
       epsilon_rho(:,:)    = 0.d0
       back_potential(:,:) = 0.0d0
       Z_star_bar(:,:)     = 0.0d0
       tke_2(:,:)     = 0.0d0
      ENDIF

      area      = 0.0d0 ;
      sum_Prod  = 0.0d0 ;
      sum_dissp = 0.0d0 ;
      sum_advec = 0.0d0 ;
      
      do k=1,NZ
       area  = area + DZ(k)
       sum_Prod  = sum_Prod  + phi_b2(k)*DZ(k)
       sum_advec = sum_advec + Z_star_bar(k,1)*(-alpha_T)*dble(CRTHX(0,k,2,1)-CRTHX(0,k,1,1))/(DY(2))*DZ(k)
      end do

      sum_Prod = sum_Prod/area ;
      sum_advec= -(NU/PR(1))*Gravity*sum_advec/area ;

      area = 0.0d0 ;
      sum_dissp = 0.0d0 ;
      sum_dtkdt = 0.0d0 ;
      sum_buoyF = 0.0d0 ;
      sum_vis   = 0.0d0 ;
 
      do j=1,NY
      do k=1,NZ
       area  = area + DY(j)*DZ(k)
       sum_dissp  = sum_dissp  + epsilon_rho(k,j)*DY(j)*DZ(k)
       sum_dtkdt  = sum_dtkdt  + tke_1(k,j)*DY(j)*DZ(k)
       sum_vis    = sum_vis    + tke_2(k,j)*DY(j)*DZ(k)
       sum_buoyF  = sum_buoyF  + back_potential(k,j)*DY(j)*DZ(k)
      enddo
      enddo

       sum_dissp  = sum_dissp/area
       sum_dtkdt   = sum_dtkdt/area
       sum_buoyF  = sum_buoyF/area
       sum_vis    = sum_vis/area


       if (rank .eq. 0) then
        deallocate (rho_1D_tot,drho_1D_totdz_s,drho_1D_totdz_s_check )
       open(66,file='plane_energy/background_budget.dat',form='formatted', &
       status='unknown' , position='append' )

       write(66,566) TIME/(60.0d0*60.0d0),sum_buoyF, sum_dtkdt, sum_Prod, sum_dissp, &
                     sum_vis,sum_advec
566    format(f10.5,6E18.9)
       close(66)

        write(6,*) 'Done with Potential_energy_budget', rank
       endif


       deallocate (rho_1D )

       return
       end

! function heaviside,x
! y = (x gt 0)*1D + ((x lt 0))*0D + 0.5D*(x eq 0)
! return,y
! end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION Heavi_side(r)

use ntypes
use Domain

IMPLICIT NONE

REAL(r8) :: Heavi_side
REAL(r8), INTENT(IN) :: r



Heavi_side = (r/abs(r)+1.d0)/2.d0
if( abs(r) == 0.0d0) Heavi_side = 0.50d0

END FUNCTION Heavi_side 
