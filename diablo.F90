!C******************************************************************************|
!C diablo.F90 -> DNS In A Box, Laptop Optimized                       VERSION 0.9
!C


!C----|--.---------.---------.---------.---------.---------.---------.-|-------|
      PROGRAM DIABLO
  !    INCLUDE 'header_duct'
   
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use omp_lib
use mpi_var, only: rank, IERROR
use variable_stat, only : energy_mean_old

implicit none

      INTEGER N,j,k
      REAL*8  int_time 
      

      WRITE(6,*) 
      WRITE(6,*) '             ****** WELCOME TO DIABLO ******'
      WRITE(6,*)

      CALL DATE_AND_TIME (VALUES=TIME_ARRAY)
      START_TIME=TIME_ARRAY(5)*3600+TIME_ARRAY(6)*60  &
        +TIME_ARRAY(7)+TIME_ARRAY(8)*0.001

      CALL INITIALIZE
! Initialize START_TIME for run timing


!C A flag to determine if we are considering the first time-step
      FIRST_TIME=.TRUE. 
     
        
      int_time = TIME
      
      DO TIME_STEP = TIME_STEP+1, TIME_STEP+N_TIME_STEPS
        IF (rank .eq. 0) then
        WRITE(6,*) 'Now beginning TIME_STEP = ',TIME_STEP, wtime
        ENDIF

!        DO N=1,N_TH 
         IF (INT_TREAT ) THEN
           If ( (TIME-int_time) .LT. 8.0*PI) then
             RI_TAU = RI_FINAL/5.0d0 + (RI_FINAL-RI_FINAL/5.0d0)*(TIME-int_time)/(8.0d0*PI)
             IF (rank .eq. 0) then
              write(6,*) 'Gravity is ', RI_TAU
             ENDIF
           ELSE
             RI_TAU = RI_FINAL
           ENDIF    
         ELSE
           RI_TAU = RI_FINAL
         ENDIF   
!     
!         ENDDO 
 
        IF (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0) THEN
! pre calculation of mean_energy and tke for getting right transient terms
!         call   precalculation_for_budget
        ENDIF

        wtime =  omp_get_wtime ( )
        DO RK_STEP=1,3
          IF (NUM_PER_DIR.EQ.3) THEN
!c            IF (TIME_AD_METH.EQ.1) CALL RK_PER_1
!c            IF (TIME_AD_METH.EQ.2) CALL RK_PER_2            
          ELSEIF (NUM_PER_DIR.EQ.2) THEN
!c            IF (TIME_AD_METH.EQ.1) CALL RK_CHAN_1
!c            IF (TIME_AD_METH.EQ.2) CALL RK_CHAN_2            
!c            IF (TIME_AD_METH.EQ.3) CALL RK_CHAN_3            
          ELSEIF (NUM_PER_DIR.EQ.1) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_DUCT_1
            IF (TIME_AD_METH.EQ.2) CALL RK_DUCT_2
            IF (TIME_AD_METH.EQ.3) CALL RK_DUCT_3            
          ELSEIF (NUM_PER_DIR.EQ.0) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CAV_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CAV_2            
          END IF
        ENDDO 


        

        TIME=TIME+DELTA_T
        FIRST_TIME=.FALSE.

        wtime = -(wtime - omp_get_wtime ( ))
! Save statistics to an output file

!c        IF (MOD(TIME_STEP,100).EQ.0) THEN
!c          IF (USE_MPI) THEN
!c           CALL GHOST_CHAN_MPI
!c          END IF
!c          CALL POISSON_P_DUCT
!c        ENDIF        

        IF (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0) THEN
            CALL SAVE_STATS(.FALSE.)
        END IF
! Save the flow to a restart file
        IF (MOD(TIME_STEP,SAVE_FLOW_INT).EQ.0) THEN
          CALL SAVE_FLOW(.FALSE.)
        END IF
! Filter the scalar field
!c        DO N=1,N_TH
!c          IF (FILTER_TH(N)
!c     &       .AND.(MOD(TIME_STEP,FILTER_INT(N)).EQ.0)) THEN
!c          write(*,*) 'Filtering...'
!c          CALL FILTER(N)
!c          END IF 
!c        END DO
       
! If we are considering a Near Wall Model, it may be necessary to
! filter the velocity field.  If so, use the following:
!c          IF (FILTER_VEL
!c     &       .AND.(MOD(TIME_STEP,FILTER_VEL_INT).EQ.0)) THEN
!c            write(*,*) 'Filtering Velocity'
!c           CALL APPLY_FILTER_VEL(CU1,1,NY)
!c           CALL APPLY_FILTER_VEL(CU2,1,NY)
!c           CALL APPLY_FILTER_VEL(CU3,1,NY)
!c        END IF

!c        IF (MOD(TIME_STEP,200).EQ.0) THEN
!c          IF (USE_MPI) THEN
!c           CALL GHOST_CHAN_MPI
!c          END IF
!c          CALL POISSON_P_DUCT
!c         ENDIF
 
      END DO

! Calculate and display the runtime for the simulation
      CALL DATE_AND_TIME (VALUES=TIME_ARRAY)
      END_TIME=TIME_ARRAY(5)*3600+TIME_ARRAY(6)*60  &
        +TIME_ARRAY(7)+TIME_ARRAY(8)*0.001
      WRITE(*,*) 'Elapsed Time (sec): ',end_time-start_time
      WRITE(*,*) 'Seconds per Iteration: '      &
          ,(end_time-start_time)/N_TIME_STEPS

      TIME_STEP=TIME_STEP-1
      CALL SAVE_FLOW(.TRUE.)
      CALL SAVE_STATS(.TRUE.)
      if (rank .eq. 0 ) then
      WRITE(6,*)
      WRITE(6,*) '        ****** Hello world!  Have a nice day! ******'
      WRITE(6,*)
      endif

      call MPI_FINALIZE(IERROR)
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INITIALIZE
!      INCLUDE 'header_duct'

use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
!use les_chan_var , only : C_DYN     
use mpi_var 
implicit none


      REAL(r8)  ::   VERSION, CURRENT_VERSION, DAMP_FACT,local_tmp
      logical   ::   RESET_TIME
      INTEGER   ::   I,II, J, K, N
 
        OPEN (11,file='signal_plane_data',form='formatted',status='old')
        READ(11,*)
        READ(11,*)
        READ(11,*)
        READ(11,*) ipn, jpn1, jpn2, kpn1, kpn2 
       close(11)

      OPEN (11,file='input.dat',form='formatted',status='old')      

!      WRITE(6,*) 'Note that this code is distributed under the ', &
!                'GNU General Public License.'
!      WRITE(6,*) 'No warranty is expressed or implied.'
!      WRITE(6,*)


!ALLOCATE VARIABLES

      call   allocate_var
      call   allocate_temps 
      

!C Read input file.
!C   (Note - if you change the following section of code, update the
!C    CURRENT_VERSION number to make obsolete previous input files!)

      CURRENT_VERSION=0.90d0
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*) VERSION
      IF (VERSION .NE. CURRENT_VERSION) STOP 'Wrong input data format. '
      READ(11,*)
      READ(11,*) USE_MPI
      READ(11,*)
      READ(11,*) NU, LX, LY, LZ
      READ(11,*)
      READ(11,*) N_TIME_STEPS, DELTA_T, RESET_TIME, VARIABLE_DT, CFL  &
                 ,DFN, UPDATE_DT
      READ(11,*)
      READ(11,*) NUM_PER_DIR, TIME_AD_METH, LES, LES_MODEL_TYPE, LES_MODEL_TYPE_TH
      READ(11,*)
      READ(11,*) VERBOSITY, SAVE_FLOW_INT, SAVE_STATS_INT, MOVIE
      READ(11,*)
      READ(11,*) CREATE_NEW_FLOW, IC_TYPE, KICK, INT_TREAT
      READ(11,*)
      READ(11,*) F_TYPE, UBULK0, PX0, OMEGA0, AMP_OMEGA0, ANG_BETA
      IF (NUM_PER_DIR.eq.3) THEN
      READ(11,*)
      READ(11,*) U_BC_XMIN, U_BC_XMIN_C1, U_BC_XMIN_C2, U_BC_XMIN_C3
      READ(11,*) 
      READ(11,*) V_BC_XMIN, V_BC_XMIN_C1, V_BC_XMIN_C2, V_BC_XMIN_C3
      READ(11,*)
      READ(11,*) W_BC_XMIN, W_BC_XMIN_C1, W_BC_XMIN_C2, W_BC_XMIN_C3
      READ(11,*)
      READ(11,*) U_BC_XMAX, U_BC_XMAX_C1, U_BC_XMAX_C2, U_BC_XMAX_C3
      READ(11,*)
      READ(11,*) V_BC_XMAX, V_BC_XMAX_C1, V_BC_XMAX_C2, V_BC_XMAX_C3
      READ(11,*)
      READ(11,*) W_BC_XMAX, W_BC_XMAX_C1, W_BC_XMAX_C2, W_BC_XMAX_C3
      END IF
      IF (NUM_PER_DIR.gt.0) THEN
      READ(11,*)
      READ(11,*) U_BC_YMIN, U_BC_YMIN_C1, U_BC_YMIN_C2, U_BC_YMIN_C3
      READ(11,*) 
      READ(11,*) V_BC_YMIN, V_BC_YMIN_C1, V_BC_YMIN_C2, V_BC_YMIN_C3
      READ(11,*)
      READ(11,*) W_BC_YMIN, W_BC_YMIN_C1, W_BC_YMIN_C2, W_BC_YMIN_C3
      READ(11,*)
      READ(11,*) U_BC_YMAX, U_BC_YMAX_C1, U_BC_YMAX_C2, U_BC_YMAX_C3
      READ(11,*)
      READ(11,*) V_BC_YMAX, V_BC_YMAX_C1, V_BC_YMAX_C2, V_BC_YMAX_C3
      READ(11,*)
      READ(11,*) W_BC_YMAX, W_BC_YMAX_C1, W_BC_YMAX_C2, W_BC_YMAX_C3
      END IF
      IF (NUM_PER_DIR.lt.2) THEN
      READ(11,*)
      READ(11,*) U_BC_ZMIN, U_BC_ZMIN_C1, U_BC_ZMIN_C2, U_BC_ZMIN_C3
      READ(11,*) 
      READ(11,*) V_BC_ZMIN, V_BC_ZMIN_C1, V_BC_ZMIN_C2, V_BC_ZMIN_C3
      READ(11,*)
      READ(11,*) W_BC_ZMIN, W_BC_ZMIN_C1, W_BC_ZMIN_C2, W_BC_ZMIN_C3
      READ(11,*)
      READ(11,*) U_BC_ZMAX, U_BC_ZMAX_C1, U_BC_ZMAX_C2, U_BC_ZMAX_C3
      READ(11,*)
      READ(11,*) V_BC_ZMAX, V_BC_ZMAX_C1, V_BC_ZMAX_C2, V_BC_ZMAX_C3
      READ(11,*)
      READ(11,*) W_BC_ZMAX, W_BC_ZMAX_C1, W_BC_ZMAX_C2, W_BC_ZMAX_C3
      END IF
! Read Stochastic forcing parameters
      READ(11,*)
      READ(11,*) STOCHASTIC_FORCING
! Filter for the velocity field
      READ(11,*)
      READ(11,*) FILTER_VEL, FILTER_VEL_INT
      
      READ(11,*)
! Read in the parameters for the N_TH scalars
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) CREATE_NEW_TH(N)
        READ(11,*)
        READ(11,*) FILTER_TH(N), FILTER_INT(N)
        READ(11,*)
        READ(11,*) RI_TAU(N), PR(N), BACKGROUND_GRAD(N), DEV_BACK_TH
        READ(11,*)
        READ(11,*) TH_BC_YMIN(N),TH_BC_YMIN_C1(N),TH_BC_YMIN_C2(N)  &
                  ,TH_BC_YMIN_C3(N)  
        READ(11,*)
        READ(11,*) TH_BC_YMAX(N),TH_BC_YMAX_C1(N),TH_BC_YMAX_C2(N)  &
                  ,TH_BC_YMAX_C3(N)
        IF (NUM_PER_DIR.lt.2) THEN
         READ(11,*)
         READ(11,*) TH_BC_ZMIN(N),TH_BC_ZMIN_C1(N),TH_BC_ZMIN_C2(N) &
                  ,TH_BC_ZMIN_C3(N)
         READ(11,*)
         READ(11,*) TH_BC_ZMAX(N),TH_BC_ZMAX_C1(N),TH_BC_ZMAX_C2(N) &
                  ,TH_BC_ZMAX_C3(N)
        ENDIF
      END DO

!C If we are using MPI, then Initialize the MPI Variables
      IF (USE_MPI) THEN
        CALL INT_MPI
      END IF

      IF (LES)THEN
       call allocation_les_var
       call   allocate_les_tmp
       NU_T(:,:,:)      = 0.0d0
       KAPPA_T(:,:,:,:) = 0.0d0
      ENDIF

!       DO N=1,N_TH
!        RI_TAU(N) = RI_FINAL
!       ENDDO
  

      DELTA_T_in = DELTA_T
      ANG_BETA = ANG_BETA*3.14159265/180.0 
      Q_H0 = 1.d0
      H0    = 10.0*LY     
      In_H0 = 1.d0/H0 
      count_data = 0
     
        

 
      if (N_TH .gt. 0) then
       rho_0   = 1000.0d0 ;
       Gravity = 10.0
       alpha_T = 2.0*10**(-4.0)
       RI_TAU(1) = RI_TAU(1)*Gravity*alpha_T
       RI_FINAL  = RI_TAU(1)

       theta_0   = 30.0d0 ;
       TH_1      = -20.0d0 ;
       TH_2      = 10.0d0 ;
       Lh_s      = 0.07   ;
       Loc_s     = 0.50d0/3.0d0 ;

       f_0       = -2.0d0 
       beta_f    = 0.6

       IBM       = .true.
      endif

     IF (rank .eq. 0) then
!C Initialize grid
      WRITE(6,*) 'Grid size: NX =',NX,', NY =',NY,', NZ =',NZ,'.'
      WRITE(6,*) 'Grid size sup : NX =',NXV,', NY =',NY,', NZ =',NZV,'.'
      WRITE(6,*) 'Domain size: LX =',LX,', LY =',LY,', LZ =',LZ,'.' 
      WRITE(6,*)'NU',NU,'DELTA_T',DELTA_T,'Delta_s',sqrt(2.0*NU/OMEGA0) 
      WRITE(6,*)'CFL',CFL, 'Diffu No' , DFN
      WRITE(6,*) 'Variable delta_t applited',  VARIABLE_DT,U_BC_ZMAX_C1
      WRITE(6,*) 'NX= ', NX, ' NXV= ', NXV, 'NX2P = ', NX2P  
       
      

      DO N=1,N_TH
        WRITE(6,*) 'Scalar number: ',N
        WRITE(6,*) ' Final Richardson number: ',RI_TAU(N)
        WRITE(6,*) '  Prandlt number: ',PR(N)
      END DO
     endif
      

      NXM=NX-1
      NYM=NY-1
      NZM=NZ-1

      NXP_L=NX-(NXP+1)*rank-1
      NX2P_L=NKX+1-(NX2P+1)*rank-1
      write(6,*)'Lower bound rank', rank, min(NXP_L,NXP),min(NX2P_L,NX2P)

      IF (NUM_PER_DIR .GT. 0) THEN

         if (rank ==0) then
          WRITE (6,*) 'Fourier in X'
         endif
         
         DO I=0,NX
           GX(I)=(I*LX)/NX
           DX(I)=LX/NX
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
         END DO
      ELSE
         WRITE (6,*) 'Finite-difference in X'
         OPEN (30,file='xgrid.txt',form='formatted',status='old')
         READ (30,*) NX_T
!C Check to make sure that grid file is the correct dimensions
         IF (NX_T.ne.NX) THEN
           WRITE(6,*) 'NX, NX_T',NX,NX_T
           STOP 'Error: xgrid.txt wrong dimensions'
         END IF
         DO I=1,NX+1
           READ(30,*) GX(I)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
         END DO 
         DO I=1,NX
           READ(30,*) GXF(I)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GXF(',I,') = ',GXF(I)
         END DO
!C Define ghost cells, if needed for this grid...
         GXF(0)=2.d0*GXF(1)-GXF(2)
         GXF(NX+1)=2.d0*GXF(NX)-GXF(NXM)
         GX(0)=2.d0*GX(1)-GX(2)
!C Define the grid spacings 
         DO I=1,NX+1
           DX(I)=(GXF(I)-GXF(I-1))
         END DO
         DO I=1,NX
           DXF(I)=(GX(I+1)-GX(I))
         END DO
         CLOSE(30)
      END IF

      IF (NUM_PER_DIR .GT. 1) THEN
         WRITE (6,*) 'Fourier in Z'
         DO K=0,NZ
           GZ(K)=(K*LZ)/NZ
           DZ(K)=LZ/NZ
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ(',K,') = ',GZ(K)
         END DO
      ELSE
         if (rank ==0) then
          write(6,*) 'Finite-difference in Z'
!         endif

         OPEN (30,file='zgrid.txt',form='formatted',status='old')
         READ (30,*) NZ_T
!C Check to make sure that grid file is the correct dimensions
         IF (NZ_T.ne.NZ) THEN
           WRITE(6,*) 'NZ, NZ_T',NZ,NZ_T
           STOP 'Error: zgrid.txt wrong dimensions'
         END IF
         DO K=1,NZ+1
           READ(30,*) GZ(k)
!           write(6,*) GZ(k), rank
!           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ(',K,') = ',GZ(K)
         END DO 
         DO K=1,NZ
           READ(30,*) GZF(k)
!           write(6,*) GZF(k)
!           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZF(',K,') = ',GZF(K)
         END DO
         CLOSE(30)
         endif

         do k=1,NZ+1
           CALL MPI_BCAST_REAL(GZ(k),1,1)
         enddo 
         do k=1,NZ
           CALL MPI_BCAST_REAL(GZF(k),1,1)
         enddo

!C Define ghost cells, if needed for this grid...
         GZF(0)=2.d0*GZF(1)-GZF(2)
         GZF(NZ+1)=2.d0*GZF(NZ)-GZF(NZM)
         GZ(0)=2.d0*GZ(1)-GZ(2)
!C Define grid spacing 
         DO K=1,NZ+1
           DZ(K)=(GZF(K)-GZF(K-1))
         END DO
         DO K=1,NZ
           DZF(K)=(GZ(K+1)-GZ(K))
!           write(400,*) DZF(K)
         END DO
         DZ(0)=DZ(1)
         DZF(NZ+1)=DZF(NZ) 
!         CLOSE(30)
      END IF
      
      IF (NUM_PER_DIR .GT. 2) THEN
         WRITE (6,*) 'Fourier in Y'
         DO J=0,NY
           GY(J)=(J*LY)/NY
           DY(J)=LY/NY
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GY(',J,') = ',GY(J)
         END DO
      ELSE

         if (rank ==0) then
          write(6,*) 'Finite-difference in Y'
!         endif
    
         OPEN (30,file='./ygrid.txt',form='formatted',status='old')
         READ (30,*) NY_T
!C Check to make sure that grid file is the correct dimensions
         IF (NY_T.ne.NY) THEN
           WRITE(6,*) 'NY, NY_T',NY,NY_T
           STOP 'Error: ygrid.txt wrong dimensions'
         END IF
         DO J=1,NY+1
           READ(30,*) GY(j)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GY(',J,') = ',GY(J)
         END DO
         DO J=1,NY
           READ(30,*) GYF(j)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GYF(',J,') = ',GYF(J)
         END DO
         CLOSE(30)
         endif

         do j=1,NY+1
           CALL MPI_BCAST_REAL(GY(j),1,1)
         enddo
         do j=1,NY
           CALL MPI_BCAST_REAL(GYF(j),1,1)
         enddo
          
!          IF (USE_MPI) THEN
!            CALL GHOST_GRID_MPI
!          ELSE
!C Define ghost cells
           GYF(0)=2.d0*GYF(1)-GYF(2)
           GYF(NY+1)=2.d0*GYF(NY)-GYF(NYM)
           GY(0)=2.d0*GY(1)-GY(2)
!         END IF

!C Define grid spacing
         DO J=1,NY+1
           DY(J)=(GYF(J)-GYF(J-1))
         END DO
         DO J=1,NY
           DYF(J)=(GY(J+1)-GY(J))
         END DO
         DY(0)=DY(1)
         DYF(NY+1)=DYF(NY)
      END IF


!    Done with grid information
      if (rank ==0) then
       write(*,*) 'Grid reading is done'
      endif
 
      DO J=0,NY+1
       DO K=0,NZ+1
        U2_BAR(K,J)=0.d0
        U1_BAR(J)=0.d0
        U3_BAR(K,J)=0.d0
       ENDDO 
      ENDDO
!C   Boundary condition prescribed by the Blasius for IC_TYPE = 5
      IF ( IC_TYPE .EQ. 5) THEN
      OPEN (80,file='./int_vel_xy.txt',form='formatted',status='old')     
       DO J=1,NY
        DO K=1,NZ
        read(80,*) U3_BAR(K,J), U2_BAR(K,J)
       ENDDO
       ENDDO

       DO J=1,NY 
        DO K=2,NZ
          S1X(0,K,J) = 0.5*(U3_BAR(K,J) + U3_BAR(K-1,J))
        ENDDO
       ENDDO

       DO J=1,NY
        DO K=2,NZ
          U3_BAR(K,J) =S1X(0,K,J)
        ENDDO
       ENDDO 

      ENDIF         


!     stop
     
     call allocation_u (.true.)
     call allocation_v (.true.)
     call allocation_w (.true.)  
     call allocation_p (.true.)

     call allocation_R1 (.true.) 
     call allocation_R2 (.true.)
     call allocation_R3 (.true.) 
     call allocation_F1 (.true.)
     call allocation_F2 (.true.)
     call allocation_F3 (.true.)
 


      IF(N_TH .GT. 0) THEN
       call allocation_th (.true.)
       call allocation_Rth (.true.)
       call allocation_Fth (.true.)
      ENDIF


! C Initialize storage arrays.
      DO K=0,NZV-1
        DO I=0,NXP
          DO J=0,NY+1
            U1X(I,K,J)=0.
            U3X(I,K,J)=0.
            U2X(I,K,J)=0.
            PX (I,K,J)=0.
            R1X(I,K,J)=0.
            R2X(I,K,J)=0.
            R3X(I,K,J)=0.
            F1X(I,K,J)=0.
            F2X(I,K,J)=0.
            F3X(I,K,J)=0.
            DO N=1,N_TH
             THX(I,K,J,N)=0.
             RTHX(I,K,J,N)=0. 
             FTHX(I,K,J,N)=0.
            END DO
            
! Array for LES subgrid model
            IF (LES) THEN
   !          NU_T(I,K,J)=0.
             DO N=1,N_TH
   !            KAPPA_T(I,K,J,N)=0.
             END DO
            ENDIF
          END DO
        END DO
      END DO

!C Initialize FFT package (includes defining the wavenumber vectors).
      CALL INIT_FFT
     
      IF (IBM) THEN
       DO I=0,NXP
        H_hill(I) = 0.020d0*(1.0d0 - cos(6.0d0*PI*DX(I)*dble(I+RANK*(NXP+1))/LX))
       ENDDO

       DO I=0,NXP
        ii = 0;
        DO J=0,NY
         IF ( (GYF(J) >= H_hill(I)) .AND. (ii==0)) THEN
          hill_ind(I) = J + 1;
          ii = 1 ;
         ENDIF
        ENDDO
        write(6,*)rank, hill_ind(I),H_hill(I)
       ENDDO

       call MPI_INST_PROF (dble(hill_ind(0:NXP)),H_hill_tot)
       hill_ind_tot=int(H_hill_tot)
       call MPI_INST_PROF (H_hill(0:NXP),H_hill_tot)

       IF(RANK==0)THEN
        open(202,file='hill_prof.dat',form='formatted',status='unknown')
        DO I=0,(NXP+1)*NP-1
         write(202,*) DX(1)*dble(I),H_hill_tot(I),hill_ind_tot(I)
        ENDDO
       ENDIF
      ENDIF


 
!C Initialize RKW3 parameters.
      H_BAR(1)=DELTA_T*(8.0/15.0)
      H_BAR(2)=DELTA_T*(2.0/15.0)
      H_BAR(3)=DELTA_T*(5.0/15.0)
      BETA_BAR(1)=1.0
      BETA_BAR(2)=25.0/8.0
      BETA_BAR(3)=9.0/4.0
      ZETA_BAR(1)=0.0
      ZETA_BAR(2)=-17.0/8.0
      ZETA_BAR(3)=-5.0/4.0
      
!C Initialize case-specific packages.
      IF (NUM_PER_DIR.EQ.3) THEN
!C        CALL INIT_PER
      ELSEIF (NUM_PER_DIR.EQ.2) THEN 
        CALL INIT_CHAN
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        CALL INIT_DUCT
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        CALL INIT_CAV
      END IF

!C Initialize values for reading of scalars
      IF (CREATE_NEW_TH(1)) THEN
       IF (NUM_PER_DIR.EQ.1) THEN
         CALL CREATE_TH_DUCT
       ELSE IF (NUM_PER_DIR.EQ.2) THEN
        CALL CREATE_TH_CHAN
       ELSE IF (NUM_PER_DIR.EQ.3) THEN
! C        CALL CREATE_TH_PER
       END IF
      ENDIF
      


      IF (CREATE_NEW_FLOW) THEN
        IF (NUM_PER_DIR.EQ.3) THEN
!C          CALL CREATE_FLOW_PER
        ELSEIF (NUM_PER_DIR.EQ.2) THEN
          CALL CREATE_FLOW_CHAN
        ELSEIF (NUM_PER_DIR.EQ.1) THEN
          CALL CREATE_FLOW_DUCT
        ELSEIF (NUM_PER_DIR.EQ.0) THEN
          CALL CREATE_FLOW_CAV
        END IF

!        write(*,*) 'A new flowfield has been created'
!c        CALL SAVE_STATS(.FALSE.)
        CALL SAVE_FLOW(.FALSE.)
      ELSE

! C     Convert velocity back to Fourier space
                CALL REAL_FOURIER_TRANS_U1 (.true.)
                CALL REAL_FOURIER_TRANS_U2 (.true.)
                CALL REAL_FOURIER_TRANS_U3 (.true.)
                CALL REAL_FOURIER_TRANS_P (.true.)
        
                CALL REAL_FOURIER_TRANS_R1 (.true.)
                CALL REAL_FOURIER_TRANS_R2 (.true.)
                CALL REAL_FOURIER_TRANS_R3 (.true.)
        
                call allocation_F1 (.false.)
                call allocation_F2 (.false.)
                call allocation_F3 (.false.)

        CU1X(:,:,:) = (0.d0,0.d0)
        CU2X(:,:,:) = (0.d0,0.d0)
        CU3X(:,:,:) = (0.d0,0.d0)
        CPX(:,:,:) = (0.d0,0.d0)

         CALL REAL_FOURIER_TRANS_TH (.true.)
         CALL REAL_FOURIER_TRANS_Rth (.true.)
         call allocation_Fth(.false.)

        
        DO N=1,N_TH
          CTHX(:,:,:,N) = (0.d0,0.d0)
        ENDDO


        if(rank.EQ.0) write(*,*) 'Reading flow...'
        CALL READ_FLOW
        if(rank.EQ.0) write(*,*) 'Done reading flow'

        IF( N_TH .gt. 0) then
         DO J=0,NY+1
          THBAR(J,1) = theta_0
         ENDDO
         DO J=0,NY+1
          THBAR(J,N_TH+1) = rho_0/1000.0d0 ;
         ENDDO
        ENDIF

        IF ( INT_TREAT ) THEN
!!  Need to add later

        ENDIF

!C Initialize flow.
      IF (RESET_TIME .OR. CREATE_NEW_FLOW) THEN
        PREVIOUS_TIME_STEP=0
        TIME_STEP=0
        TIME=0
      END IF

        CALL SAVE_STATS(.FALSE.)
        INIT_FLAG = .TRUE.       
        
!         IF (NUM_PER_DIR.EQ.3) THEN
! !C          CALL POISSON_P_PER
!         ELSEIF (NUM_PER_DIR.EQ.2) THEN
!           CALL POISSON_P_CHAN
!         ELSEIF (NUM_PER_DIR.EQ.1) THEN
!           CALL POISSON_P_DUCT
!         ELSEIF (NUM_PER_DIR.EQ.0) THEN
!           CALL POISSON_P_CAV
!         END IF

        CALL SAVE_STATS(.FALSE.) 

      END IF


      RETURN
      END





!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS(FINAL)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
 !     INCLUDE 'header_duct'
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mpi_var, only: rank      
implicit none


     LOGICAL FINAL

      IF (NUM_PER_DIR.EQ.3) THEN
!C        CALL SAVE_STATS_PER(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
!C        CALL SAVE_STATS_CHAN(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        CALL SAVE_STATS_DUCT(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        CALL SAVE_STATS_CAV(FINAL)          
      END IF

      if (rank ==0) then
      write(*,*) 'done save_stats diablo'
      endif

!       IF (FINAL) THEN
!         IF (NUM_PER_DIR.EQ.3) THEN
! !C          CALL VIS_FLOW_PER         
!         ELSEIF (NUM_PER_DIR.EQ.2) THEN
! !          CALL VIS_FLOW_CHAN         
!         ELSEIF (NUM_PER_DIR.EQ.1) THEN
!           CALL VIS_FLOW_DUCT          
!         ELSEIF (NUM_PER_DIR.EQ.0) THEN
!           CALL VIS_FLOW_CAV         
!         END IF
!       END IF

      RETURN
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE READ_FLOW
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
 !     INCLUDE 'header_duct'
  
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mpi_var, only: rank      
implicit none


      CHARACTER*35 FNAME
      CHARACTER*35 FNAME_TH(N_TH)
      INTEGER I, J, K, N, NUM_PER_DIR_T


     IF (USE_MPI) THEN
       I = RANK+1

        FNAME ='last_saved/diablo_'                 &
              //CHAR(MOD(I,1000)/100+48)            &
              //CHAR(MOD(I,100)/10+48)              &
              //CHAR(MOD(I,10)+48) // '.start'
        DO N=1,N_TH
           FNAME_TH(N)='last_saved/diablo_th'       &
             //CHAR(MOD(N,100)/10+48)               &
             //CHAR(MOD(N,10)+48) //'_'             &
             //CHAR(MOD(I,1000)/100+48)             &
             //CHAR(MOD(I,100)/10+48)               &
             //CHAR(MOD(I,10)+48) // '.start'
        ENDDO
      ELSE
      FNAME='diablo.start'
      DO N=1,N_TH
        FNAME_TH(N)='diablo_th'          &
             //CHAR(MOD(N,100)/10+48)    &
             //CHAR(MOD(N,10)+48) // '.start'
      END DO
      ENDIF

 
      if(rank==0)then
      WRITE(6,*)   'Reading flow from ',FNAME, NX2P
      WRITE(6,*)   'Reading flow from theta data ', FNAME_TH(1)
      endif

      OPEN(UNIT=10,FILE=FNAME,STATUS="OLD",FORM="UNFORMATTED")
      READ (10) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP


      write(*,*) 'NX_T, NY_T, NZ_T: ',NX_T,NY_T,NZ_T,NUM_PER_DIR_T, &
                  NUM_PER_DIR

      IF ((NX .NE. NX_T) .OR. (NY .NE. NY_T) .OR. (NZ .NE. NZ_T))THEN
           WRITE(*,*) NX_T, NY_T, NZ_T
           STOP 'Error: old flowfield wrong dimensions. '
      ENDIF
         
      IF (NUM_PER_DIR .NE. NUM_PER_DIR_T) &
          STOP 'Error: old flowfield wrong NUM_PER_DIR. '

      write(6,*) 'READING FLOW for rank =', rank

      IF (NUM_PER_DIR.EQ.3) THEN

      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        
      ELSEIF (NUM_PER_DIR.EQ.1) THEN

        READ (10) (((CU1X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1), &
                 (((CU2X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1),  &
                 (((CU3X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1),  &
                 (((CPX(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1) 
         Write(6,*) 'Done with velocity field'
        DO N=1,N_TH
! Specify in input.dat which scalars are to be read
          OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="OLD" &
                ,FORM="UNFORMATTED")
          READ (11) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP
          READ (11) (((CTHX(I,K,J,N) &
                ,I=0,NX2P),K=0,NZ+1),J=0,NY+1)
         CLOSE(11)
         Write(6,*) 'Done with theta field', N,' for rank = ', rank
        END DO

      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        

      END IF
      CLOSE(10)
      CLOSE(11)

!C Apply initial boundary conditions, set ghost cells
!       IF (USE_MPI) THEN
!         call APPLY_BC_VEL_MPI
!       ELSE
        call APPLY_BC_VEL_LOWER
        call APPLY_BC_VEL_UPPER
        call APPLY_BC_VEL_LEFT
        call APPLY_BC_VEL_RIGHT
!      END IF

      RETURN
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_FLOW(FINAL)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!      INCLUDE 'header_duct'
 
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mpi_var, only : rank      
implicit none


      CHARACTER*35 FNAME, frame(4),frame_th(4)
      CHARACTER*35 FNAME_TH(N_TH)
      INTEGER      I, J, K, N, m, NY_min, NY_max,no_p,ny_p
      LOGICAL      FINAL
      PARAMETER  (no_p=4, ny_p=66)
      
      IF (USE_MPI) THEN
       I=1+RANK
       IF (FINAL) THEN
         FNAME='last_saved/diablo_'      &
             //CHAR(MOD(I,1000)/100+48)  &
             //CHAR(MOD(I,100)/10+48)    &
             //CHAR(MOD(I,10)+48) // '.res'


         DO N=1,N_TH
            FNAME_TH(N)='last_saved/diablo_th'  &
              //CHAR(MOD(N,100)/10+48)          &
              //CHAR(MOD(N,10)+48) //'_'        &
              //CHAR(MOD(I,1000)/100+48)        &
              //CHAR(MOD(I,100)/10+48)          &
              //CHAR(MOD(I,10)+48) // '.res'
         END DO

       ELSE
         FNAME='last_saved/diablo_'        &
              //CHAR(MOD(I,1000)/100+48)   &
              //CHAR(MOD(I,100)/10+48)     &
              //CHAR(MOD(I,10)+48) // '.saved'

         DO N=1,N_TH
            FNAME_TH(N)='last_saved/diablo_th'  &
              //CHAR(MOD(N,100)/10+48)          &
              //CHAR(MOD(N,10)+48) //'_'        &
              //CHAR(MOD(I,1000)/100+48)        &
              //CHAR(MOD(I,100)/10+48)          &
              //CHAR(MOD(I,10)+48) // '.saved'
         END DO
       ENDIF
      ELSE
       IF (FINAL) THEN
         FNAME='diablo.res'
         DO N=1,N_TH
           FNAME_TH(N)='diablo_th'           &
             //CHAR(MOD(N,100)/10+48)        &
             //CHAR(MOD(N,10)+48) // '.res'
         END DO
       ELSE
         FNAME='diablo.saved'
         DO N=1,N_TH
           FNAME_TH(N)='diablo_th'          &
             //CHAR(MOD(N,100)/10+48)       &
             //CHAR(MOD(N,10)+48) // '.saved'
         END DO

       END IF
      ENDIF

      if (rank .eq. 0) then  
       WRITE(6,*) 'Writing flow to ',FNAME, 'Rank = ', RANK
      endif

      OPEN(UNIT=10,FILE=FNAME,STATUS="UNKNOWN",FORM="UNFORMATTED")
      WRITE(10) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP


      IF (NUM_PER_DIR.EQ.3) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      ELSEIF (NUM_PER_DIR.EQ.2) THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        WRITE(10) (((CU1X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1), &
                 (((CU2X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1),  &
                 (((CU3X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1), &
                 (((CPX(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1)
        DO N=1,N_TH
          OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="UNKNOWN" &
            ,FORM="UNFORMATTED")
         WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
         WRITE(11) (((CTHX(I,K,J,N),I=0,NX2P),K=0,NZ+1),J=0,NY+1)
         CLOSE(11)
        END DO
 
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      END IF
      CLOSE(10)
      CLOSE(11)

      RETURN
      END


!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FILTER(n)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
     ! INCLUDE 'header_duct'

use TIME_STEP_VAR

      integer n

      IF (NUM_PER_DIR.EQ.3) THEN
!C        CALL FILTER_PER
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
!C        CALL FILTER_CHAN(n)
      END IF

      RETURN
      END

