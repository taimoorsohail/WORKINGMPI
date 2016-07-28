
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_TH_DUCT
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C Initialize the scalar fields
!C In this subroutine, you should initialize each scalar field for the
!C particular problem of interest

!      INCLUDE 'header_duct'

use ntypes
use Domain
use Grid
use Fft_var, only : pi
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use mpi_var, only : rank
      
implicit none

      INTEGER I,J,K,N,MM,NN,M
      REAL*8 RNUM1,RNUM2,RNUM3,DAMP_FACT,SUM1, &
            FACT1,FACT2,FACT3,FACT4,FACT5,FACT
      REAL*8 NGY(NY), NGZ(NZ)

      DO N=1,N_TH
       IF (CREATE_NEW_TH(N)) THEN
        if (rank .eq. 0) then    
         write(6,*) 'A new thfield has been created for ', 'rank = ', rank
        end if     
        IF ( IC_TYPE.eq.0 ) THEN
        ELSE IF ( (IC_TYPE.eq.5).OR.(IC_TYPE.eq.6) )then
         DO J=0,NY+1
          DO K=0,NZ+1
           DO I=0,NXP
             THX(I,K,J,N)=U3_BAR(1,J) 
           END DO
          END DO
         ENDDO
        ELSE IF ( IC_TYPE .eq. 7 )then
         DO J=0,NY+1
          DO K=0,NZ+1
           DO I=0,NXP
             THX(I,K,J,N)=0.0d0 
           END DO
          END DO
          IF (N .eq. 1) THEN
!        setting bacground for temperature
           THBAR(J,N) = theta_0
          ELSEIF (N .eq. 2) THEN
!        setting bacground for Salinity
          ENDIF
         ENDDO
        ENDIF
        ENDIF
       ENDDO

        IF (N_TH .gt. 0) THEN
!        setting bacground for density
         DO J=0,NY+1 
          THBAR(J,N_TH+1) = rho_0/1000.0d0 ;
         ENDDO 
        ENDIF

        IF (CREATE_NEW_TH(1)) THEN
         CALL REAL_FOURIER_TRANS_TH (.true.)
         CALL REAL_FOURIER_TRANS_Rth (.true.)
         call allocation_Fth(.false.)
        ENDIF     

      RETURN
      END 







!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_DUCT
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!      INCLUDE 'header_duct'

use ntypes
use Domain
use Grid
use Fft_var, only : pi
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use mpi_var, only : rank      
implicit none

      INTEGER I,J,K,MM,NN,M,N
      REAL*8 RNUM1,RNUM2,RNUM3,DAMP_FACT,SUM1,  &
            FACT1,FACT2,FACT3,FACT4,FACT5,FACT
      REAL*8 NGY(NY), NGZ(NZ)
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

!C Initialize the random number generator
      CALL RANDOM_SEED(SIZE = K)
      Allocate (seed(1:K))
      seed(1:K)=10
      CALL RANDOM_SEED(PUT = seed)

!C UBULK0 and KICK should be set in input.dat

!C Set the laminar velocity profile in physical space
      if (rank .eq. 0) then  
       write(6,*) 'A new Velocity Field has been created for ', 'rank = ', rank
      end if

       IF (IC_TYPE.eq.0) then
!C For closed duct flow
      mm = 5 
      nn = 5 
      fact  = LZ/LY

      DO K =1,NZ
           NGZ(K) = 2.0*GZF(K)/LZ ;
      ENDDO 
      DO J =1,NY
           NGY(J) = 2.0*GYF(J)/LY ;
      ENDDO

       DO J=1,NY
         DO K=1,NZ
          SUM1 = 0.0 ;

         DO m=0,mm
           DO n=0,nn

           fact1 = (2.0*m+1) ;
           fact2 = (2.0*n+1) ;
           fact3 = fact1*fact2**3.0 + fact**2.0*fact2*fact1**3.0 ;
           fact5 = ((-1.0)**(m+n))*fact/fact3 ;
           sum1 = sum1 + fact5*cos(fact1*pi*NGZ(K)/2)*  &
                  cos(fact2*pi*NGY(J)/2) ;
          ENDDO
         ENDDO
        
           DO I=0,NXP
             U1X(I,K,J)=UBULK0*sum1 
             U2X(I,K,J)=0.
             U3X(I,K,J)=0.
           END DO
         END DO
      END DO
      else if ((IC_TYPE.eq.1) .OR. (IC_TYPE.eq.4) ) then 
!C For open channel flow :
       DO K=0,NZM
         DO I=0,NXP
           DO J=1,NY
!            U1(I,K,J)=-(3./2.)*UBULK0*GYF(J)**2.+3.*UBULK0*GYF(J)
!            U1(I,K,J)=(-GYF(J)**2.d0+(NU+LY**2.d0)*GYF(J)/LY)/NU
             U1X(I,K,J)= (-GYF(J)**2.d0+ 2.0*LY*GYF(J))/LY**2.0
!             U1(I,K,J)=0.d0
             U2X(I,K,J)=0.
             U3X(I,K,J)=0.
           END DO
           U1X(I,K,0)=0.
           U3X(I,K,0)=0.
           U1X(I,K,NY+1)=0.
           U3X(I,K,NY+1)=0.
         END DO
      END DO
      else if (IC_TYPE.eq.2) then
!C For Couette flow:
       DO J=0,NY
         DO K=0,NZM
           DO I=0,NXP
             U1X(I,K,J)=gyf(j)
             U2X(I,K,J)=0.
             U3X(I,K,J)=0.
           END DO
         END DO
      END DO
      else if (IC_TYPE.eq.3) then
! Shear layer
       DO J=0,NY+1
         DO K=0,NZ+1
           DO I=0,NX+1
             U1X(I,K,J)=TANH(GYF(J)*20.d0)
             U2X(I,K,J)=0.d0
             U3X(I,K,J)=0.d0
            END DO
          END DO
        END DO
       else if ( (IC_TYPE.eq.5).OR.(IC_TYPE.eq.6) )then        
         DO J=0,NY+1
          DO K=0,NZ+1
           DO I=0,NXP
             U1X(I,K,J)=U1_BAR(J)
             U2X(I,K,J)=0.0!U2_BAR(K,J)
             U3X(I,K,J)=U3_BAR(1,J) 
           END DO
         END DO
        ENDDO
      else if ( (IC_TYPE.eq.7) )then        
         DO J=0,NY+1
          DO K=0,NZ+1
           DO I=0,NXP
             U1X(I,K,J)=0.0d0
             U2X(I,K,J)=0.0d0
             U3X(I,K,J)=0.0d0
           END DO
         END DO
        ENDDO
      end if

!C Zero the ghost cells
       IF (.NOT.USE_MPI) THEN       
       DO K=0,NZM
         DO I=0,NXP
           U1X(I,K,0)=0.
           U2X(I,K,0)=0.
           U3X(I,K,0)=0.
           U1X(I,K,NY+1)=0.
           U2X(I,K,NY+1)=0.
           U3X(I,K,NY+1)=0.
         END DO
      END DO
      END IF
   
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
!c        WRITE(556,*) CU1(0,NY/2,NZ/2)
    
      if (rank .eq. 0) then
      WRITE(*,*) 'KICK is : ',KICK
      endif 

      if (rank .eq. 0) then 
      DO I=1,min(NX2P,NX2P_L)
        DO J=1,NY
          DO K=1,NZ
!C Now, give the velocity field a random perturbation
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            CALL RANDOM_NUMBER(RNUM3)

            IF (IC_TYPE.eq.3) THEN
!C If we are initializing with a shear layer 
              CU1X(I,K,J)=CU1X(I,K,J)  &
                  +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
              CU2X(I,K,J)=CU2X(I,K,J)  &
                  +(RNUM1-0.5)*KICK*EXP(-(GY(J)*20.d0)**2.d0)
              CU3X(I,K,J)=CU3X(I,K,J)  &
                  +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
            ELSEIF (IC_TYPE.eq.7) THEN
!C If we are initializing with bottom boundary layer
              CU1X(I,K,J)=CU1X(I,K,J)  &
                  +(RNUM1-0.5)*KICK*EXP(-((GYF(J)-0.2)*50.d0)**2.d0)
              CU2X(I,K,J)=CU2X(I,K,J)  &
                  +(RNUM1-0.5)*KICK*EXP(-((GY(J)-0.2)*50.d0)**2.d0)
              CU3X(I,K,J)=CU3X(I,K,J)  &
                  +(RNUM1-0.5)*KICK*EXP(-((GYF(J)-0.2)*50.d0)**2.d0)
             ELSE
              CU1X(I,K,J)=CU1X(I,K,J)+(RNUM1-0.5)*KICK
              CU2X(I,K,J)=CU2X(I,K,J)+(RNUM2-0.5)*KICK
              CU3X(I,K,J)=CU3X(I,K,J)+(RNUM3-0.5)*KICK
            END IF
          END DO

        END DO
      END DO      
      else
      DO I=0,min(NX2P,NX2P_L)
        DO J=1,NY
          DO K=1,NZ
!C Now, give the velocity field a random perturbation
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            CALL RANDOM_NUMBER(RNUM3)

            IF (IC_TYPE.eq.3) THEN
!C If we are initializing with a shear layer 
              CU1X(I,K,J)=CU1X(I,K,J)  &
                  +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
              CU2X(I,K,J)=CU2X(I,K,J)  &
                  +(RNUM1-0.5)*KICK*EXP(-(GY(J)*20.d0)**2.d0)
              CU3X(I,K,J)=CU3X(I,K,J)  &
                  +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
            ELSEIF (IC_TYPE.eq.7) THEN
!C If we are initializing with bottom boundary layer
              CU1X(I,K,J)=CU1X(I,K,J)  &
                  +(RNUM1-0.5)*KICK*EXP(-((GYF(J)-0.2)*50.d0)**2.d0)
              CU2X(I,K,J)=CU2X(I,K,J)  &
                  +(RNUM1-0.5)*KICK*EXP(-((GY(J)-0.2)*50.d0)**2.d0)
              CU3X(I,K,J)=CU3X(I,K,J)  &
                  +(RNUM1-0.5)*KICK*EXP(-((GYF(J)-0.2)*50.d0)**2.d0)
             ELSE
              CU1X(I,K,J)=CU1X(I,K,J)+(RNUM1-0.5)*KICK
              CU2X(I,K,J)=CU2X(I,K,J)+(RNUM2-0.5)*KICK
              CU3X(I,K,J)=CU3X(I,K,J)+(RNUM3-0.5)*KICK
            END IF
          END DO

        END DO
      END DO
      endif  

     

!       IF (USE_MPI) THEN
!         CALL GHOST_CHAN_MPI
!       END IF

!c      CALL SAVE_STATS_DUCT(.FALSE.)


!C Apply Boundary conditions to velocity field
!       IF (USE_MPI) THEN
!         CALL APPLY_BC_VEL_MPI
!       ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
        CALL APPLY_BC_VEL_LEFT
        CALL APPLY_BC_VEL_RIGHT
!      END IF

!C Remove the divergence of the velocity field

      CALL SAVE_STATS_DUCT(.FALSE.)

       write(6,*) 'Done with flow field creation for rank = ', rank
!C Apply Boundary conditions to velocity field

!       IF (USE_MPI) THEN
!         CALL APPLY_BC_VEL_MPI
!       ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
        CALL APPLY_BC_VEL_LEFT
        CALL APPLY_BC_VEL_RIGHT
!      END IF


       INIT_FLAG = .TRUE.
       CALL REM_DIV_DUCT

!       IF (USE_MPI) THEN
!         CALL GHOST_CHAN_MPI
!       END IF

!C Get the pressure from the poisson equation
!      CALL POISSON_P_DUCT

!       IF (USE_MPI) THEN
!         CALL GHOST_CHAN_MPI
!       END IF

        
      CALL SAVE_STATS_DUCT(.FALSE.)
      
      RETURN
      END
