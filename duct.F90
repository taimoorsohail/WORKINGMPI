!C******************************************************************************|
!C duct.f, the duct-flow solvers for diablo.                        VERSION 0.9
!C These solvers were written by ? and ? (spring 2001).
!C******************************************************************************|

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_DUCT
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C Initialize any constants here
!      INCLUDE 'header_duct'
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG, bc
use mpi_var

implicit none

      INTEGER J, K,N
      
      PI=4.D0*ATAN(1.D0)

! At YMIN Location
!c        write(*,*) 'U_BC_YMIN: ',U_BC_YMIN
        IF (W_BC_YMIN.EQ.0) THEN
          JSTART=2
        ELSE IF (W_BC_YMIN.EQ.1) THEN
          JSTART=1
        ELSE IF (W_BC_YMIN.EQ.7) THEN
          JSTART=1
        ELSE
          JSTART=2
        END IF
! At ZMIN Location  
!c        write(*,*) 'U_BC_ZMIN: ',U_BC_ZMIN
        IF (U_BC_ZMIN.EQ.0) THEN
          ZSTART=2
        ELSE IF (U_BC_ZMIN.EQ.1) THEN
          ZSTART=1
        ELSE
          ZSTART=2
        END IF
! Now, set the indexing for the scalar equations
! At YMIN Location
        DO N=1,N_TH
          IF (TH_BC_YMIN(N).EQ.0) THEN
            JSTART_TH(N)=2
          ELSE IF (TH_BC_YMIN(N).EQ.1) THEN
            JSTART_TH(N)=1
          ELSE IF (TH_BC_YMIN(N).EQ.7) THEN
            JSTART_TH(N)=1
          ELSE
            JSTART_TH(N)=2
          END IF
        END DO
! At ZMIN Location    
        DO N=1,N_TH
          IF (TH_BC_ZMIN(N).EQ.0) THEN
            ZSTART_TH(N)=2
          ELSE IF (TH_BC_ZMIN(N).EQ.1) THEN
            ZSTART_TH(N)=1
          ELSE
            ZSTART_TH(N)=2
          END IF
        END DO      

! At YMAX Location
!c        write(*,*) 'U_BC_YMAX: ',U_BC_YMAX
        IF (U_BC_YMAX.EQ.0) THEN
          JEND=NY-1
        ELSE IF (U_BC_YMAX.EQ.1) THEN
          JEND=NY
        ELSE
          JEND=NY-1
        END IF

! At ZMAX Location
!c        write(*,*) 'U_BC_ZMAX: ',U_BC_ZMAX
        IF (U_BC_ZMAX.EQ.0) THEN
          ZEND=NZ-1
        ELSE IF (U_BC_ZMAX.EQ.1) THEN
          ZEND=NZ
        ELSE
          ZEND=NZ-1
        END IF

! Set the upper and lower limits of timestepping of the scalar equations
! At YMAX Location 
        DO N=1,N_TH
        IF(TH_BC_YMAX(N).EQ.0)THEN
          JEND_TH(N)=NY-1
        ELSE IF (TH_BC_YMAX(N).EQ.1) THEN
          JEND_TH(N)=NY
        ELSE IF (TH_BC_YMAX(N).EQ.5) THEN
          JEND_TH(N)=NY
        ELSE IF (TH_BC_YMAX(N).EQ.7) THEN
          JEND_TH(N)=NY
        ELSE
          JEND_TH(N)=NY-1
        END IF
        END DO
! At ZMAX Location
        DO N=1,N_TH
        IF (TH_BC_ZMAX(N).EQ.0) THEN
          ZEND_TH(N)=NZ-1
        ELSE IF (TH_BC_ZMAX(N).EQ.1) THEN
          ZEND_TH(N)=NZ
        ELSE IF (TH_BC_ZMAX(N).EQ.5) THEN
          ZEND_TH(N)=NZ
        ELSE
          ZEND_TH(N)=NZ-1
        END IF
        END DO

       IF (rank .eq. 0 ) then
       WRITE(6,*) '###################################################'
       WRITE(6,*)'Boundary condition for U'
       WRITE(6,*)'U_BC_YMIN',U_BC_YMIN,'U_BC_YMAX',U_BC_YMAX
       WRITE(6,*)'U_BC_ZMIN',U_BC_ZMIN,'U_BC_ZMAX',U_BC_ZMAX
       WRITE(6,*) '###################################################'
       WRITE(6,*)'U_BC_YMIN_C1',U_BC_YMIN_C1,'U_BC_YMAX_C1',U_BC_YMAX_C1
       WRITE(6,*)'U_BC_ZMIN_C1',U_BC_ZMIN_C1,'U_BC_ZMAX_C1',U_BC_ZMAX_C1

       WRITE(6,*)'Boundary condition for W'
       WRITE(6,*)'W_BC_YMIN',W_BC_YMIN,'W_BC_YMAX', W_BC_YMAX
       WRITE(6,*)'W_BC_ZMIN',W_BC_ZMIN,'W_BC_ZMAX', W_BC_ZMAX
       WRITE(6,*) '###################################################'
       WRITE(6,*)'W_BC_YMIN_C1',W_BC_YMIN_C1,'W_BC_YMAX_C1',W_BC_YMAX_C1
       WRITE(6,*)'W_BC_ZMIN_C1',W_BC_ZMIN_C1,'W_BC_ZMAX_C1',W_BC_ZMAX_C1


       WRITE(6,*)'Boundary condition for V'
       WRITE(6,*)'V_BC_YMIN',V_BC_YMIN,'V_BC_YMAX',V_BC_YMAX
       WRITE(6,*)'V_BC_ZMIN',V_BC_ZMIN,'V_BC_ZMAX',V_BC_ZMAX
       WRITE(6,*) '###################################################' 
       WRITE(6,*)'V_BC_YMIN_C1',V_BC_YMIN_C1,'V_BC_YMAX_C1',V_BC_YMAX_C1
       WRITE(6,*)'V_BC_ZMIN_C1',V_BC_ZMIN_C1,'V_BC_ZMAX_C1',V_BC_ZMAX_C1


       WRITE(6,*)'#####################################################'
       WRITE(6,*) 'JSATRT', JSTART, 'JEND', JEND
       WRITE(6,*) 'ZSATRT', ZSTART, 'ZEND', ZEND 
       WRITE(6,*)'####################################################'
       DO N=1,N_TH
       WRITE(6,*)'Boundary condition for TH'
       WRITE(6,*)'TH_BC_YMIN',TH_BC_YMIN(N),'TH_BC_YMAX', TH_BC_YMAX(N)
       WRITE(6,*)'TH_BC_ZMIN',TH_BC_ZMIN(N),'TH_BC_ZMAX', TH_BC_ZMAX(N)
       WRITE(6,*) '###################################################'
      WRITE(6,*)'TH_BC_YMIN_C1',TH_BC_YMIN_C1(N),                       &
                'TH_BC_YMAX_C1',TH_BC_YMAX_C1(N)
      WRITE(6,*)'TH_BC_ZMIN_C1',TH_BC_ZMIN_C1(N),                       &
                 'TH_BC_ZMAX_C1',TH_BC_ZMAX_C1(N)
      WRITE(6,*)'#####################################################'
       WRITE(6,*) 'JSATRT_TH', JSTART_TH(N), 'JEND_TH', JEND_TH(N)
       WRITE(6,*) 'ZSATRT_TH', ZSTART_TH(N), 'ZEND_TH', ZEND_TH(N) 
       WRITE(6,*)'#####################################################' 
          
       ENDDO
       ENDIF


      RETURN
      END




!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_DUCT_1
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C Main time-stepping algorithm for the duct-flow case.
!C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
!C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_DUCT_2
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C Alternative time-stepping algorithm for the duct-flow case.
!C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
!C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END


!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_DUCT_3
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C Alternative time-stepping algorithm for the duct-flow case with ADI.
!C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
!C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
use omp_lib      
use les_chan_var
use mpi_var
implicit none


!      INCLUDE 'header_duct'

      INTEGER I,J,K,N      
      REAL(r8) TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, UBULK, &
            MEAN_U1,D
!C Communicate the information between ghost cells


!c      CALL GHOST_CHAN_MPI

!C Define the constants that are used in the time-stepping
!C For reference, see Numerical Renaissance
      TEMP1=NU * H_BAR(RK_STEP) / 2.0
      TEMP2=H_BAR(RK_STEP) / 2.0
      TEMP3=ZETA_BAR(RK_STEP) * H_BAR(RK_STEP)
      TEMP4=H_BAR(RK_STEP)
      TEMP5=BETA_BAR(RK_STEP) * H_BAR(RK_STEP)

!C First, we will compute the explicit RHS terms and store in Ri
!C Note, Momentum equation and hence the RHS is evaluated at the
!C corresponding velocity points.
     CR1X(:,:,:) = (0.d0,0.d0)
     CR2X(:,:,:) = (0.d0,0.d0)
     CR3X(:,:,:) = (0.d0,0.d0) 


      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CR1X(I,K,J)=CU1X(I,K,J)
          END DO
        END DO
      END DO

      DO J=2,NY 
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CR2X(I,K,J)=CU2X(I,K,J)
          END DO
        END DO
      END DO

      DO J=JSTART,JEND
        DO K=2,NZ
          DO I=0,NX2P
            CR3X(I,K,J)=CU3X(I,K,J)
          END DO
        END DO
      END DO


!C Add the R-K term from the rk-1 step

      IF (RK_STEP .GT. 1) THEN
        DO J=JSTART,JEND
          DO K=ZSTART,ZEND
            DO I=0,NX2P
              CR1X(I,K,J)=CR1X(I,K,J)+TEMP3*CF1X(I,K,J)
            END DO
          END DO
        END DO
        DO J=2,NY
          DO K=ZSTART,ZEND
            DO I=0,NX2P
              CR2X(I,K,J)=CR2X(I,K,J)+TEMP3*CF2X(I,K,J)
            END DO
          END DO
        END DO
        DO J=JSTART,JEND
          DO K=2,NZ
            DO I=0,NX2P
              CR3X(I,K,J)=CR3X(I,K,J)+TEMP3*CF3X(I,K,J)
            END DO
          END DO
        END DO
      END IF
          
!C Take the y-derivative of the pressure at GY points in Fourier space
      DO J=2,NY
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CS1X(I,K,J)=(CPX(I,K,J) - CPX(I,K,J-1)) / DY(J)
          END DO
        END DO
      END DO
    
!C Take the z-derivative of the pressure at GZ points in Fourier space
      DO J=JSTART,JEND
        DO K=2,NZ
          DO I=0,NX2P
             CS2X(I,K,J)=(CPX(I,K,J) - CPX(I,K-1,J)) / DZ(K)
          END DO
         END DO
       END DO

!C Add the pressure gradient to the RHS as explicit Euler
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CR1X(I,K,J)=CR1X(I,K,J)-TEMP4*(CIKXP(I)*CPX(I,K,J))
          END DO
        END DO
      END DO

      DO J=2,NY
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CR2X(I,K,J)=CR2X(I,K,J)-TEMP4*CS1X(I,K,J)
          END DO
        END DO
      END DO

      DO J=JSTART,JEND
        DO K=2,NZ
          DO I=0,NX2P
            CR3X(I,K,J)=CR3X(I,K,J)-TEMP4*CS2X(I,K,J)
          END DO
        END DO
      END DO 

!C Here, add the constant, forcing pressure gradient
!C There are several ways of doing this
      IF (F_TYPE.EQ.1) THEN 
!C Add forcing for a constant pressure gradient
        DO J=JSTART,JEND
         DO K=ZSTART,ZEND 
          CR1X(0,K,J)=CR1X(0,K,J)-TEMP4*PX0
         END DO
        END DO
      ELSE IF (F_TYPE.EQ.0) THEN  
!C Add the mean pressure gradient to keep Ubulk constant
!C This section needs to be parallelized 

      ELSE IF (F_TYPE.EQ.2) THEN
!C If oscillatory pressure gradient
        DO J=JSTART,JEND
          DO K=ZSTART,ZEND
            CR1X(0,K,J)=CR1X(0,K,J)-TEMP4*(PX0 +  &
                        AMP_OMEGA0*cos(OMEGA0*TIME))
          END DO
        END DO

      ELSE IF (F_TYPE.EQ.4) THEN
!C If oscillatory pressure gradient
        DO J=JSTART,JEND
          DO K=ZSTART,ZEND
           CR1X(0,K,J)=CR1X(0,K,J)-TEMP4*(PX0  &
                    +AMP_OMEGA0*cos(OMEGA0*TIME))* cos(ANG_BETA)
          END DO
        END DO
      
      ELSE IF (F_TYPE .EQ. 5) THEN
              
!C End if forcing type
      END IF

!C Now compute the term R-K term Ai
!C Compile terms of Ai in CFi which will be saved for next time step
!C First, store the horizontal viscous terms in CFi
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF1X(I,K,J)=-NU * KX2P(I) * CU1X(I,K,J) 
          END DO
        END DO
      END DO

      DO J=2,NY
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF2X(I,K,J)=-NU * KX2P(I) * CU2X(I,K,J) 
          END DO 
        END DO
      END DO

      DO J=JSTART,JEND
        DO K=2,NZ
          DO I=0,NX2P
            CF3X(I,K,J)=-NU * KX2P(I) * CU3X(I,K,J)
          END DO
        END DO
      END DO
      

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF1X(I,K,J)=CF1X(I,K,J) + (f_0 + beta_f*GZF(K))*0.50d0*(CU3X(I,K+1,J)+CU3X(I,K,J))
          END DO
        END DO
      END DO

      DO J=JSTART,JEND
        DO K=2,NZ
          DO I=0,NX2P
            CF3X(I,K,J)=CF3X(I,K,J) - (f_0 + beta_f*GZ(K))*0.50d0*(CU1X(I,K,J)+CU1X(I,K-1,J))
          END DO
        END DO
      END DO


! Do for each scalar
      DO N=1,N_TH
! If a scalar contributes to the denisty, RI_TAU is not equal to zero and
! add the buoyancy term as explicit R-K.  Don't add the 0,0 mode which 
! corresponds to a plane average.  The plane averaged density balances
! the hydrostratic pressure component.
      IF ((F_TYPE .EQ. 4).OR.(F_TYPE.EQ.5)) THEN        
        DO J=2,NY 
         DO K=ZSTART,ZEND
           DO I=0,NX2P 
! Use second order interpolation
              CF2X(I,K,J)=CF2X(I,K,J) - RI_TAU(N)*cos(ANG_BETA)* &
            (CTHX(I,K,J,N)*DYF(J-1)+CTHX(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
           END DO
         END DO
        END DO        
       
        DO J=JSTART,JEND
         DO K=ZSTART,ZEND
          DO I=0,NX2P
             CF1X(I,K,J)=CF1X(I,K,J)-RI_TAU(N)*sin(ANG_BETA)*CTHX(I,K,J,N)
          END DO
         END DO
        END DO

      ELSE

      IF (DEV_BACK_TH) THEN
       DO J=2,NY
         DO K=ZSTART,ZEND
           DO I=0,NX2P
! Use second order interpolation
             CF2X(I,K,J)=CF2X(I,K,J)+RI_TAU(N)* &
           (CTHX(I,K,J,N)*DYF(J-1)+CTHX(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
          END DO
        END DO
       END DO
       ELSE
        DO J=2,NY
         DO K=ZSTART,ZEND
           DO I=0,NX2P
! Use second order interpolation
             CF2X(I,K,J)=CF2X(I,K,J)+RI_TAU(N)*  &
           (CTHX(I,K,J,N)*DYF(J-1)+CTHX(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
          END DO
        END DO
       END DO
      ENDIF

      END IF
       
    
! Now, compute the RHS vector for the scalar equations
! Since TH is defined at horizontal velocity points, the
! scalar update equation will be very similar to the horizontal
! velocity update.

! We will store the RHS scalar terms in CRTH, RTH
! The k-1 term for the R-K stepping is saved in FTH, CFTH
    
      CRTHX(:,:,:,N) = (0.d0,0.d0)
! First, build the RHS vector, use CRTH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
          CRTHX(I,K,J,N)=CTHX(I,K,J,N)
         ENDDO
       END DO
      END DO
! Add term from k-2 step to free up CFTH variable
      IF (RK_STEP .GT. 1) THEN
        DO J=JSTART_TH(N),JEND_TH(N)
          DO K=ZSTART_TH(N),ZEND_TH(N)
            DO I=0,NX2P
              CRTHX(I,K,J,N)=CRTHX(I,K,J,N)+TEMP3*CFTHX(I,K,J,N)
            END DO
          END DO
        END DO
       END IF

! Now compute the explicit R-K term Ai
! Compile terms of Ai in CFi which will be saved for next time step
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
            CFTHX(I,K,J,N)=-(NU/PR(N)) * KX2P(I) * CTHX(I,K,J,N)
          END DO
        END DO
      END DO
      
!C Need to be forced by  background temperature gradient    
      IF (DEV_BACK_TH) THEN
       DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
            CFTHX(I,K,J,N)=CFTHX(I,K,J,N) - 0.5*(CU2X(I,K,J+1)+CU2X(I,K,J))*  &
          (THBAR(j+1,n)-THBAR(j-1,n))/(2.*DYF(j))
          END DO
        END DO
      END DO
      END IF

! End of  loop for passive scalars (N_TH)   
      END DO

!      IF (IBM) THEN
!        CALL sink_momentum
!      ENDIF 


      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
!C If we have created new flow with random perturbations, wait for a
!C spinup before applying the subgrid model for stability purposes
!C In the process, Ui is converted to physical space
!          write(*,*) 'calling les'
          call les_chan
!C Add the subgrid scale scalar flux to the scalar equations
          IF (N_TH .gt. 0) then
           DO N=1,N_TH
 !           call les_chan_th(N)
           ENDDO
          ENDIF

!         write(*,*) 'Done with les'
      !C convert to physical space.
!      CALL REAL_FOURIER_TRANS_U1 (.false.)
!      CALL REAL_FOURIER_TRANS_U2 (.false.)
!      CALL REAL_FOURIER_TRANS_U3 (.false.)
! Transform THETA to physical space for computation of nonlinear terms
! Here pass the first location in memory of the array for scalar n
      IF (N_TH .gt. 0) then
       CALL REAL_FOURIER_TRANS_TH (.false.)
      ENDIF


      ELSE 
!C If the subgrid model hasn't been called, then it is necessary to 
!C convert to physical space.
      CALL REAL_FOURIER_TRANS_U1 (.false.)
      CALL REAL_FOURIER_TRANS_U2 (.false.)
      CALL REAL_FOURIER_TRANS_U3 (.false.)
! Transform THETA to physical space for computation of nonlinear terms
! Here pass the first location in memory of the array for scalar n
      IF (N_TH .gt. 0) then
       CALL REAL_FOURIER_TRANS_TH (.false.)
      ENDIF
      
      END IF

!c      IF ( F_TYPE .EQ. 1 ) THEN

      S1X(:,:,:) =0.d0

!C U1*U3
      DO J=JSTART,JEND
        DO K=2,NZ
          DO I=0,NXP
            S1X(I,K,J)=U3X(I,K,J)*((DZF(K)*U1X(I,K,J)            &
                     +DZF(K-1)*U1X(I,K-1,J))/(2.*DZ(K)))
!C           S1(I,K,J)=U1(I,K,J)*((DZ(K)*U3(I,K,J)
!C     &                +DZ(K+1)*U3(I,K+1,J))/(2.D0*DZF(K)))
          END DO
        END DO
      END DO
      
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

      
      DO J=JSTART,JEND
        DO K=2,NZ
          DO I=0,NX2P
            CF3X(I,K,J)=CF3X(I,K,J) - CIKXP(I) * CS1X(I,K,J) 
          END DO
        END DO
      END DO

!C U1*U1
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            S1X(I,K,J)=U1X(I,K,J)*U1X(I,K,J)
          END DO
        END DO
      END DO
      
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

      
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF1X(I,K,J)=CF1X(I,K,J) -  CIKXP(I) * CS1X(I,K,J) 
          END DO
        END DO
      END DO

!C U1*U2
      DO J=2,NY
        DO K=ZSTART,ZEND
          DO I=0,NXP
            S1X(I,K,J)=((DYF(J)*U1X(I,K,J)                  &
                     +DYF(J-1)*U1X(I,K,J-1))/(2.*DY(J)))   &
                     *U2X(I,K,J)
!c            S1(I,K,J)=((DY(J+1)*U2(I,K,J+1)
!c     &                +DY(J)*U2(I,K,J))/(2.D0*DYF(J)))
          END DO
        END DO
      END DO
      
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)
      
      DO J=2,NY
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF2X(I,K,J)=CF2X(I,K,J) - CIKXP(I) * CS1X(I,K,J) 
          END DO
        END DO
      END DO

! Add the vertical (y) derivative and spanwise (z) derivative term explicitly 

!C d(U1*U2)/dY  +  d(U1*U3)/dZ
                
       DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
           S1X(I,K,J) = (U2X(I,K,J+1)*((U1X(I,K,J+1)*DYF(J+1)+U1X(I,K,J)*             &  !C d(U1*U2)/dY
                DYF(J))/(2.d0*DY(J+1))) - U2X(I,K,J)*((U1X(I,K,J)*DYF(J)+             &
                U1X(I,K,J-1)*DYF(J-1))/(2.d0*DY(J))))/DYF(J)                          &
                +                                                                     & ! d(U1*U3)/dZ
                      (U3X(I,K+1,J)*((U1X(I,K+1,J)*DZF(K+1)+U1X(I,K,J)               &
                      *DZF(K))/(2.D0*DZ(K+1)))-U3X(I,K,J)*((U1X(I,K,J)                & 
                      *DZF(K)+U1X(I,K-1,J)*DZF(K-1))/(2.D0*DZ(K))))                   &
                      /DZF(K)
          END DO
        END DO
      END DO

      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF1X(I,K,J)=CF1X(I,K,J) - CS1X(I,K,J) 
         END DO
        END DO
      END DO

! !C d(U1*U3)/dZ
! 
!       DO J=JSTART,JEND
!         DO K=ZSTART,ZEND
!           DO I=0,NXP
! !c            S1(I,K,J)= ((U1(I,K+1,J)*U3(I,K+1,J) + U1(I,K,J)*U3(I,K+1,J)
! !c     &     - U1(I,K,J)*U3(I,K,J) - U1(I,K-1,J)*U3(I,K,J))/(2.d0*DZF(K)))
!            S1X(I,K,J) = (U3X(I,K+1,J)*((U1X(I,K+1,J)*DZF(K+1)+U1X(I,K,J)               &
!                       *DZF(K))/(2.D0*DZ(K+1)))-U3X(I,K,J)*((U1X(I,K,J)               & 
!                       *DZF(K)+U1X(I,K-1,J)*DZF(K-1))/(2.D0*DZ(K))))                 &
!                       /DZF(K)             
!           END DO
!         END DO
!       END DO
! 
!       CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
!       varp(:,:,:) = 0.d0
!       DO I=0,NXM
!        varp(I,:,:)=S1Z(I,:,:)
!       ENDDO
!       CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
!       DO I=0,NKX
!       CS1Z(I,:,:)=cvarp(I,:,:)
!       ENDDO
!       CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
!       CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)
! 
!       DO J=JSTART,JEND
!         DO K=ZSTART,ZEND
!           DO I=0,NX2P
!             CF1X(I,K,J)=CF1X(I,K,J) - CS1X(I,K,J)
!          END DO
!         END DO
!       END DO

! d(U2*U2)/dY + d(U2*U3)/dZ

      DO J=2,NY
        DO K=ZSTART,ZEND
          DO I=0,NXP
            S1X(I,K,J)=                                                           &   ! d(U2*U2)/dY
         (0.25d0*(U2X(I,K,J)+U2X(I,K,J+1))**2.d0                                  & 
         -0.25d0*(U2X(I,K,J)+U2X(I,K,J-1))**2.d0)/DY(J)                           &
         +                                                                        &   ! d(U2*U3)/dZ
         ( (U2X(I,K+1,J)*DZF(K+1)+ U2X(I,K,J)*DZF(K))/(2.*DZ(K+1))*              &
                 (U3X(I,K+1,J)*DYF(J)+U3X(I,K+1,J-1)*DYF(J-1))/(2.*DY(J))         &
               - (U2X(I,K-1,J)*DZF(K-1) + U2X(I,K,J)*DZF(K))/(2.*DZ(K))*          &
                 (U3X(I,K,J-1)*DYF(J-1)+U3X(I,K,J)*DYF(J))/(2.*DY(J)) )           &
                 /DZF(K)   
          END DO
        END DO
      END DO

      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

      DO J=2,NY
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF2X(I,K,J)=CF2X(I,K,J) - CS1X(I,K,J)
          END DO
        END DO
      END DO

! !C d(U2*U3)/dZ
! 
!       DO J=2,NY
!         DO K=ZSTART,ZEND
!          DO I=0,NXP
! 
!       S1X(I,K,J)=( (U2X(I,K+1,J)*DZF(K+1)+ U2X(I,K,J)*DZF(K))/(2.*DZ(K+1))*        &
!                  (U3X(I,K+1,J)*DYF(J)+U3X(I,K+1,J-1)*DYF(J-1))/(2.*DY(J))         &
!                - (U2X(I,K-1,J)*DZF(K-1) + U2X(I,K,J)*DZF(K))/(2.*DZ(K))*          &
!                  (U3X(I,K,J-1)*DYF(J-1)+U3X(I,K,J)*DYF(J))/(2.*DY(J)) )           &
!                  /DZF(K)
! 
!           END DO
!         END DO
!       END DO
!    
! 
! 
! 
!       CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
!       varp(:,:,:) = 0.d0
!       DO I=0,NXM
!        varp(I,:,:)=S1Z(I,:,:)
!       ENDDO
!       CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
!       DO I=0,NKX
!       CS1Z(I,:,:)=cvarp(I,:,:)
!       ENDDO
!       CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
!       CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)
! 
!       DO J=2,NY
!         DO K=ZSTART,ZEND
!           DO I=0,NX2P
!             CF2X(I,K,J)=CF2X(I,K,J) - CS1X(I,K,J)
!           END DO
!         END DO
!       END DO
! !CCC4

! d(U2*U3)/dY +  d(U3*U3)/dZ


      DO J=JSTART,JEND
        DO K=2,NZ
          DO I=0,NXP
            S1X(I,K,J)=( ((DYF(J+1)*U3X(I,K,J+1)              &             ! d(U2*U3)/dY
                     +DYF(J)*U3X(I,K,J))/(2.*DY(J+1)))       &
                     *((DZF(K)*U2X(I,K,J+1)                  &
                     +DZF(K-1)*U2X(I,K-1,J+1))/(2.*DZ(K)))   &
                     -((DYF(J)*U3X(I,K,J)                    &
                     +DYF(J-1)*U3X(I,K,J-1))/(2.*DY(J)))     &
                     *((DZF(K)*U2X(I,K,J)                    &   
                     +DZF(K-1)*U2X(I,K-1,J))/(2.*DZ(K)))  )/DYF(J) &
                     +                                                &     ! d(U3*U3)/dZ
                      (0.25d0*(U3X(I,K,J)+U3X(I,K+1,J))**2.d0         &
         -0.25d0*(U3X(I,K,J)+U3X(I,K-1,J))**2.d0)/DZ(K)
          END DO
        END DO
      END DO

      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

      DO J=JSTART,JEND
        DO K=2,NZ
          DO I=0,NX2P
            CF3X(I,K,J)=CF3X(I,K,J) - CS1X(I,K,J) 
          END DO
        END DO
      END DO

!C d(U3*U3)/dZ
!       
!       DO J=JSTART,JEND
!        DO K=2,NZ
!          DO I=0,NXP
!             S1X(I,K,J)=                                 &
!          (0.25d0*(U3X(I,K,J)+U3X(I,K+1,J))**2.d0         &
!          -0.25d0*(U3X(I,K,J)+U3X(I,K-1,J))**2.d0)/DZ(K)  
!          END DO
!         END DO
!       END DO
! 
!       CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
!       varp(:,:,:) = 0.d0
!       DO I=0,NXM
!        varp(I,:,:)=S1Z(I,:,:)
!       ENDDO
!       CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
!       DO I=0,NKX
!       CS1Z(I,:,:)=cvarp(I,:,:)
!       ENDDO
!       CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
!       CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)
! 
!       DO J=JSTART,JEND
!         DO K=2,NZ
!           DO I=0,NX2P
!             CF3X(I,K,J)=CF3X(I,K,J) - CS1X(I,K,J)
!           END DO
!         END DO
!       END DO

!c      ENDIF

!C -- At this point, we are done computing the nonlinear terms --

!C Finally, Add CFi to CRi

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
           DO I=0,NX2P
            CR1X(I,K,J)=CR1X(I,K,J) + TEMP5 * CF1X(I,K,J)
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=ZSTART,ZEND
           DO I=0,NX2P
            CR2X(I,K,J)=CR2X(I,K,J) + TEMP5 * CF2X(I,K,J)
          END DO
        END DO
      END DO
 
      DO J=JSTART,JEND
        DO K=2,NZ
           DO I=0,NX2P
            CR3X(I,K,J)=CR3X(I,K,J) + TEMP5 * CF3X(I,K,J)
          END DO
        END DO
      END DO

!C Convert RHS terms to physical space
!C made a change in FFT_X_TO_PHYSICAL(CR2,R2,0,NY+1,0,NZ+1)
!C

      CALL REAL_FOURIER_TRANS_R1 (.false.)
      CALL REAL_FOURIER_TRANS_R2 (.false.)
      CALL REAL_FOURIER_TRANS_R3 (.false.)


    
!C Now, build the explicit RHS terms for the passive scalar(s)

      DO N=1,N_TH
! Do for each scalar:

! Compute the nonlinear terms that are present in the explicit term A
! U1*TH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            S1X(I,K,J)=THX(I,K,J,N)*U1X(I,K,J)
          END DO
        END DO
      END DO
      
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
            CFTHX(I,K,J,N)=CFTHX(I,K,J,N) - CIKXP(I) * CS1X(I,K,J)
          END DO
        END DO
      END DO

! d(U3*TH)/dz + d(U2TH)/dy
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
         S1X(I,K,J)=((THX(I,K+1,J,N)*U3X(I,K+1,J) + THX(I,K,J,N)*U3X(I,K+1,J) &       ! d(U3*TH)/dz
         -THX(I,K,J,N)*U3X(I,K,J)-THX(I,K-1,J,N)*U3X(I,K,J))/(2.d0*DZF(K)))  &
                   +                                                      &      !  d(U2TH)/dy
                    ((THX(I,K,J+1,N)*U2X(I,K,J+1) + THX(I,K,J,N)*U2X(I,K,J+1)  &
         -THX(I,K,J,N)*U2X(I,K,J)-THX(I,K,J-1,N)*U2X(I,K,J))/(2.d0*DYF(J))) 
          END DO
        END DO
      END DO

      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
            CFTHX(I,K,J,N)=CFTHX(I,K,J,N) -  CS1X(I,K,J)
          END DO
        END DO
      END DO

! ! U2*TH
!       DO J=JSTART_TH(N),JEND_TH(N)
!         DO K=ZSTART_TH(N),ZEND_TH(N)
!           DO I=0,NXM
!          S1(I,K,J)=((TH(I,K,J+1,N)*U2(I,K,J+1) + TH(I,K,J,N)*U2(I,K,J+1)  &
!          -TH(I,K,J,N)*U2(I,K,J)-TH(I,K,J-1,N)*U2(I,K,J))/(2.d0*DYF(J)))
!           END DO
!         END DO
!       END DO
!       CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)
!       DO J=JSTART_TH(N),JEND_TH(N)
!         DO K=ZSTART_TH(N),ZEND_TH(N)
!           DO I=0,NKX
!             CFTH(I,K,J,N)=CFTH(I,K,J,N) - CS1(I,K,J)
!           END DO
!         END DO
!       END DO

! We are done with the horizontal derivatives of the nonlinear terms

! Add CFTH to the RHS vector CRTH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
            CRTHX(I,K,J,N)=CRTHX(I,K,J,N) + TEMP5 * CFTHX(I,K,J,N)
          END DO
        END DO
      END DO
! Done with computation of RHS, explicit terms for the THETA equation
! Transform back to physical space
      END DO   ! end of the scalar variables 

      CALL REAL_FOURIER_TRANS_Rth (.false.)


!C Before going to  ADI1 we have to store RHS after subtracting previous TH_n 
!C at previous time step to bulid new RHS for  ADI2 based on TH_n+1/2 at intermidiate 

       DO N=1,N_TH
! ! Do for each scalar:

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            S1X(I,K,J)= RTHX(I,K,J,N)-THX(I,K,J,N)
          END DO
        END DO
      END DO


!C Compute the vertical viscous term in physical space and add to RHS
!C  at  ADI1 step .... Important in this step we have to add TH_n
!C at previous time step so that after multiplying a factor 1/2 into the 
!C RHS side during ADI1 operation RHS will be taken care TH_n not 1/2*TH_n
 

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            RTHX(I,K,J,N)=0.5*S1X(I,K,J) + THX(I,K,J,N) +  (TEMP1/PR(N))  &
              *(  ((THX(I,K+1,J,N) - THX(I,K,J,N)) / DZ(K+1)             &
                -(THX(I,K,J,N)   - THX(I,K-1,J,N)) / DZ(K)) /DZF(K)  )
          END DO
        END DO
      END DO

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
       DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            RTHX(I,K,J,N)=RTHX(I,K,J,N) +  TEMP2  &
              *(( 0.5d0*(KAPPA_T(I,K+1,J,N)+KAPPA_T(I,K,J,N))*(THX(I,K+1,J,N) - THX(I,K,J,N)) / DZ(K+1)             &
                -  0.5d0*(KAPPA_T(I,K,J,N)+KAPPA_T(I,K-1,J,N))*(THX(I,K,J,N)   - THX(I,K-1,J,N)) / DZ(K)) /DZF(K)  )
          END DO
        END DO
      END DO
      ENDIF 
!C Solve for U1
!C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXP
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO 

!C Build the implicit system of equations for U1 
      DO K=ZSTART_TH(N)-1,ZEND_TH(N)+1

        DO J=JSTART_TH(N),JEND_TH(N)
          DO I=0,NXP
            MATL(I,J)=-(TEMP1/PR(N))/(DY(J)*DYF(J))
            MATD(I,J)=1.-(TEMP1/PR(N))*(-1./(DY(J+1)*DYF(J))  &
            -1./(DY(J)*DYF(J))) 
            MATU(I,J)=-(TEMP1/PR(N))/(DY(J+1)*DYF(J))
            VEC(I,J)=RTHX(I,K,J,N)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly

!        IF (LES) THEN
!        DO J=JSTART,JEND

         IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
          DO J=JSTART_TH(N),JEND_TH(N)
          DO I=0,NXP
            MATL(I,J) = MATL(I,J)-TEMP2*0.5d0*(KAPPA_T(I,K,J,N) + KAPPA_T(I,K,J-1,N))/(DY(J)*DYF(J))
            MATD(I,J) = MATD(I,J) + TEMP2*(0.5d0*(KAPPA_T(I,K,J,N) + KAPPA_T(I,K,J+1,N))/(DY(J+1)*DYF(J))  &
            +0.5d0*(KAPPA_T(I,K,J,N) + KAPPA_T(I,K,J-1,N))/(DY(J)*DYF(J))) 
            MATU(I,J)= MATU(I,J) -TEMP2*0.5d0*(KAPPA_T(I,K,J,N) + KAPPA_T(I,K,J+1,N))/(DY(J+1)*DYF(J))
          END DO
          END DO


!c          DO J=J1,J2
!          DO I=0,NXM
!            MATL(I,J) = MATL(I,J) - TEMP2 * NU_T(I,K,J) 
!     &                               / (DY(J)*DYF(J))
!            MATD(I,J) = MATD(I,J) + TEMP2 * NU_T(I,K,J+1)
!     &                              / (DY(J+1)*DYF(J))
!     &                            + TEMP2 * NU_T(I,K,J)
!     &                              / (DY(J)*DYF(J))
!            MATU(I,J) = MATU(I,J) - TEMP2 * NU_T(I,K,J+1)
!     &                             / (DY(J+1)*DYF(J))
!          END DO
!        END DO
        END IF


!C Set the boundary conditions for U1
          CALL APPLY_BC_TH_LOWER(N,K)
          CALL APPLY_BC_TH_UPPER(N,K)
!C Now, solve the tridiagonal system for U1(:,k,:)
          CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXP)

        DO J=JSTART_TH(N)-1,JEND_TH(N)+1
          DO I=0,NXP
            THX(I,K,J,N)=VEC(I,J)
          END DO
        END DO
! End do k
      END DO



!C Compute the horizontal viscous term in physical space and add to new RHS
!C already store in S1 
!C Important we have to multiply 1/2 with S1 before adding to RHS due
!C half time splitting during AD1

     
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            RTHX(I,K,J,N)=0.5*S1X(I,K,J) + THX(I,K,J,N) + (TEMP1/PR(N))  &
              *(  ((THX(I,K,J+1,N) - THX(I,K,J,N)) / DY(J+1)            &
              -(THX(I,K,J,N)   - THX(I,K,J-1,N)) / DY(J)) /DYF(J)  )

          END DO
        END DO
      END DO

     IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
     DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            RTHX(I,K,J,N)=RTHX(I,K,J,N) + TEMP2  &
              *(  (0.5d0*(KAPPA_T(I,K,J,N) + KAPPA_T(I,K,J+1,N))*(THX(I,K,J+1,N) - THX(I,K,J,N)) / DY(J+1)            &
              -0.5d0*(KAPPA_T(I,K,J,N) + KAPPA_T(I,K,J-1,N))*(THX(I,K,J,N)   - THX(I,K,J-1,N)) / DY(J)) /DYF(J)  )

          END DO
        END DO
      END DO 
      ENDIF

! Initialize the matrix used to store implicit coefficients
      DO K=0,NZ+1
        DO I=0,NXP
          MATL_Z(I,K)=0.
          MATD_Z(I,K)=1.
          MATU_Z(I,K)=0.
          VEC_Z(I,K)=0.
        END DO
      END DO 

!C Build the implicit system of equations for TH 

      DO J=JSTART_TH(N),JEND_TH(N)

        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            MATL_Z(I,K)=-(TEMP1/PR(N))/(DZ(K)*DZF(K))
            MATD_Z(I,K)=1.-(TEMP1/PR(N))*(-1./(DZ(K+1)*DZF(K)) &
              -1./(DZ(K)*DZF(K))) 
            MATU_Z(I,K)=-(TEMP1/PR(N))/(DZ(K+1)*DZF(K))
            VEC_Z(I,K)=RTHX(I,K,J,N)
          END DO
        END DO
! IF using a subgrid model (LES) then add the eddy viscosity part implicitly

        IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
         DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            MATL_Z(I,K)=MATL_Z(I,K) - TEMP2*0.5d0*(KAPPA_T(I,K-1,J,N)+KAPPA_T(I,K,J,N))/(DZ(K)*DZF(K))
            MATD_Z(I,K)=MATD_Z(I,K) + TEMP2*(0.5d0*(KAPPA_T(I,K+1,J,N)+KAPPA_T(I,K,J,N))/(DZ(K+1)*DZF(K)) &
              + 0.5d0*(KAPPA_T(I,K-1,J,N)+KAPPA_T(I,K,J,N))/(DZ(K)*DZF(K))) 
            MATU_Z(I,K)=MATU_Z(I,K) - TEMP2*0.5d0*(KAPPA_T(I,K+1,J,N)+KAPPA_T(I,K,J,N))/(DZ(K+1)*DZF(K))
          END DO
         END DO
        ENDIF
!        IF (LES) THEN
!        DO J=JSTART,JEND

!c         IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
!c          DO J=J1,J2
!          DO I=0,NXM
!            MATL(I,J) = MATL(I,J) - TEMP2 * NU_T(I,K,J) 
!     &                               / (DY(J)*DYF(J))
!            MATD(I,J) = MATD(I,J) + TEMP2 * NU_T(I,K,J+1)
!     &                              / (DY(J+1)*DYF(J))
!     &                            + TEMP2 * NU_T(I,K,J)
!     &                              / (DY(J)*DYF(J))
!            MATU(I,J) = MATU(I,J) - TEMP2 * NU_T(I,K,J+1)
!     &                             / (DY(J+1)*DYF(J))
!          END DO
!        END DO
!        END IF


!C Set the boundary conditions for U1
          CALL APPLY_BC_TH_LEFT(N,J)
          CALL APPLY_BC_TH_RIGHT(N,J)
!C Now, solve the tridiagonal system for TH(:,k,:)
        IF ( TH_BC_ZMAX(N)  .EQ. 5 ) THEN
         D=1.0/DZ(NZ)
         CALL THOMAS_REAL_SP(MATL_Z,MATD_Z,MATU_Z,VEC_Z,D,NZ,NZ+1,NXP)
        ELSE   
         CALL THOMAS_REAL(MATL_Z,MATD_Z,MATU_Z,VEC_Z,NZ+1,NXP)
        ENDIF
 
        DO K=ZSTART_TH(N)-1,ZEND_TH(N)+1
          DO I=0,NXP
            THX(I,K,J,N)=VEC_Z(I,K)
          END DO
        END DO
! End do J
      END DO


! End do number of passive scalars
        END DO





!C Before going to  ADI1 we have to store RHS after subtracting previous  U_n 
!C at previous time step to bulid new RHS for  ADI2 based on U_n+1/2 at intermidiate 



      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            S1X(I,K,J)= R1X(I,K,J)-U1X(I,K,J)
          END DO
        END DO
      END DO


!C Compute the vertical viscous term in physical space and add to RHS
!C  at  ADI1 step .... Important in this step we have to add U_n
!C at previous time step so that after multiplying a factor 1/2 into the 
!C RHS side during ADI1 operation RHS will be taken care U_n not 1/2*U_n
 

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R1X(I,K,J)=0.5*S1X(I,K,J) + U1X(I,K,J) +  TEMP1*  &
               (  ((U1X(I,K+1,J) - U1X(I,K,J)) / DZ(K+1)     &
                -(U1X(I,K,J)   - U1X(I,K-1,J)) / DZ(K)) /DZF(K)  )
          END DO
        END DO
      END DO

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R1X(I,K,J)=R1X(I,K,J) +  TEMP2*  &
               ((0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J))*(U1X(I,K+1,J) - U1X(I,K,J)) / DZ(K+1)     &
                -0.5d0*(NU_T(I,K-1,J)+NU_T(I,K,J))*(U1X(I,K,J)   - U1X(I,K-1,J)) / DZ(K)) /DZF(K)  )
          END DO
        END DO
      END DO
      ENDIF

!C Solve for U1
!C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXP
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO 

!C Build the implicit system of equations for U1 
      DO K=ZSTART,ZEND

        DO J=JSTART,JEND
          DO I=0,NXP
            MATL(I,J)=-TEMP1/(DY(J)*DYF(J))
            MATD(I,J)=1.-TEMP1*(-1./(DY(J+1)*DYF(J))  &
             -1./(DY(J)*DYF(J))) 
            MATU(I,J)=-TEMP1/(DY(J+1)*DYF(J))
            VEC(I,J)=R1X(I,K,J)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
        IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
         DO J=JSTART,JEND
          DO I=0,NXP
            MATL(I,J)=MATL(I,J) - TEMP2*0.5d0*(NU_T(I,K,J-1)+NU_T(I,K,J))/(DY(J)*DYF(J))
            MATD(I,J)=MATD(I,J) + TEMP2*(0.5d0*(NU_T(I,K,J+1)+NU_T(I,K,J))/(DY(J+1)*DYF(J))  &
             +0.5d0*(NU_T(I,K,J-1)+NU_T(I,K,J))/(DY(J)*DYF(J))) 
            MATU(I,J)=MATU(I,J) - TEMP2*0.5d0*(NU_T(I,K,J+1)+NU_T(I,K,J))/(DY(J+1)*DYF(J))
          END DO
        END DO
       ENDIF

!        IF (LES) THEN
!        DO J=JSTART,JEND

!c         IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
!c          DO J=J1,J2
!          DO I=0,NXM
!            MATL(I,J) = MATL(I,J) - TEMP2 * NU_T(I,K,J) 
!     &                               / (DY(J)*DYF(J))
!            MATD(I,J) = MATD(I,J) + TEMP2 * NU_T(I,K,J+1)
!     &                              / (DY(J+1)*DYF(J))
!     &                            + TEMP2 * NU_T(I,K,J)
!     &                              / (DY(J)*DYF(J))
!            MATU(I,J) = MATU(I,J) - TEMP2 * NU_T(I,K,J+1)
!     &                             / (DY(J+1)*DYF(J))
!          END DO
!        END DO
!        END IF


!C Else, we are running in serial mode
!C Set the boundary conditions for U1
          CALL APPLY_BC_1_LOWER(K)
          CALL APPLY_BC_1_UPPER(K)
!C Now, solve the tridiagonal system for U1(:,k,:)
          CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXP)

        DO J=JSTART-1,JEND+1
          DO I=0,NXP
            U1X(I,K,J)=VEC(I,J)
          END DO
        END DO
! End do k
      END DO



!C Compute the horizontal viscous term in physical space and add to new RHS
!C already store in S1 
!C Important we have to multiply 1/2 with S1 before adding to RHS due
!C half time splitting during AD1

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R1X(I,K,J)=0.5*S1X(I,K,J) + U1X(I,K,J) + TEMP1*        & 
            (  ((U1X(I,K,J+1) - U1X(I,K,J)) / DY(J+1)             &
              -(U1X(I,K,J)   - U1X(I,K,J-1)) / DY(J)) /DYF(J)  )

          END DO
        END DO
      END DO

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
       DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R1X(I,K,J)=R1X(I,K,J) + TEMP2*        & 
            (( 0.5d0*(NU_T(I,K,J+1)+NU_T(I,K,J))*(U1X(I,K,J+1) - U1X(I,K,J)) / DY(J+1)             &
              -0.5d0*(NU_T(I,K,J-1)+NU_T(I,K,J))*(U1X(I,K,J)   - U1X(I,K,J-1)) / DY(J)) /DYF(J)  )

          END DO
        END DO
      END DO
     ENDIF 

! Initialize the matrix used to store implicit coefficients
      DO K=0,NZ+1
        DO I=0,NXP
          MATL_Z(I,K)=0.
          MATD_Z(I,K)=1.
          MATU_Z(I,K)=0.
          VEC_Z(I,K)=0.
        END DO
      END DO 

!C Build the implicit system of equations for U1 

      DO J=JSTART,JEND

        DO K=ZSTART,ZEND
          DO I=0,NXP
            MATL_Z(I,K)=-TEMP1/(DZ(K)*DZF(K))
            MATD_Z(I,K)=1.-TEMP1*(-1./(DZ(K+1)*DZF(K))   &
              -1./(DZ(K)*DZF(K))) 
            MATU_Z(I,K)=-TEMP1/(DZ(K+1)*DZF(K))
            VEC_Z(I,K)=R1X(I,K,J)
          END DO
        END DO
! IF using a subgrid model (LES) then add the eddy viscosity part implicitly

        IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN          
         DO K=ZSTART,ZEND
          DO I=0,NXP
            MATL_Z(I,K)=MATL_Z(I,K) -TEMP2* 0.5d0*(NU_T(I,K-1,J)+NU_T(I,K,J))/(DZ(K)*DZF(K))
            MATD_Z(I,K)=MATD_Z(I,K) + TEMP2*(0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J))/(DZ(K+1)*DZF(K))   &
              + 0.5d0*(NU_T(I,K-1,J)+NU_T(I,K,J))/(DZ(K)*DZF(K))) 
            MATU_Z(I,K)=MATU_Z(I,K)-TEMP2*0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J))/(DZ(K+1)*DZF(K))
          END DO
         END DO
        ENDIF

!        IF (LES) THEN
!        DO J=JSTART,JEND

!c         IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
!c          DO J=J1,J2
!          DO I=0,NXM
!            MATL(I,J) = MATL(I,J) - TEMP2 * NU_T(I,K,J) 
!     &                               / (DY(J)*DYF(J))
!            MATD(I,J) = MATD(I,J) + TEMP2 * NU_T(I,K,J+1)
!     &                              / (DY(J+1)*DYF(J))
!     &                            + TEMP2 * NU_T(I,K,J)
!     &                              / (DY(J)*DYF(J))
!            MATU(I,J) = MATU(I,J) - TEMP2 * NU_T(I,K,J+1)
!     &                             / (DY(J+1)*DYF(J))
!          END DO
!        END DO
!        END IF


!C Set the boundary conditions for U1
          CALL APPLY_BC_1_LEFT(J)
          CALL APPLY_BC_1_RIGHT(J)
!C Now, solve the tridiagonal system for U1(:,k,:)
          CALL THOMAS_REAL(MATL_Z,MATD_Z,MATU_Z,VEC_Z,NZ+1,NXP)

        DO K=ZSTART-1,ZEND+1
          DO I=0,NXP
            U1X(I,K,J)=VEC_Z(I,K)
          END DO
        END DO
! End do J
      END DO


!C Before going to  ADI1 we have to store RHS after subtracting previous  U_n 
!C at previous time step to bulid new RHS for  ADI2 based on U_n+1/2 at intermidiate 

      DO J=2,NY
        DO K=ZSTART,ZEND
          DO I=0,NXP
            S1X(I,K,J)= R2X(I,K,J)-U2X(I,K,J)
          END DO
        END DO
      END DO


!C Compute the vertical viscous term in physical space and add to RHS
!C  at  ADI1 step .... Important in this step we have to add U_n
!C at previous time step so that after multiplying a factor 1/2 into the 
!C RHS side during ADI1 operation RHS will be taken care U_n not 1/2*U_n

      DO J=2,NY
        DO K=ZSTART,ZEND
          DO I=0,NXP
             R2X(I,K,J)=0.5*S1X(I,K,J) + U2X(I,K,J) + TEMP1*            &
             (  ((U2X(I,K,J+1) - U2X(I,K,J))  / DYF(J)                 & 
                -(U2X(I,K,J)   - U2X(I,K,J-1))/ DYF(J-1))/DY(J)  )
           ENDDO
         ENDDO
      ENDDO

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
       DO J=2,NY
        DO K=ZSTART,ZEND
          DO I=0,NXP
             R2X(I,K,J)=R2X(I,K,J) + TEMP2*            &
             (( NU_T(I,K,J)*(U2X(I,K,J+1) - U2X(I,K,J))  / DYF(J)                 & 
                -NU_T(I,K,J-1)*(U2X(I,K,J)   - U2X(I,K,J-1))/ DYF(J-1))/DY(J)  )
           ENDDO
         ENDDO
      ENDDO
      ENDIF

! Initialize the matrix used to store implicit coefficients
      DO K=0,NZ+1
        DO I=0,NXP
          MATL_Z(I,K)=0.
          MATD_Z(I,K)=1.
          MATU_Z(I,K)=0.
          VEC_Z(I,K)=0.
        END DO
      END DO

      DO J=2,NY

        DO K=ZSTART,ZEND
          DO I=0,NXP
            MATL_Z(I,K)=-TEMP1/(DZ(K)*DZF(K))
            MATD_Z(I,K)=1.-TEMP1*(-1./(DZ(K+1)*DZF(K))   &
              -1./(DZ(K)*DZF(K)))
            MATU_Z(I,K)=-TEMP1/(DZ(K+1)*DZF(K))
            VEC_Z(I,K)=R2X(I,K,J)
          END DO
        END DO
        
       IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
       	DO K=ZSTART,ZEND
          DO I=0,NXP
            MATL_Z(I,K)=MATL_Z(I,K) -TEMP2*0.250d0*(NU_T(I,K-1,J)+NU_T(I,K-1,J-1)+NU_T(I,K,J)+NU_T(I,K,J-1))/(DZ(K)*DZF(K))
            MATD_Z(I,K)=MATD_Z(I,K) + TEMP2*(0.250d0*(NU_T(I,K+1,J)+NU_T(I,K+1,J-1)+NU_T(I,K,J)+NU_T(I,K,J-1))/(DZ(K+1)*DZF(K))   &
              + 0.250d0*(NU_T(I,K-1,J)+NU_T(I,K-1,J-1)+NU_T(I,K,J)+NU_T(I,K,J-1))/(DZ(K)*DZF(K)))
            MATU_Z(I,K)=MATU_Z(I,K)-TEMP2*0.250d0*(NU_T(I,K+1,J)+NU_T(I,K+1,J-1)+NU_T(I,K,J)+NU_T(I,K,J-1))/(DZ(K+1)*DZF(K))
          END DO
        END DO          	       
       ENDIF	       

!       IF (LES) THEN
!c       IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
!        DO J=2,NY
!          DO I=0,NXM
!            MATL(I,J) = MATL(I,J) 
!     &      - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))/(DYF(J-1)*DY(J))
!            MATD(I,J) = MATD(I,J) 
!     &      + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))/(DYF(J)*DY(J))
!     &      + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))/(DYF(J-1)*DY(J))
!            MATU(I,J) = MATU(I,J) 
!     &      - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))/(DYF(J)*DY(J))
!          END DO
!        END DO
!        END IF


!C Set the boundary conditions for U2
          CALL APPLY_BC_2_LEFT(J)
          CALL APPLY_BC_2_RIGHT(J)

!C Now, solve the tridiagonal system for U2(i,:,k)
        IF ( V_BC_ZMAX.EQ.5 ) THEN
         D=1.0/DZ(NZ)
         CALL THOMAS_REAL_SP(MATL_Z,MATD_Z,MATU_Z,VEC_Z,D,NZ,NZ+1,NXP)
        ELSE
         CALL THOMAS_REAL(MATL_Z,MATD_Z,MATU_Z,VEC_Z,NZ+1,NXP)
        ENDIF


        DO K=ZSTART-1,ZEND+1
          DO I=0,NXP
            U2X(I,K,J)=VEC_Z(I,K)
          END DO
        END DO
! End do J
      END DO


      DO J=2,NY
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R2X(I,K,J)=0.5*S1X(I,K,J) + U2X(I,K,J) +  TEMP1*   &
              (  ((U2X(I,K+1,J) - U2X(I,K,J)) / DZ(K+1)       &
                -(U2X(I,K,J)   - U2X(I,K-1,J)) / DZ(K)) /DZF(K)  )
          END DO
        END DO
      END DO

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
      DO J=2,NY
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R2X(I,K,J)=R2X(I,K,J) +  TEMP2*   &
              (( 0.250d0*(NU_T(I,K+1,J)+NU_T(I,K+1,J-1)+NU_T(I,K,J)+NU_T(I,K,J-1))*(U2X(I,K+1,J) - U2X(I,K,J)) / DZ(K+1)       &
                -0.250d0*(NU_T(I,K-1,J)+NU_T(I,K-1,J-1)+NU_T(I,K,J)+NU_T(I,K,J-1))*(U2X(I,K,J)   - U2X(I,K-1,J)) / DZ(K)) /DZF(K)  )
          END DO
        END DO
      END DO
      
      ENDIF  	      
!C Solve for U2
!C Note, here the matrix will be indexed from 1...NY+1 corresponding to U2(0:NY)

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXP
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO

!C Build implicit matrix for U2
      DO K=ZSTART,ZEND

        DO J=2,NY
          DO I=0,NXP
            MATL(I,J)= -TEMP1/(DYF(J-1)*DY(J))
            MATD(I,J)=1.+TEMP1/(DYF(J)*DY(J)) + TEMP1/(DYF(J-1)*DY(J))
            MATU(I,J)= -TEMP1/(DYF(J)*DY(J))
            VEC(I,J)=R2X(I,K,J)
          END DO
        END DO
        
       IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
        DO J=2,NY
          DO I=0,NXP
            MATL(I,J)= MATL(I,J) -TEMP2*NU_T(I,K,J-1)/(DYF(J-1)*DY(J))
            MATD(I,J)= MATD(I,J) +TEMP2*(NU_T(I,K,J)/(DYF(J)*DY(J)) + NU_T(I,K,J-1)/(DYF(J-1)*DY(J)))
            MATU(I,J)= MATU(I,J) -TEMP2*NU_T(I,K,J)/(DYF(J)*DY(J))
          END DO
        END DO       
       ENDIF     

!       IF (LES) THEN
!c       IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
!        DO J=2,NY
!          DO I=0,NXM
!            MATL(I,J) = MATL(I,J) 
!     &      - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))/(DYF(J-1)*DY(J))
!            MATD(I,J) = MATD(I,J) 
!     &      + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))/(DYF(J)*DY(J))
!     &      + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))/(DYF(J-1)*DY(J))
!            MATU(I,J) = MATU(I,J) 
!     &      - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))/(DYF(J)*DY(J))
!          END DO
!        END DO
!        END IF


!C Set the boundary conditions for U2
          CALL APPLY_BC_2_LOWER(K)
          CALL APPLY_BC_2_UPPER(K)

!C Now, solve the tridiagonal system for U2(i,:,k)
          CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXP)

        DO J=1,NY+1
          DO I=0,NXP
            U2X(I,K,J)=VEC(I,J)
          END DO
        ENDDO
! End do k
      END DO

      

 
!C Before going to  ADI1 we have to store RHS after subtracting previous  U_n 
!C at previous time step to bulid new RHS for  ADI2 based on U_n+1/2 at intermidiate 

      DO J=JSTART,JEND
        DO K=2,NZ
          DO I=0,NXP
            S1X(I,K,J)= R3X(I,K,J)-U3X(I,K,J)
          END DO
        END DO
      END DO





      DO J=JSTART,JEND
        DO K=2,NZ
          DO I=0,NXP
             R3X(I,K,J)=0.5*S1X(I,K,J) + U3X(I,K,J) + TEMP1*            &
            (  ((U3X(I,K,J+1) - U3X(I,K,J)) / DY(J+1)                  &
                -(U3X(I,K,J)   - U3X(I,K,J-1)) / DY(J)) /DYF(J)  )
           END DO
        END DO
      END DO

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
       DO J=JSTART,JEND
        DO K=2,NZ
          DO I=0,NXP
             R3X(I,K,J)= R3X(I,K,J) + TEMP2*            &
            (( 0.25*(NU_T(I,K,J+1)+NU_T(I,K-1,J+1)+NU_T(I,K,J)+NU_T(I,K-1,J))*(U3X(I,K,J+1) - U3X(I,K,J)) / DY(J+1)  &
              -0.25*(NU_T(I,K,J-1)+NU_T(I,K-1,J-1)+NU_T(I,K,J)+NU_T(I,K-1,J))*(U3X(I,K,J)   - U3X(I,K,J-1)) / DY(J)) &
                    /DYF(J)  )
           END DO
        END DO
      END DO
      
      ENDIF   
      

      DO K=0,NZ+1
        DO I=0,NXP
          MATL_Z(I,K)=0.
          MATD_Z(I,K)=1.
          MATU_Z(I,K)=0.
          VEC_Z(I,K)=0.
        END DO
      END DO


!C Build implicit matrix for U3
      DO J=JSTART-1,JEND+1

        DO K=2,NZ
          DO I=0,NXP
            MATL_Z(I,K)= -TEMP1/(DZF(K-1)*DZ(K))
            MATD_Z(I,K)=1.+TEMP1/(DZF(K)*DZ(K))+TEMP1/(DZF(K-1)*DZ(K))
            MATU_Z(I,K)= -TEMP1/(DZF(K)*DZ(K))
            VEC_Z(I,K)=R3X(I,K,J)
          END DO
        END DO

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
      	DO K=2,NZ
          DO I=0,NXP
            MATL_Z(I,K)= MATL_Z(I,K) - TEMP2*NU_T(I,K-1,J)/(DZF(K-1)*DZ(K))
            MATD_Z(I,K)= MATD_Z(I,K) + TEMP2*(NU_T(I,K,J)/(DZF(K)*DZ(K))+NU_T(I,K-1,J)/(DZF(K-1)*DZ(K)))
            MATU_Z(I,K)= MATU_Z(I,K) - TEMP2*NU_T(I,K,J)/(DZF(K)*DZ(K))
          END DO
        END DO      
      ENDIF	      

!C Set the boundary conditions for U3
          CALL APPLY_BC_3_LEFT(J)
          CALL APPLY_BC_3_RIGHT(J)
!C Now, solve the tridiagonal system for U3(i,:,k)
         IF ( W_BC_ZMAX.EQ. 5 ) THEN 
          D=1.0/DZF(NZ-1)   
          CALL THOMAS_REAL_SP(MATL_Z,MATD_Z,MATU_Z,VEC_Z,D,NZ,NZ+1,NXP)
         ELSE
          CALL THOMAS_REAL(MATL_Z,MATD_Z,MATU_Z,VEC_Z,NZ+1,NXP)
         ENDIF


        DO K=1,NZ+1
          DO I=0,NXP
            U3X(I,K,J)=VEC_Z(I,K)
          END DO
        END DO
! End do J
      END DO
          

!C Compute the vertical viscous term in physical space and add to RHS
!C  at  ADI1 step .... Important in this step we have to add U3_n
!C at previous time step so that after multiplying a factor 1/2 into the 
!C RHS side during ADI1 operation RHS will be taken care U3_n not 1/2*U_n


      DO J=JSTART,JEND
        DO K=2,NZ
          DO I=0,NXP
            R3X(I,K,J)=0.5*S1X(I,K,J) + U3X(I,K,J) + TEMP1*   &
             (  ((U3X(I,K+1,J) - U3X(I,K,J))  / DZF(K)       &
                -(U3X(I,K,J)   - U3X(I,K-1,J))/ DZF(K-1))/DZ(K)  )
          END DO
        END DO
      END DO
      
     IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
      DO J=JSTART,JEND
        DO K=2,NZ
          DO I=0,NXP
            R3X(I,K,J)=R3X(I,K,J) + TEMP2*   &
             (( NU_T(I,K,J)*(U3X(I,K+1,J) - U3X(I,K,J))  / DZF(K)       &
                -NU_T(I,K-1,J)*(U3X(I,K,J)   - U3X(I,K-1,J))/ DZF(K-1))/DZ(K)  )
          END DO
        END DO
      END DO
     ENDIF 

!C Build the implicit system of equations for U3 

      DO J=0,NY+1
        DO I=0,NXP
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO      

!C Build the implicit system of equations for U3
      DO K=2,NZ

        DO J=JSTART,JEND
          DO I=0,NXP
            MATL(I,J)=-TEMP1/(DY(J)*DYF(J))
            MATD(I,J)=1.-TEMP1*(-1./(DY(J+1)*DYF(J))  &
              -1./(DY(J)*DYF(J))) 
            MATU(I,J)=-TEMP1/(DY(J+1)*DYF(J))
            VEC(I,J)=R3X(I,K,J)
          END DO
        END DO

       IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
       	 DO J=JSTART,JEND
          DO I=0,NXP
            MATL(I,J)=MATL(I,J) - TEMP2*0.25*(NU_T(I,K,J-1)+NU_T(I,K-1,J-1)+NU_T(I,K,J)+NU_T(I,K-1,J))/(DY(J)*DYF(J))
            MATD(I,J)=MATD(I,J) + TEMP2*(0.25*(NU_T(I,K,J+1)+NU_T(I,K-1,J+1)+NU_T(I,K,J)+NU_T(I,K-1,J))/(DY(J+1)*DYF(J))  &
              +0.25*(NU_T(I,K,J-1)+NU_T(I,K-1,J-1)+NU_T(I,K,J)+NU_T(I,K-1,J))/(DY(J)*DYF(J))) 
            MATU(I,J)=MATU(I,J) - TEMP2*0.25*(NU_T(I,K,J+1)+NU_T(I,K-1,J+1)+NU_T(I,K,J)+NU_T(I,K-1,J))/(DY(J+1)*DYF(J))
          END DO
        END DO      
       ENDIF	       
! IF using a subgrid model (LES) then add the eddy viscosity part implicitly

!        IF (LES) THEN
!        DO J=JSTART,JEND

!c         IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.10))) THEN
!c          DO J=J1,J2
!          DO I=0,NXM
!            MATL(I,J) = MATL(I,J) - TEMP2 * NU_T(I,K,J) 
!     &                               / (DY(J)*DYF(J))
!            MATD(I,J) = MATD(I,J) + TEMP2 * NU_T(I,K,J+1)
!     &                              / (DY(J+1)*DYF(J))
!     &                            + TEMP2 * NU_T(I,K,J)
!     &                              / (DY(J)*DYF(J))
!            MATU(I,J) = MATU(I,J) - TEMP2 * NU_T(I,K,J+1)
!     &                             / (DY(J+1)*DYF(J))
!          END DO
!        END DO
!        END IF


!C Set the boundary conditions for U3
          CALL APPLY_BC_3_LOWER(K)
          CALL APPLY_BC_3_UPPER(K)
!C Now, solve the tridiagonal system for U3(:,k,:)
          CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXP)

        DO J=JSTART-1,JEND+1
          DO I=0,NXP
            U3X(I,K,J)=VEC(I,J)
          END DO
        END DO
! End do k
      END DO
    

!C If Variable timestepping and done with one full R-K step, update
!C DELTA_T based on the specified CFL number
!C This is not parallelized and should be used only in the serial
!C version to ensure that each process uses the same timestep
!       IF ((.NOT.USE_MPI).and.(VARIABLE_DT).and.(RK_STEP.eq.3)  &
!             .and.(MOD(TIME_STEP,UPDATE_DT).EQ.0)) THEN

     IF ((VARIABLE_DT).and.(RK_STEP.eq.3) &
             .and.(MOD(TIME_STEP,UPDATE_DT).EQ.0)) THEN  

        IF (USE_MPI)THEN
          CALL COURANT_DUCT_MPI
         ELSE
          CALL COURANT_DUCT
         ENDIF
      END IF 

      ! Transform TH and U to Fourier Space 
      CALL REAL_FOURIER_TRANS_U1 (.true.)
      CALL REAL_FOURIER_TRANS_U2 (.true.) 
      CALL REAL_FOURIER_TRANS_U3 (.true.) 

      call allocation_R1 (.false.)
      call allocation_R2 (.false.)
      call allocation_R3 (.false.)  

!C      CALL FFT_X_TO_FOURIER(U1,CU1,0,NY+1,0,NZ+1)
!C      CALL FFT_X_TO_FOURIER(U2,CU2,0,NY+1,0,NZ+1)
!C      CALL FFT_X_TO_FOURIER(U3,CU3,0,NY+1,0,NZ+1)
      IF (N_TH .gt. 0 ) THEN
       CALL REAL_FOURIER_TRANS_TH (.true.)
       CALL allocation_Rth (.false.) 
!C        CALL FFT_X_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N),0,NY+1,0,NZ+1)
      ENDIF



!C Begin second step of the Fractional Step algorithm, making u divergence free
!C The following subroutine projects Uhat onto divergence free space
      CALL REM_DIV_DUCT
  
 
!C Now, phi is stored in CR1, use this to update the pressure field
!C Note, here we divide by H_BAR since it was absorbed into PHI in REM_DIV
      DO J=JSTART-1,JEND+1
        DO K=ZSTART-1,ZEND+1
          DO I=0,NX2P
            CPX(I,K,J)=CPX(I,K,J)+CR1X(I,K,J)/TEMP4
          END DO
        END DO
      END DO


      RETURN
      END

!C--.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_DUCT
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!      INCLUDE 'header_duct'

use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG, BC, CALL_CHEEK_DIV
use pr_rem
use mpi_var      
implicit none


!       REAL(r8) RTMP(0:NX+1,0:NZ+1,0:NY+1)
!       COMPLEX*16 CRTMP(0:NX/2,0:NZ+1,0:NY+1)
!       EQUIVALENCE (RTMP,CRTMP)

      REAl*8 DIV  
      INTEGER I,J,K

!c     Set BCs for phi
!c     bc(1) : left, bc(2) : right, bc(3) : bottom, bc(4) : top
!c     For periodic and dirchlet set RHS to corresponding values
!c     For neumann, set RHS to the derivative times the grid spacing
!c     normal to the wall
      bc(:)=1

      CR1X(:,:,:)  = (0.0d0,0.0d0)
      CS1X(:,:,:)  = (0.0d0,0.0d0) 
      CALL_CHEEK_DIV = .false.
 

!C Now, create the RHS vector
      DO J=2,NYM
       DO K=2,NZM
        DO I=0,NX2P
           CR1X(I,K,J)=CIKXP(I)*CU1X(I,K,J)               &
                + (CU2X(I,K,J+1)-CU2X(I,K,J))/DYF(J)     &
                + (CU3X(I,K+1,J)-CU3X(I,K,J))/DZF(K)
           CS1X(I,K,J)=CR1X(I,K,J)
        END DO
       END DO
      END DO

      IF (CALL_CHEEK_DIV) THEN 

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


      DIV = 0.0D0
      DO J = 2,NYM
       DO K = 2,NZM
        DO I = 1,NXP
         DIV = DIV + DABS(S1X(I,K,J))
        ENDDO
       ENDDO
      ENDDO
      DIV = DIV/(DBLE(NXP*NP)*DBLE(NY)*DBLE(NZ))
      CALL MPI_COMBINE_STATS(DIV,1,1)

!C      CALL FFT_X_TO_FOURIER(RTMP,CRTMP,0,NY+1,0,NZ+1)
      IF (rank .eq. 0) THEN

      write(6,*) "THE DIVERGENCE IS ", DIV, &
                " BEFORE REMOVING THE DIVERGENCE"
     
      ENDIF
      ENDIF



      IF ( INIT_FLAG ) THEN
       DO I = 0,NX2P
        P_TMP(:,:) = 0.0d0; RHS_TMP(:,:) = 0.0d0
        DO K = 2,NZM
         DO J = 2,NYM
          RHS_TMP(K,J) = dble(CR1X(I,K,J))
!c           RHS_TMP(K,J) =1.0
         ENDDO
        ENDDO
        CALL MULTIGRID(P_TMP,RHS_TMP,I)
       ENDDO
      ENDIF
      INIT_FLAG = .FALSE.

      DO I = 0,NX2P
!c     Solve for complex part of phi
       RHS_TMP(:,:) = 0.0d0
       DO K = 2,NZM
        DO J = 2,NYM
         RHS_TMP(K,J) =-dimag(CR1X(I,K,J))
!c          RHS_TMP(K,J) = 1.0 
        ENDDO
       ENDDO
       P_TMP(:,:) = 0.0d0



       CALL MULTIGRID(P_TMP,RHS_TMP,I)


!C     Solve for real part of phi
       RHS_TMP(:,:) = 0.0d0
       DO K = 2,NZM
        DO J = 2,NYM
         RHS_TMP(K,J) = -dble(CR1X(I,K,J))
!c          RHS_TMP(K,J) = 1.0
        ENDDO
       ENDDO

       P_TMP2(:,:) = 0.0d0
 
       CALL MULTIGRID(P_TMP2,RHS_TMP,I)
       
 
       DO J = 1,NY
        DO K = 1,NZ
         CR1X(I,K,J) = CMPLX(P_TMP2(K,J),P_TMP(K,J))
        ENDDO
       ENDDO 
      ENDDO 

!C Now, Solve for CUi, the divergenceless velocity field
      DO J=1,NY
        DO K=1,NZ
          DO I = 0,NX2P
            CU1X(I,K,J)=CU1X(I,K,J)-CIKXP(I)*CR1X(I,K,J)
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=1,NZ
          DO I=0,NX2P
            CU2X(I,K,J)=CU2X(I,K,J)-(CR1X(I,K,J) &
                  -CR1X(I,K,J-1))/DY(J)
          END DO
        END DO
      END DO

      DO J=1,NY
        DO K=2,NZ
          DO I=0,NX2P
            CU3X(I,K,J)=CU3X(I,K,J)-(CR1X(I,K,J)  &
                  -CR1X(I,K-1,J))/DZ(K)
          END DO
        END DO
      END DO

     IF (CALL_CHEEK_DIV) THEN 
      CS1X(:,:,:)  = (0.0d0,0.0d0)
! C  Now, create the RHS vector
      DO J=2,NYM
       DO K=2,NZM
        DO I=0,NX2P
           CS1X(I,K,J)=CIKXP(I)*CU1X(I,K,J)          &
                + (CU2X(I,K,J+1)-CU2X(I,K,J))/DYF(J)  &
                + (CU3X(I,K+1,J)-CU3X(I,K,J))/DZF(K)
        END DO
       END DO
      END DO

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

!      CALL FFT_X_TO_PHYSICAL(CRTMP,RTMP,0,NY+1,0,NZ+1)

      DIV = 0.0D0
      DO J = 2,NYM
       DO K = 2,NZM
        DO I = 1,NXP
         DIV = DIV + DABS(S1X(I,K,J))
        ENDDO
       ENDDO
      ENDDO
      DIV = DIV/(DBLE(NXP*NP)*DBLE(NY)*DBLE(NZ))

      CALL MPI_COMBINE_STATS(DIV,1,1)
!C       CALL FFT_X_TO_FOURIER(RTMP,CRTMP,0,NY+1,0,NZ+1)

      IF (rank .eq. 0) THEN
      write(6,*) "THE DIVERGENCE IS ", DIV,  &
                " AFTER REMOVING THE DIVERGENCE"

      ENDIF
      ENDIF


      

      RETURN
      END


!C----*|r-.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE POISSON_P_DUCT
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C We have CUi, need to compute CP.  Solve using multigrid

use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG, BC, CALL_CHEEK_DIV
use pr_rem
use mpi_var      
implicit none


      

      RETURN
      END



!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE VIS_FLOW_DUCT
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

           
      subroutine courant_duct
! This subroutine sets the timestep based on the specified CFL number
! The subroutine should be called with the velocity in physical space

use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG, BC, CALL_CHEEK_DIV
use pr_rem
use mpi_var      
implicit none
!      include 'header_duct'

      real*8 vel
      real*8 dt
      real*8 dt_x,dt_y,dt_z,min_x,dt_dif
      integer i,j,k
      integer imin,jmin,kmin

! Set the initial dt to some arbitrary large number
      dt=999.d0

!    dt based on Diffusion no 
      min_x = LX/dble(NX)       
      dt_dif    = DFN*min_x**2/NU
     
      do j=1,NY
        do k=0,NZM
          do i=0,NXP
            dt_x=cfl*dx(i)/abs(U1X(i,k,j))
            dt_y=cfl*dy(j)/abs(U2X(i,k,j))
            dt_z=cfl*dz(k)/abs(U3X(i,k,j))
!c            if( dt_z .le. 0 ) then
!c             write(6,*) 'dt=', dt_x, 'i,j,k=',i,j,k,dz(k),abs(U3(i,k,j))
!c            endif       
            dt=min(dt,dt_x,dt_y,dt_z,dt_dif)
          end do
        end do
      end do

      if (dt.le.0) then
        write(*,*) 'Error: dt<=0 in courant'
! Set DELTA_T to some small default value
        write(*,*) dt
        DELTA_T=0.0001d0
      else if (dt.ge.1.) then
        write(*,*) 'WARNING: DELTA_T > 0.1, value capped at 0.1'
        DELTA_T=0.10d0
        H_BAR(1)=DELTA_T*(8.0/15.0)
        H_BAR(2)=DELTA_T*(2.0/15.0)
        H_BAR(3)=DELTA_T*(5.0/15.0)
      else
        DELTA_T=dt
        H_BAR(1)=DELTA_T*(8.0/15.0)
        H_BAR(2)=DELTA_T*(2.0/15.0)
        H_BAR(3)=DELTA_T*(5.0/15.0)
      end if

      return
      end


subroutine courant_duct_mpi
! This subroutine sets the timestep based on the specified CFL number
! The subroutine should be called with the velocity in physical space

use Domain
use Grid
!use Fft_var, only : NKX
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use mpi_var, only : rank      
implicit none

      real*8 vel
      real*8 dt
      real*8 dt_x,dt_y,dt_z,min_x,dt_dif
      integer i,j,k
      integer imin,jmin,kmin

! Set the initial dt to some arbitrary large number
      dt=999.d0

!    dt based on Diffusion no 
      min_x = LX/dble(NX)       
      dt_dif    = DFN*min_x**2/NU
     
      do j=1,NY
        do k=0,NZM
          do i=0,NXP     
            dt_x=cfl*dx(i)/abs(U1X(i,k,j))
            dt_y=cfl*dy(j)/abs(U2X(i,k,j))
            dt_z=cfl*dz(k)/abs(U3X(i,k,j))
! c            if( dt_z .le. 0 ) then
! c             write(6,*) 'dt=', dt_x, 'i,j,k=',i,j,k,dz(k),abs(U3(i,k,j))
! c            endif       
            dt=min(dt,dt_x,dt_y,dt_z)
          end do
        end do
      end do

! c      stop
      CALL MPI_COURANT(dt)

      if (dt.le.0) then
        write(*,*) 'Error: dt<=0 in courant'
! Set DELTA_T to some small default value
        write(*,*) dt
! c        DELTA_T=0.0001d0
      else if (dt.ge.DELTA_T_in) then
       if (rank .eq. 0) then  
        write(*,*) 'WARNING: DELTA_T > ',DELTA_T_in, ', value capped at',DELTA_T_in,',dt=', dt
       endif
        DELTA_T= DELTA_T_in
        H_BAR(1)=DELTA_T*(8.0/15.0)
        H_BAR(2)=DELTA_T*(2.0/15.0)
        H_BAR(3)=DELTA_T*(5.0/15.0)
! c        DELTA_T=0.0001d0
      else
        DELTA_T=dt
        H_BAR(1)=DELTA_T*(8.0/15.0)
        H_BAR(2)=DELTA_T*(2.0/15.0)
        H_BAR(3)=DELTA_T*(5.0/15.0)
      end if

      return
      end


       subroutine sink_momentum

! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state for the velocity field
! The intention is to allow an open boundary
use ntypes
use Domain
use Grid
!use Fft_var, only : NKX
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use mpi_var
implicit none

! The following variables will store the background state
      real*8 TH_0(0:NY+1)
      real*8 TH_TOP,b,co_b
      integer i,j,k,N,NZ_C,NY_C

! Damp fluctuation flowfield

      

      do j=jstart,jend
       do k=zstart,zend
        do i=1,NX2P
        CF1X(i,k,j)=CF1X(i,k,j)-SINK_PROF(k,j)*(CU1X(i,k,j)-0.d0)
        CF2X(i,k,j)=CF2X(i,k,j)-SINK_PROF(k,j)*(CU2X(i,k,j)-0.d0)
        CF3X(i,k,j)=CF3X(i,k,j)-SINK_PROF(k,j)*(CU3X(i,k,j)-0.d0)

!        CF1X(i,k,j)=CF1X(i,k,j)-SPONGE_TEMP(k,j)*(CU1X(i,k,j)-0.d0) &
!                           *INT_JACOB(K,J)
!        CF2X(i,k,j)=CF2X(i,k,j)-SPONGE_TEMP(k,j)*(CU2X(i,k,j)-0.d0) &
!                           *INT_JACOB(K,J)
!        CF3X(i,k,j)=CF3X(i,k,j)-SPONGE_TEMP(k,j)*(CU3X(i,k,j)-0.d0) &
!                           *INT_JACOB(K,J)
        end do
       end do
      end do


! Damp mean flow
      IF (RANK .EQ. 0) THEN

      do j=jstart,jend
       do k=zstart,zend
        CF1X(0,k,j)=CF1X(0,k,j)-SINK_PROF(k,j)*(CU1X(0,k,j)-0.d0)
        CF2X(0,k,j)=CF2X(0,k,j)-SINK_PROF(k,j)*(CU2X(0,k,j)-0.d0)
        CF3X(0,k,j)=CF3X(0,k,j)-SINK_PROF(k,j)*(CU3X(0,k,j)-0.d0)


!        CF1X(0,k,j)=CF1X(0,k,j)-SPONGE_SIGMA(k,j)*(CU1X(0,k,j)-0.d0)  &
!                           *INT_JACOB(K,J)
!        CF2X(0,k,j)=CF2X(0,k,j)-SPONGE_SIGMA(k,j)*(CU2X(0,k,j)-0.d0)  &
!                           *INT_JACOB(K,J)
!        CF3X(0,k,j)=CF3X(0,k,j)-SPONGE_SIGMA(k,j)*(CU3X(0,k,j)-0.d0)  &
!                           *INT_JACOB(K,J)
       end do
      end do

      ELSE
        i = 0

      do j=jstart,jend
       do k=zstart,zend
        CF1X(i,k,j)=CF1X(i,k,j)-SINK_PROF(k,j)*(CU1X(i,k,j)-0.d0)
        CF2X(i,k,j)=CF2X(i,k,j)-SINK_PROF(k,j)*(CU2X(i,k,j)-0.d0)
        CF3X(i,k,j)=CF3X(i,k,j)-SINK_PROF(k,j)*(CU3X(i,k,j)-0.d0)

!        CF1X(i,k,j)=CF1X(i,k,j)-SPONGE_TEMP(k,j)*(CU1X(i,k,j)-0.d0) &
!                           *INT_JACOB(K,J)
!        CF2X(i,k,j)=CF2X(i,k,j)-SPONGE_TEMP(k,j)*(CU2X(i,k,j)-0.d0) &
!                           *INT_JACOB(K,J)
!        CF3X(i,k,j)=CF3X(i,k,j)-SPONGE_TEMP(k,j)*(CU3X(i,k,j)-0.d0) &
!                           *INT_JACOB(K,J)
       end do
      end do

      ENDIF

      return
      end
!------------------------------------------------------------------.----|--------|


!CCC BLANK SUBROUTINE FOR COMPILATION
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_CHAN_1(FINAL)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      return
      end

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_CHAN_2(FINAL)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      return
      end

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_CHAN_3(FINAL)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      return
      end

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_CHAN
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       

      RETURN
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE POISSON_P_CHAN
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE create_flow_chan
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_CHAN

      RETURN
      END
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
        SUBROUTINE filter_chan
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE create_TH_chan
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS_CHAN(FINAL)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      return
      end 
