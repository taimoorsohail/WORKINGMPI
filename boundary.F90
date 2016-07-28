
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_REAL_SP(A,B,C,G,D,K,NY,NX)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
!C The RHS vector and solution are real
!C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
!C Returns solution in x
!C The indexing should be done by ROW, ie.
!C [ b1  c1   0   0   0 ...
!C [ a2  b2  c2   0   0 ...
!C [  0  a3  b3   c3  0 ...

      INTEGER I, J, NX, NY,k
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY), G(0:NX,0:NY)
      REAL*8 D

       J= K
!c       write(6,*) 'I am here', J
       DO I=0,NX
!c          D=-D/A(I,J+1)
          A(I,J+1)=A(I,J+1)-B(I,J)*D/A(I,J)
          B(I,J+1)=B(I,J+1)-C(I,J)*D/A(I,J)
          G(I,J+1)=G(I,J+1)-G(I,J)*D/A(I,J)
       END DO

!c       write(6,*) 'I am here after', J

      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO J=NY-1,0,-1
        DO I=0,NX
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END




!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_REAL(A,B,C,G,NY,NX)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
!C The RHS vector and solution are real
!C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
!C Returns solution in x
!C The indexing should be done by ROW, ie.
!C [ b1  c1   0   0   0 ...
!C [ a2  b2  c2   0   0 ...
!C [  0  a3  b3   c3  0 ...

      INTEGER I, J, NX, NY
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY), G(0:NX,0:NY)

      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO J=NY-1,0,-1
        DO I=0,NX
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|    
      SUBROUTINE THOMAS_COMPLEX(A,B,C,G,NY,NX)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

!C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
!C The RHS vector and solution is complex
!C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
!C Returns solution in x
!C The indexing should be done by ROW, ie.
!C [ b1  c1   0   0   0 ...
!C [ a2  b2  c2   0   0 ...
!C [  0  a3  b3   c3  0 ...

      INTEGER I, J, NY, NX
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY)
      COMPLEX*16 G(0:NX,0:NY)

      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO I=0,NX
        DO J=NY-1,0,-1
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END


!C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_TH_LOWER(N,K)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-----
!      INCLUDE 'header_duct'
 
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

     INTEGER I,N,K,J

!C Bottom Wall:
      IF ( TH_BC_YMIN(N) .EQ.0) THEN
!C Dirichlet
        DO I=0,NXP
          MATL(I,0)=0. 
          MATD(I,0)=1.
          MATU(I,0)=0.                   
          VEC(I,0)=0.

          MATL(I,1)=0. 
          MATD(I,1)=1.
          MATU(I,1)=0.                   
          VEC(I,1)=TH_BC_YMIN_C1(N) 
        END DO
      ELSE IF (TH_BC_YMIN(N) .eq. 1) THEN
!C Neumann
        DO I=0,NXP
          MATL(I,1)=0.
          MATD(I,1)=-1.
          MATU(I,1)=1.
          VEC(I,1)=DY(2)*TH_BC_YMIN_C1(N)
        END DO

        DO I=0,NXP
          MATL(I,0)=0.
          MATD(I,0)=-1.
          MATU(I,0)=1.
          VEC(I,0)=DY(1)*TH_BC_YMIN_C1(N)
        END DO 

        IF(IBM)THEN
         DO I=0,NXP
          DO J=0,hill_ind(I)
           MATL(I,J)=0.
           MATD(I,J)=-1.
           MATU(I,J)=1.
           VEC(I,J)=DY(J+1)*TH_BC_YMIN_C1(N) 
          ENDDO
         ENDDO
        ENDIF

      ElSE IF (TH_BC_YMIN(N) .EQ. 7) THEN
!C Mixed boundary condition for w at the bottom wall 
       
       DO I=0,NXP
          MATL(I,1)=0.
          MATD(I,1)=1.
          MATU(I,1)=0.
          VEC(I,1)= TH_1 + (TH_2-TH_1)*sqrt(0.5*(1.0**2.0- 0.d0**2.0) &
         *tanh((GZF(k)-Loc_s)/Lh_s) + 0.5*(1.0d0**2.0 + 0.0d0**2.0)) ;
          MATL(I,0)=0.
          MATD(I,0)=1.
          MATU(I,0)=0.
          VEC(I,0) = TH_1 + (TH_2-TH_1)*sqrt(0.5*(1.0**2.0- 0.d0**2.0) &
         *tanh((GZF(k)-Loc_s)/Lh_s) + 0.5*(1.0d0**2.0 + 0.0d0**2.0)) ;
        END DO
! 
!       IF ( K .lE. 80 ) THEN
!         DO I=0,NXP
!           MATL(I,1)=0.
!           MATD(I,1)=-1.
!           MATU(I,1)=1.
!           VEC(I,1)=DY(2)*TH_BC_YMIN_C1(N)          
! 
!           MATL(I,0)=0.
!           MATD(I,0)=-1.
!           MATU(I,0)=1.
!           VEC(I,0)=DY(1)*TH_BC_YMIN_C1(N)
!         END DO
!        ELSE
!         DO I=0,NXP
!           MATL(I,0)=0.
!           MATD(I,0)=1.
!           MATU(I,0)=0.
!           VEC(I,0)=0.
! 
!           MATL(I,1)=0.
!           MATD(I,1)=1.
!           MATU(I,1)=0.
!           VEC(I,1)=TH_BC_YMIN_C1(N)
!         END DO   
!        END IF   
      ELSE
        write(*,*) 'WARNING: TH_BC_LOWER is of unknown type'
      END IF

      RETURN 
      END



!C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_TH_UPPER(N,K)
!C----*|--.---------.---------.---------.---------.---------.---------.-|--
 !     INCLUDE 'header_duct'
  
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

    INTEGER I, N,K

!C Top wall
      IF (TH_BC_YMAX(N) .EQ. 0) THEN
!C Dirichlet
        DO I=0,NXP
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=TH_BC_YMAX_C1(N)

          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=TH_BC_YMAX_C1(N)
        END DO
      ELSE IF (TH_BC_YMAX(N) .eq. 1) THEN
!C Neumann
        DO I=0,NXP
          MATL(I,NY)=-1.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=DY(NY)*TH_BC_YMAX_C1(N)
        END DO
        DO I=0,NXP
          MATL(I,NY+1)=-1
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=DY(NY+1)*TH_BC_YMAX_C1(N)
        END DO 
      ElSE IF (TH_BC_YMAX(N) .EQ. 7) THEN
!C Mixed boundary condition for w at the bottom wall 

        DO I=0,NXP
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)= TH_1 + (TH_2-TH_1)*sqrt(0.5*(1.0**2.0- 0.d0**2.0) &
         *tanh((GZF(k)-Loc_s)/Lh_s) + 0.5*(1.0d0**2.0 + 0.0d0**2.0)) ;

          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY) = TH_1 + (TH_2-TH_1)*sqrt(0.5*(1.0**2.0- 0.d0**2.0) &
         *tanh((GZF(k)-Loc_s)/Lh_s) + 0.5*(1.0d0**2.0 + 0.0d0**2.0)) ;
        END DO
      END IF

      RETURN
      END   

!C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_TH_LEFT(N,J)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-----
 !     INCLUDE 'header_duct'
  
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

    INTEGER I,N,J

!C Left Wall:
      IF (TH_BC_ZMIN(N) .EQ.0) THEN
!C Dirichlet
       IF (IC_TYPE .EQ. 5) THEN
        DO I=0,NXP
          MATL_Z(I,0)=0.
          MATD_Z(I,0)=1.
          MATU_Z(I,0)=0.
          VEC_Z(I,0)=0.

          MATL_Z(I,1)=0.
          MATD_Z(I,1)=1.
          MATU_Z(I,1)=0.
          VEC_Z(I,1)=U3_BAR(1,J)
        END DO
       ElSE  
        DO I=0,NXP          
          MATL_Z(I,0)=0. 
          MATD_Z(I,0)=1.
          MATU_Z(I,0)=0.                   
          VEC_Z(I,0)=TH_BC_ZMIN_C1(N)

          MATL_Z(I,1)=0. 
          MATD_Z(I,1)=1.
          MATU_Z(I,1)=0.                   
          VEC_Z(I,1)=TH_BC_ZMIN_C1(N) 
        END DO
       ENDIF
      ELSE IF (TH_BC_ZMIN(N).eq.1) THEN
!C Neumann
        DO I=0,NXP
          MATL_Z(I,1)=0.
          MATD_Z(I,1)=-1.
          MATU_Z(I,1)=1.
          VEC_Z(I,1)=DZ(2)*TH_BC_ZMIN_C1(N)
        END DO

        DO I=0,NXP
          MATL_Z(I,0)=0.
          MATD_Z(I,0)=-1.
          MATU_Z(I,0)=1.
          VEC_Z(I,0)=DZ(1)*TH_BC_ZMIN_C1(N)
        END DO

      ELSE
        write(*,*) 'WARNING: TH_BC_LEFT is of unknown type'
      END IF

      RETURN 
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_TH_RIGHT(N,J)
!C----*|--.---------.---------.---------.---------.---------.---------.-|--
  !    INCLUDE 'header_duct'

use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

   INTEGER I,N,J

!C Right wall
      IF (TH_BC_ZMAX(N).EQ.0) THEN
!C Dirichlet
        DO I=0,NXP
          MATL_Z(I,NZ+1)=0.
          MATD_Z(I,NZ+1)=1.
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1)=TH_BC_ZMAX_C1(N)

          MATL_Z(I,NZ)=0.
          MATD_Z(I,NZ)=1.
          MATU_Z(I,NZ)=0.
          VEC_Z(I,NZ)=TH_BC_ZMAX_C1(N)
        END DO
      ELSE IF (TH_BC_ZMAX(N) .eq. 1) THEN
!C Neumann
        DO I=0,NXP
          MATL_Z(I,NZ)=-1.
          MATD_Z(I,NZ)=1.
          MATU_Z(I,NZ)=0.
          VEC_Z(I,NZ)=DZ(NZ)*TH_BC_ZMAX_C1(N)
        END DO
        DO I=0,NXP
          MATL_Z(I,NZ+1)=-1.
          MATD_Z(I,NZ+1)=1.
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1)=DZ(NZ+1)*TH_BC_ZMAX_C1(N)
        END DO
      ELSE IF (TH_BC_ZMAX(N) .eq. 5) THEN
        DO I=0,NXP
          MATL_Z(I,NZ+1)=-(1.d0/DZ(NZ)+1.d0/DZ(NZ+1))
          MATD_Z(I,NZ+1)=1.d0/DZ(NZ+1)
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1) =0.
        END DO
      END IF

      RETURN
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_1_LOWER(K)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-----
!      INCLUDE 'header_duct'
 
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

     INTEGER I,J,K

!C Bottom Wall:
      IF (U_BC_YMIN.EQ.0) THEN
!C Dirichlet
       If (F_TYPE.EQ.5)THEN
        DO I=0,NXP
          MATL(I,0)=0.
          MATD(I,0)=1.
          MATU(I,0)=0.
          VEC(I,0)=0.

          MATL(I,1)=0.
          MATD(I,1)=1.
          MATU(I,1)=0.
          VEC(I,1)=-Q_H0*(1.d0/cos(ANG_BETA)+GX(NX/2)*tan(ANG_BETA) &
                  *In_H0)*sin(OMEGA0*TIME)           
!c           VEC(I,1)=0.0
        END DO
       ELSE

        DO I=0,NXP
          MATL(I,0)=0. 
          MATD(I,0)=1.
          MATU(I,0)=0.                   
          VEC(I,0)=U_BC_YMIN_C1

          MATL(I,1)=0. 
          MATD(I,1)=1.
          MATU(I,1)=0.                   
          VEC(I,1)=U_BC_YMIN_C1 
        END DO

        IF(IBM)THEN
        DO I=0,NXP
         DO J=0,hill_ind(I)
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=U_BC_YMIN_C1

          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=U_BC_YMIN_C1
         ENDDO
        END DO
        ENDIF

       ENDIF
      ELSE IF (U_BC_YMIN.eq.1) THEN
!C Neumann
        DO I=0,NXP
          MATL(I,0)=0.
          MATD(I,0)=-1.
          MATU(I,0)=1.
          VEC(I,0)=DY(1)*U_BC_YMIN_C1
        END DO
      ELSE IF (U_BC_YMIN.eq.3) THEN
!C Wall model artificial slip boundary condition
!C Neumann
        DO I=0,NXP
          MATL(I,1)=0.
          MATD(I,1)=-1.
          MATU(I,1)=1.
          VEC(I,1)=DY(2)*U_BC_LOWER_NWM(I,K)
        END DO
!C This gridpoint is not used here
        DO I=0,NXP
          MATL(I,0)=0.
          MATD(I,0)=1.
          MATU(I,0)=0.
          VEC(I,0)=0.
        END DO
      ELSE
        write(*,*) 'WARNING: U_BC_LOWER is of unknown type'
      END IF

      RETURN 
      END



!C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_1_UPPER(K)
!C----*|--.---------.---------.---------.---------.---------.---------.-|--
 !     INCLUDE 'header_duct'
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use fft_var, only : PI
use mg_vari, only : INIT_FLAG
      
implicit none


    INTEGER I,K

!C Top wall
      IF (U_BC_YMAX.EQ.0) THEN
!C Dirichlet
        DO I=0,NXP
        IF (GZF(K) <= 0.75*GZF(NZ)) THEN
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=U_BC_YMAX_C1*(DSIN(PI*((4*GZF(K))/(3*(GZF(NZ)-GZF(1))))))**2

          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=U_BC_YMAX_C1*(DSIN(PI*((4*GZF(K))/(3*(GZF(NZ)-GZF(1))))))**2
       
       ELSE
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.

          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=0.
        END IF
         END DO
      ELSE IF (U_BC_YMAX.eq.1) THEN
!C Neumann
        DO I=0,NXP
        IF (GZF(K) <= 0.75*GZF(NZ)) THEN
          MATL(I,NY)=-1.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=DY(NY)*U_BC_YMAX_C1*(DSIN(PI*((4*GZF(K))/(3*(GZF(NZ)-GZF(1))))))**2
         ELSE

          MATL(I,NY)=-1.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=0.
        END IF

        END DO
        DO I=0,NXP
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.
        END DO
      ELSE IF (U_BC_YMAX.EQ.3) THEN
!C Use wall model artificial boundary condition
!C Neumann
        DO I=0,NXP
          MATL(I,NY)=-1.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=DY(NY)*U_BC_UPPER_NWM(I,K)
        END DO
!C This gridpoint is not used here
        DO I=0,NXP
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=U1X(I,K,NY)
        END DO
      END IF

      RETURN
      END   

!C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_1_LEFT(J)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-----
 !     INCLUDE 'header_duct'
  
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

    INTEGER I,J

!C Left Wall:
      IF (U_BC_ZMIN.EQ.0) THEN
!C Dirichlet
        DO I=0,NXP
          MATL_Z(I,0)=0. 
          MATD_Z(I,0)=1.
          MATU_Z(I,0)=0.                   
          VEC_Z(I,0)=0.
        IF (IC_TYPE .EQ. 5) THEN
          MATL_Z(I,1)=0.
          MATD_Z(I,1)=1.
          MATU_Z(I,1)=0.
          VEC_Z(I,1)= U1_BAR(J)
        ELSE  
          MATL_Z(I,1)=0. 
          MATD_Z(I,1)=1.
          MATU_Z(I,1)=0.                   
          VEC_Z(I,1)=U_BC_ZMIN_C1
         ENDIF 
        END DO
      ELSE IF (U_BC_ZMIN.eq.1) THEN
!C Neumann
        DO I=0,NXP
          MATL_Z(I,0)=0.
          MATD_Z(I,0)=-1.
          MATU_Z(I,0)=1.
          VEC_Z(I,0)=DZ(1)*U_BC_ZMIN_C1
        END DO 
      ELSE IF (U_BC_ZMIN.eq.3) THEN
!C Wall model artificial slip boundary condition
!C Neumann
        DO I=0,NXP
          MATL_Z(I,1)=0.
          MATD_Z(I,1)=-1.
          MATU_Z(I,1)=1.
          VEC_Z(I,1)=DZ(2)*U_BC_LOWER_NWM(I,J)
        END DO
!C This gridpoint is not used here
        DO I=0,NXP
          MATL_Z(I,0)=0.
          MATD_Z(I,0)=1.
          MATU_Z(I,0)=0.
          VEC_Z(I,0)=0.
        END DO
      ELSE
        write(*,*) 'WARNING: U_BC_LEFT is of unknown type'
      END IF

      RETURN 
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_1_RIGHT(J)
!C----*|--.---------.---------.---------.---------.---------.---------.-|--
!      INCLUDE 'header_duct'
 
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

     INTEGER I,J

!C Right wall
      IF (U_BC_ZMAX.EQ.0) THEN
!C Dirichlet
        DO I=0,NXP
          MATL_Z(I,NZ+1)=0.
          MATD_Z(I,NZ+1)=1.
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1)=U_BC_ZMAX_C1

          MATL_Z(I,NZ)=0.
          MATD_Z(I,NZ)=1.
          MATU_Z(I,NZ)=0.
          VEC_Z(I,NZ)=U_BC_ZMAX_C1
        END DO
      ELSE IF (U_BC_ZMAX.eq.1) THEN
!C Neumann
        DO I=0,NXP
          MATL_Z(I,NZ)=-1.
          MATD_Z(I,NZ)=1.
          MATU_Z(I,NZ)=0.
          VEC_Z(I,NZ)=DZ(NZ)*U_BC_ZMAX_C1
        END DO
        DO I=0,NXP
          MATL_Z(I,NZ+1)=0.
          MATD_Z(I,NZ+1)=1.
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1)=0.
        END DO
      ELSE IF (U_BC_ZMAX.EQ.3) THEN
!C Use wall model artificial boundary condition
!C Neumann
        DO I=0,NXP
          MATL_Z(I,NZ)=-1.
          MATD_Z(I,NZ)=1.
          MATU_Z(I,NZ)=0.
          VEC_Z(I,NZ)=DZ(NZ)*U_BC_UPPER_NWM(I,J)
        END DO
!C This gridpoint is not used here
        DO I=0,NXP
          MATL_Z(I,NZ+1)=0.
          MATD_Z(I,NZ+1)=1.
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1)=U1X(I,NZ,J)
        END DO
      END IF

      RETURN
      END   

!C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_2_LOWER(K)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-----
   !   INCLUDE 'header_duct'
    
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none


  INTEGER I,J,K

!C Bottom Wall:
      IF (V_BC_YMIN.EQ.0) THEN
!C Dirichlet
       
        DO I=0,NXP
          MATL(I,1)=0. 
          MATD(I,1)=1.
          MATU(I,1)=0.                   
          VEC(I,1)=V_BC_YMIN_C1

          MATL(I,2)=0. 
          MATD(I,2)=1.
          MATU(I,2)=0.                   
          VEC(I,2)=V_BC_YMIN_C1 
        END DO

        IF (IBM) THEN
        DO I=0,NXP
         DO J=1,hill_ind(I)
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=V_BC_YMIN_C1

          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=V_BC_YMIN_C1
         ENDDO
        END DO
        ENDIF

      ELSE IF (V_BC_YMIN.eq.1) THEN
!C Neumann
        DO I=0,NXP
          MATL(I,1)=0.
          MATD(I,1)=-1.
          MATU(I,1)=1.
          VEC(I,1)=DYF(1)*V_BC_YMIN_C1
        END DO
      ELSE IF (V_BC_YMIN.eq.3) THEN
!C Wall model artificial slip boundary condition
!C Neumann
        DO I=0,NXP
          MATL(I,1)=0.
          MATD(I,1)=1.
          MATU(I,1)=0.
          VEC(I,1)=0.
  
          MATL(I,2)=0.
          MATD(I,2)=1.
          MATU(I,2)=0.
          VEC(I,2)=0.
        END DO
!C This gridpoint is not used here
        DO I=0,NXP
          MATL(I,0)=0.
          MATD(I,0)=1.
          MATU(I,0)=0.
          VEC(I,0)=0.
        END DO
      ELSE
        write(*,*) 'WARNING: V_BC_LOWER is of unknown type'
      END IF

!C The following is only a placeholder, this row is used for U1 and U3
      DO I=0,NXP
        MATL(I,0) = 0.
        MATD(I,0) = 1.
        MATU(I,0) = 0.
        VEC(I,0) = 0.
      END DO



      RETURN 
      END



!C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_2_UPPER(K)
!C----*|--.---------.---------.---------.---------.---------.---------.-|--
!      INCLUDE 'header_duct'
 use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none


     INTEGER I,K

!C Top wall
      IF (V_BC_YMAX.EQ.0) THEN
!C Dirichlet
        DO I=0,NXP
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=V_BC_YMAX_C1

          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=V_BC_YMAX_C1
        END DO
      ELSE IF (V_BC_YMAX.eq.1) THEN
!C Neumann
        DO I=0,NXP
          MATL(I,NY+1)=-1.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=DYF(NY)*V_BC_YMAX_C1
        END DO
      ELSE IF (V_BC_YMAX.EQ.3) THEN
!C Use wall model artificial boundary condition
!C Neumann
        DO I=0,NXP
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.
        END DO
!C This gridpoint is not used here
        DO I=0,NXP
          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=0.
        END DO
      END IF

      RETURN
      END   

!C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_2_LEFT(J)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-----
 !     INCLUDE 'header_duct'
  
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

    INTEGER I,J

!C Left Wall:
      IF (V_BC_ZMIN.EQ.0) THEN
!C Dirichlet
        DO I=0,NXP
          MATL_Z(I,0)=0. 
          MATD_Z(I,0)=1.
          MATU_Z(I,0)=0.                   
          VEC_Z(I,0)=V_BC_ZMIN_C1
         IF (IC_TYPE .EQ. 5) THEN
          MATL_Z(I,1)=0.
          MATD_Z(I,1)=1.
          MATU_Z(I,1)=0.
          VEC_Z(I,1)=U2_BAR(1,j) 
         ELSE
          MATL_Z(I,1)=0. 
          MATD_Z(I,1)=1.
          MATU_Z(I,1)=0.                   
          VEC_Z(I,1)=V_BC_ZMIN_C1
         ENDIF 
        END DO
      ELSE IF (V_BC_ZMIN.eq.1) THEN
!C Neumann
        DO I=0,NXP
          MATL_Z(I,0)=0.
          MATD_Z(I,0)=-1.
          MATU_Z(I,0)=1.
          VEC_Z(I,0)=DZ(1)*V_BC_ZMIN_C1
        END DO

      ELSE IF (V_BC_ZMIN.eq.3) THEN
!C Wall model artificial slip boundary condition
!C Neumann
!c        DO I=0,NXP
!c          MATL_Z(I,1)=0.
!c          MATD_Z(I,1)=-1.
!c          MATU_Z(I,1)=1.
!c          VEC_Z(I,1)=0.
!c        END DO
!C This gridpoint is not used here
        DO I=0,NXP
          MATL_Z(I,0)=0.
          MATD_Z(I,0)=1.
          MATU_Z(I,0)=0.
          VEC_Z(I,0)=DZ(1)*V_BC_ZMIN_C1
        END DO
      ELSE
        write(*,*) 'WARNING: V_BC_LEFT is of unknown type'
      END IF

      RETURN 
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_2_RIGHT(J)
!C----*|--.---------.---------.---------.---------.---------.---------.-|--
 !     INCLUDE 'header_duct'
  
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

    INTEGER I,J

!C Right wall
      IF (V_BC_ZMAX.EQ.0) THEN
!C Dirichlet
        DO I=0,NXP
          MATL_Z(I,NZ+1)=0.
          MATD_Z(I,NZ+1)=1.
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1)=V_BC_ZMAX_C1

          MATL_Z(I,NZ)=0.
          MATD_Z(I,NZ)=1.
          MATU_Z(I,NZ)=0.
          VEC_Z(I,NZ)=V_BC_ZMAX_C1
        END DO
      ELSE IF (V_BC_ZMAX.eq.1) THEN
!C Neumann
        DO I=0,NXP
          MATL_Z(I,NZ)=-1.
          MATD_Z(I,NZ)=1.
          MATU_Z(I,NZ)=0.
          VEC_Z(I,NZ)=DZ(NZ)*V_BC_ZMAX_C1
        END DO
        DO I=0,NXP
          MATL_Z(I,NZ+1)=0.
          MATD_Z(I,NZ+1)=1.
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1)=0.
        END DO
      ELSE IF (V_BC_ZMAX.EQ.4) THEN
!C  outlet Boundary condition
        DO I=0,NXP
          MATL_Z(I,NZ+1)=0.
          MATD_Z(I,NZ+1)=1.
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1)=U2X(I,NZ-1,J)+(U2X(I,NZ,J)-U2X(I,NZ-1,J)) &
                   *DZ(NZ+1)/DZ(NZ)

          MATL_Z(I,NZ)=-1.d0/DZ(NZ)
          MATD_Z(I,NZ)=1.d0/DZ(NZ)+1.d0/DZ(NZ+1)
          MATU_Z(I,NZ)=-1.d0/DZ(NZ+1)
          VEC_Z(I,NZ) = 0.
        END DO

!c        DO I=0,NXP
!c          MATL_Z(I,NZ+1)=0.
!c          MATD_Z(I,NZ+1)=1.
!c          MATU_Z(I,NZ+1)=0.
!c          VEC_Z(I,NZ+1)=U2(I,NZ-1,J)+(U2(I,NZ,J)-U2(I,NZ-1,J))
!c     +              *DZ(NZ+1)/DZ(NZ)

!c          MATL_Z(I,NZ)=-(1.d0/DZ(NZ-1)+1.d0/DZ(NZ))
!c          MATD_Z(I,NZ)=1.d0/DZ(NZ)
!c          MATU_Z(I,NZ)=0.
!c          VEC_Z(I,NZ) =0.
!c        END DO
      ELSE IF ( V_BC_ZMAX.EQ.5 ) THEN

!c        DO I=0,NXP
!c          MATL_Z(I,NZ)=-(1.d0/DZ(NZ-1)+1.d0/DZ(NZ))
!c          MATD_Z(I,NZ)=1.d0/DZ(NZ)
!c          MATU_Z(I,NZ)=0.
!c          VEC_Z(I,NZ) =0.

!c          MATL_Z(I,NZ+1)=0.
!c          MATD_Z(I,NZ+1)=1.
!c          MATU_Z(I,NZ+1)=0.
!c          VEC_Z(I,NZ+1) =0.
!c        END DO

        DO I=0,NXP
!c         MATL_Z(I,NZ)=-(1.d0/DZ(NZ-1)+1.d0/DZ(NZ))
!c          MATD_Z(I,NZ)=1.d0/DZ(NZ)
!c          MATU_Z(I,NZ)=0.
!c          VEC_Z(I,NZ) =0.

          MATL_Z(I,NZ+1)=-(1.d0/DZ(NZ)+1.d0/DZ(NZ+1))
          MATD_Z(I,NZ+1)=1.d0/DZ(NZ+1)
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1) =0.
        END DO

      ELSE IF (V_BC_ZMAX.EQ.3) THEN
!C Use wall model artificial boundary condition
!C Neumann
        DO I=0,NXP
          MATL_Z(I,NZ)=-1.
          MATD_Z(I,NZ)=1.
          MATU_Z(I,NZ)=0.
          VEC_Z(I,NZ)=0.0
        END DO
!C This gridpoint is not used here
        DO I=0,NXP
          MATL_Z(I,NZ+1)=0.
          MATD_Z(I,NZ+1)=1.
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1)=U2X(I,NZ,J)
        END DO
      END IF

      RETURN
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_3_LOWER(K)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-----
!      INCLUDE 'header_duct'
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none


      INTEGER I,J,K

!C Bottom Wall:
      IF (W_BC_YMIN.EQ.0) THEN
!C Dirichlet
        DO I=0,NXP
          MATL(I,0)=0. 
          MATD(I,0)=1.
          MATU(I,0)=0.                   
          VEC(I,0)=W_BC_YMIN_C1

          MATL(I,1)=0. 
          MATD(I,1)=1.
          MATU(I,1)=0.                   
          VEC(I,1)=W_BC_YMIN_C1 
        END DO

        IF (IBM) THEN
        DO I=0,NXP
         DO J=0,hill_ind(I)
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=W_BC_YMIN_C1

          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=W_BC_YMIN_C1
        END DO
        END DO
        END IF

      ELSE IF (W_BC_YMIN.eq.1) THEN
!C Neumann
        DO I=0,NXP
          MATL(I,0)=0.
          MATD(I,0)=-1.
          MATU(I,0)=1.
          VEC(I,0)=DY(1)*W_BC_YMIN_C1
        END DO
      ELSE IF (W_BC_YMIN.eq.3) THEN
!C Wall model artificial slip boundary condition
!C Neumann
        DO I=0,NXP
          MATL(I,1)=0.
          MATD(I,1)=-1.
          MATU(I,1)=1.
          VEC(I,1)=DY(2)*W_BC_LOWER_NWM(I,K)
        END DO
!C This gridpoint is not used here
        DO I=0,NXP
          MATL(I,0)=0.
          MATD(I,0)=1.
          MATU(I,0)=0.
          VEC(I,0)=0.
        END DO
      ElSE IF (W_BC_YMIN .EQ. 7) THEN
!C Mixed boundary condition for w at the bottom wall 
       IF ( K .lE. 80 ) THEN
        DO I=0,NXP
          MATL(I,0)=0.
          MATD(I,0)=-1.
          MATU(I,0)=1.
          VEC(I,0)=DY(1)*W_BC_YMIN_C1
        END DO
       ELSE
        DO I=0,NXP
          MATL(I,0)=0.
          MATD(I,0)=1.
          MATU(I,0)=0.
          VEC(I,0)=0.

          MATL(I,1)=0.
          MATD(I,1)=1.
          MATU(I,1)=0.
          VEC(I,1)=W_BC_YMIN_C1
        END DO
       ENDIF
 
      ELSE
        write(*,*) 'WARNING: W_BC_LOWER is of unknown type'
      END IF

      RETURN 
      END



!C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_3_UPPER(K)
!C----*|--.---------.---------.---------.---------.---------.---------.-|--
 !     INCLUDE 'header_duct'

use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none


    INTEGER I,K

!C Top wall
      IF (W_BC_YMAX.EQ.0) THEN
!C Dirichlet
        DO I=0,NXP
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=W_BC_YMAX_C1

          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=W_BC_YMAX_C1
        END DO
      ELSE IF (W_BC_YMAX.eq.1) THEN
!C Neumann
        DO I=0,NXP
          MATL(I,NY)=-1.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=DY(NY)*W_BC_YMAX_C1
        END DO
        DO I=0,NXP
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.
        END DO
      ELSE IF (W_BC_YMAX.EQ.3) THEN
!C Use wall model artificial boundary condition
!C Neumann
        DO I=0,NXP
          MATL(I,NY)=-1.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=DY(NY)*W_BC_UPPER_NWM(I,K)
        END DO
!C This gridpoint is not used here
        DO I=0,NXP
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=U3X(I,K,NY)
        END DO
      END IF

      RETURN
      END   

!C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_3_LEFT(J)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-----
 !     INCLUDE 'header_duct'
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,J,K

!C Bottom Wall:
      IF (W_BC_ZMIN.EQ.0) THEN
!C Dirichlet
       
        DO I=0,NXP
          IF (IC_TYPE .EQ. 5) THEN
          MATL_Z(I,1)=0.
          MATD_Z(I,1)=1.
          MATU_Z(I,1)=0.
          VEC_Z(I,1)= U3_BAR(1,J)

          MATL_Z(I,2)=0.
          MATD_Z(I,2)=1.
          MATU_Z(I,2)=0.
          VEC_Z(I,2)=U3_BAR(1,J)
        ELSE
          MATL_Z(I,1)=0. 
          MATD_Z(I,1)=1.
          MATU_Z(I,1)=0.                   
          VEC_Z(I,1)=W_BC_ZMIN_C1

          MATL_Z(I,2)=0. 
          MATD_Z(I,2)=1.
          MATU_Z(I,2)=0.                   
          VEC_Z(I,2)=W_BC_ZMIN_C1
         ENDIF 
        END DO
      ELSE IF (W_BC_ZMIN.eq.1) THEN
!C Neumann
        DO I=0,NXP
          MATL_Z(I,1)=0.
          MATD_Z(I,1)=-1.
          MATU_Z(I,1)=1.
          VEC_Z(I,1)=DZF(1)*W_BC_ZMIN_C1
        END DO
      ELSE IF (W_BC_ZMIN.eq.3) THEN
!C Wall model artificial slip boundary condition
!C Neumann
        DO I=0,NXP
          MATL_Z(I,1)=0.
          MATD_Z(I,1)=1.
          MATU_Z(I,1)=0.
          VEC_Z(I,1)=0.
  
          MATL_Z(I,2)=0.
          MATD_Z(I,2)=1.
          MATU_Z(I,2)=0.
          VEC_Z(I,2)=0.
        END DO
!C This gridpoint is not used here
        DO I=0,NXP
          MATL_Z(I,0)=0.
          MATD_Z(I,0)=1.
          MATU_Z(I,0)=0.
          VEC_Z(I,0)=0.
        END DO
      ELSE
        write(*,*) 'WARNING: W_BC_LEFT is of unknown type'
      END IF

!C The following is only a placeholder, this row is used for U1 and U2
      DO I=0,NXP
        MATL_Z(I,0) = 0.
        MATD_Z(I,0) = 1.
        MATU_Z(I,0) = 0.
        VEC_Z(I,0) = 0.
      END DO

      RETURN 
      END



!C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_3_RIGHT(J)
!C----*|--.---------.---------.---------.---------.---------.---------.-|--
!      INCLUDE 'header_duct'


use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none


      INTEGER I,K,J

!C Top wall
      IF (W_BC_ZMAX.EQ.0) THEN
!C Dirichlet
        DO I=0,NXP
          MATL_Z(I,NZ+1)=0.
          MATD_Z(I,NZ+1)=1.
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1)=W_BC_ZMAX_C1

          MATL_Z(I,NZ)=0.
          MATD_Z(I,NZ)=1.
          MATU_Z(I,NZ)=0.
          VEC_Z(I,NZ)=W_BC_ZMAX_C1
        END DO
      ELSE IF (W_BC_ZMAX.eq.1) THEN
!C Neumann
        DO I=0,NXP
          MATL_Z(I,NZ+1)=-1.
          MATD_Z(I,NZ+1)=1.
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1)=DZF(NZ)*W_BC_ZMAX_C1
        END DO
        ELSE IF (W_BC_ZMAX.EQ.4) THEN
!C  outlet Boundary condition
!c        DO I=0,NXP
!c          MATL_Z(I,NZ+1)=0.
!c          MATD_Z(I,NZ+1)=1.
!c          MATU_Z(I,NZ+1)=0.
!c          VEC_Z(I,NZ+1)=0.

!c          MATL_Z(I,NZ)=-1.d0/DZF(NZ-1)
!c          MATD_Z(I,NZ)=1.d0/DZF(NZ)+1.d0/DZF(NZ-1)
!c          MATU_Z(I,NZ)=-1.d0/DZF(NZ)
!c          VEC_Z(I,NZ) = 0.
!c        END DO

         DO I=0,NXP
           MATL_Z(I,NZ+1)=0.
           MATD_Z(I,NZ+1)=1.
           MATU_Z(I,NZ+1)=0.
           VEC_Z(I,NZ+1)=U3X(I,NZ-1,J)+(U3X(I,NZ,J)-U3X(I,NZ-1,J)) &
                  *DZF(NZ)/DZF(NZ-1)

           MATL_Z(I,NZ)=-1.d0/DZF(NZ-1)
           MATD_Z(I,NZ)=1.d0/DZF(NZ)+1.d0/DZF(NZ-1)
           MATU_Z(I,NZ)=-1.d0/DZF(NZ)
           VEC_Z(I,NZ) = 0.
          END DO   

!c        DO I=0,NXP
!c          MATL_Z(I,NZ+1)=0.
!c          MATD_Z(I,NZ+1)=1.
!c          MATU_Z(I,NZ+1)=0.
!c          VEC_Z(I,NZ+1)=U3(I,NZ-1,J)+(U3(I,NZ,J)-U3(I,NZ-1,J))
!c     +              *DZF(NZ)/DZF(NZ-1)

!c          MATL_Z(I,NZ)=-(1.d0/DZF(NZ-1)+1.d0/DZF(NZ-2))
!c          MATD_Z(I,NZ)=1.d0/DZF(NZ-1)
!c          MATU_Z(I,NZ)=0.
!c          VEC_Z(I,NZ) =0.
!c        END DO
       ELSE IF (W_BC_ZMAX.EQ.5) THEN

        DO I=0,NXP
          MATL_Z(I,NZ+1)=-(1.d0/DZF(NZ-1)+1.d0/DZF(NZ))
          MATD_Z(I,NZ+1)=1.d0/DZF(NZ)
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1) =0.

!c          MATL_Z(I,NZ)=-(1.d0/DZF(NZ-1)+1.d0/DZF(NZ-2))
!c          MATD_Z(I,NZ)=1.d0/DZF(NZ-1)
!c          MATU_Z(I,NZ)=0.
!c          VEC_Z(I,NZ) =0.
        END DO

      ELSE IF (W_BC_ZMAX.EQ.3) THEN
!C Use wall model artificial boundary condition
!C Neumann
        DO I=0,NXP
          MATL_Z(I,NZ+1)=-1.
          MATD_Z(I,NZ+1)=1.
          MATU_Z(I,NZ+1)=0.
          VEC_Z(I,NZ+1)=DZF(NZ)*W_BC_UPPER_NWM(I,K)
        END DO
!C This gridpoint is not used here
        DO I=0,NXP
          MATL_Z(I,NZ)=0.
          MATD_Z(I,NZ)=1.
          MATU_Z(I,NZ)=0.
          VEC_Z(I,NZ)=0.
        END DO
      END IF

      RETURN
      END


!C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_VEL_LOWER
!C----*|--.---------.---------.---------.---------.---------.---------.-|--
!C This subroutine is called after initializing the flow
!C It sets the appropriate boundary conditions including ghost cell values
!C  on the velocity field in Fourier space
!      INCLUDE 'header_duct'

use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mpi_var, only : rank
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,K      

!C Now, apply the boundary conditions depending on the type specified 
      IF (U_BC_YMIN.EQ.0) THEN
!C Dirichlet 
!C Start with zero
         DO K=ZSTART-1,ZEND+1
           DO I=0,NX2P
             CU1X(I,K,1)=0.d0
           END DO
!C Now, set only the mean
         IF (RANK .EQ. 0) THEN
          IF (F_TYPE .EQ. 5) THEN
           CU1X(0,K,1)=-Q_H0*(1.d0/cos(ANG_BETA)+GX(NX/2)*tan(ANG_BETA) &
                   *In_H0)*sin(OMEGA0*TIME)
!c          CU1X(0,0,1)=U_BC_YMIN_C1
          ELSE
           CU1X(0,K,1)=U_BC_YMIN_C1
          ENDIF
         ENDIF
!C Ghost cell not used
         CU1X(0,K,0)=0.d0
        ENDDO
      ELSE IF (U_BC_YMIN.EQ.1) THEN
!C Neumann
         DO K=ZSTART-1,ZEND+1
           DO I=0,NX2P
             CU1X(I,K,0)=CU1X(I,K,1)-DY(1)*U_BC_YMIN_C1 
           END DO
         END DO
      ELSE
         STOP 'Error: U_BC_YMIN must be 0, or 1'
      END IF

      IF (W_BC_YMIN.EQ.0) THEN
!C Dirichlet
!C Start with zero
         DO K=ZSTART-1,ZEND+1
           DO I=0,NX2P
             CU3X(I,K,1)=0.d0
           END DO
!C Now, set only the mean
           IF (RANK .EQ. 0) THEN
             CU3X(0,K,1)=W_BC_YMIN_C1
!C Ghost cell not used
             CU3X(0,K,0)=0.d0
           ENDIF
         ENDDO
      ELSE IF (W_BC_YMIN.EQ.1) THEN
!C Neumann
         DO K=ZSTART-1,ZEND+1
           DO I=0,NX2P
             CU3X(I,K,0)=CU3X(I,K,1)-DY(1)*W_BC_YMIN_C1 
           END DO
         END DO
      ELSE IF (W_BC_YMIN.EQ. 7) THEN
!C Mixed boundary condition
       DO K=ZSTART-1,ZEND+1
         IF ( K .LE. 80 ) THEN
           DO I=0,NX2P
             CU3X(I,K,0)=CU3X(I,K,1)-DY(1)*W_BC_YMIN_C1
           END DO
         ELSE

           DO I=0,NX2P
             CU3X(I,K,1)=0.d0
           END DO
!C Now, set only the mean
          IF (RANK .EQ. 0) THEN
           CU3X(0,K,1)=W_BC_YMIN_C1
!C Ghost cell not used
           CU3X(0,K,0)=0.d0
          ENDIF
         ENDIF

       END DO

      ELSE
         STOP 'Error: W_BC_YMIN must be 0, or 1' 
      END IF

      IF (V_BC_YMIN.EQ.0) THEN
!C Dirichlet
!C Set the vertical velocity at GYF(1) (halfway between GY(2) and GY(1))
         DO K=ZSTART-1,ZEND+1
           DO I=0,NX2P
             CU2X(I,K,1)=2.d0*V_BC_YMIN_C1-CU2X(I,K,2) 
           END DO
         END DO
      ELSE IF (V_BC_YMIN.EQ.1) THEN
!C Neumann
         DO K=ZSTART-1,ZEND+1
           DO I=0,NX2P
             CU2X(I,K,1)=CU2X(I,K,2)-DYF(1)*V_BC_YMIN_C1 
           END DO
         END DO
      ELSE IF (V_BC_YMIN.EQ.2) THEN
!C Upstream-travelling wave proposed by Speyer/Kim
!C (initialize as zero)
       DO K=ZSTART-1,ZEND+1
        CU2X(0,K,1)=-CU2X(0,K,2)
       ENDDO
      ELSE
         STOP 'Error: V_BC_YMIN must be 0, 1, or 2'
      END IF

      RETURN
      END

     
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE  APPLY_BC_VEL_UPPER
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C This subroutine is called after initializing the flow
!C It sets the appropriate boundary conditions including ghost cell values
!C  on the velocity field in Fourier space
  !    INCLUDE 'header_duct'

use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mpi_var, only : rank
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,K      

! Now, apply boundary conditions to the top of the domain
      IF (U_BC_YMAX.EQ.0) THEN
!C Dirichlet 
!C Start with zero
         DO K=ZSTART-1,ZEND+1
           DO I=0,NX2P
             CU1X(I,K,NY)=0.d0
           END DO
           IF (RANK .EQ. 0) THEN
!C Now, set only the mean
             CU1X(0,K,NY)=U_BC_YMAX_C1
!C Ghost cell not used
             CU1X(0,K,NY+1)=0.d0
           ENDIF
        ENDDO
      ELSE IF (U_BC_YMAX.EQ.1) THEN
!C Neumann
         DO K=ZSTART-1,ZEND+1
           DO I=0,NX2P
             CU1X(I,K,NY+1)=CU1X(I,K,NY)+DY(NY+1)*U_BC_YMAX_C1
           END DO
         END DO
      ELSE
         STOP 'Error: U_BC_YMAX must be 0, or 1'
      END IF

      IF (W_BC_YMAX.EQ.0) THEN
!C Dirichlet
!C Start with zero
         DO K=ZSTART-1,ZEND+1
           DO I=0,NX2P
             CU3X(I,K,NY)=0.d0
           END DO
           IF (RANK .EQ. 0) THEN
!C Now, set only the mean
             CU3X(0,K,NY)=W_BC_YMAX_C1
!C Ghost cell not used
             CU3X(0,K,NY+1)=0.d0
           ENDIF
        ENDDO 
      ELSE IF (W_BC_YMAX.EQ.1) THEN
!C Neumann
         DO K=ZSTART-1,ZEND+1 
           DO I=0,NX2P
             CU3X(I,K,NY+1)=CU3X(I,K,NY)+DY(NY+1)*W_BC_YMAX_C1
           END DO
         END DO
      ELSE
        STOP 'Error: W_BC_YMAX must be 0, or 1'
      END IF

      IF (V_BC_YMAX.EQ.0) THEN
!C Dirichlet
!C Set the vertical velocity at GYF(NY) (halfway between GY(NY) and GY(NY+1))
         DO K=ZSTART-1,ZEND+1
           DO I=0,NX2P
             CU2X(I,K,NY+1)=2.d0*V_BC_YMAX_C1-CU2X(I,K,NY)
           END DO
         END DO
      ELSE IF (V_BC_YMAX.EQ.1) THEN
!C Neumann
         DO K=ZSTART-1,ZEND+1
           DO I=0,NX2P
             CU2X(I,K,NY+1)=CU2X(I,K,NY)+DYF(NY)*V_BC_YMAX_C1
           END DO
         END DO
      ELSE IF (V_BC_YMAX.EQ.2) THEN
!C Upstream-travelling wave proposed by Speyer/Kim
!C (initialize as zero gradient)
        CU2X(0,0,NY+1)=-CU2X(0,0,NY)
      ELSE IF (V_BC_YMAX.EQ.4) THEN
        DO K=ZSTART-1,ZEND+1
           DO I=0,NX2P
            CU2X(I,K,NY+1)= 2.d0*AMP_OMEGA0*sin(OMEGA0*TIME)*  &
                       sin(ANG_BETA) -CU2X(I,K,NY)    
           END DO
        END DO
      ELSE
         STOP 'Error: V_BC_YMAX must be 0, 1, or 2'
      END IF
   
      RETURN
      END

!C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_VEL_LEFT
!C----*|--.---------.---------.---------.---------.---------.---------.-|--
!C This subroutine is called after initializing the flow
!C It sets the appropriate boundary conditions including ghost cell values
!C  on the velocity field in Fourier space
!      INCLUDE 'header_duct'

use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mpi_var, only : rank
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,J      

!C Now, apply the boundary conditions depending on the type specified 
      IF (U_BC_ZMIN.EQ.0) THEN
!C Dirichlet 
!C Start with zero
         DO J=JSTART-1,JEND+1
          DO I=0,NX2P
             CU1X(I,1,J)=0.d0
          END DO
          IF (RANK .EQ. 0) THEN
!C Now, set only the mean
           IF (F_TYPE .EQ. 5) THEN
            CU1X(0,1,J)=-Q_H0*(1.d0/cos(ANG_BETA)+GX(NX/2)*tan(ANG_BETA) &
                   *In_H0)*sin(OMEGA0*TIME)
!c          CU1X(0,1,J)=U_BC_ZMIN_C1
           ELSE
            IF ( IC_TYPE .EQ. 5) THEN
             CU1X(0,1,J)=U1_BAR(J)
            ELSE 
             CU1X(0,1,J)=U_BC_ZMIN_C1
            ENDIF
           ENDIF
!C Ghost cell not used
           CU1X(0,0,J)=0.d0
          ENDIF
        ENDDO
      ELSE IF (U_BC_ZMIN.EQ.1) THEN
!C Neumann
         DO J=JSTART-1,JEND+1
           DO I=0,NX2P
             CU1X(I,0,J)=CU1X(I,1,J)-DZ(1)*U_BC_ZMIN_C1 
           END DO
         END DO
      ELSE
         STOP 'Error: U_BC_ZMIN must be 0, or 1'
      END IF

      IF (W_BC_ZMIN.EQ.0) THEN
!C Dirichlet
!C Start with zero
!c         DO J=JSTART-1,JEND+1
!c           DO I=0,NX2P
!c             CU3X(I,1,J)=0.d0
!c           END DO
!C Now, set only the mean
!c         CU3X(0,1,J)=W_BC_ZMIN_C1
!C Ghost cell not used
!c         CU3X(0,0,J)=0.d0
!c         ENDDO
          DO J=JSTART-1,JEND+1
           DO I=0,NX2P
             CU3X(I,1,J)=0.d0 
           ENDDO
           IF (RANK .EQ. 0) THEN
!C Now, set only the mean
            IF ( IC_TYPE .EQ. 5) THEN           
             CU3X(0,1,J)=2.d0*U3_BAR(1,J)-CU3X(0,2,J)
            ELSE
             CU3X(0,1,J)=2.d0*W_BC_ZMIN_C1-CU3X(0,2,J)
            ENDIF 
           ENDIF 
          END DO
          
      ELSE IF (W_BC_ZMIN.EQ.1) THEN
!C Neumann
         DO J=JSTART-1,JEND+1
           DO I=0,NX2P
             CU3X(I,0,J)=CU3X(I,1,J)-DZ(1)*W_BC_ZMIN_C1 
           END DO
         END DO
      ELSE
         STOP 'Error: W_BC_ZMIN must be 0, or 1' 
      END IF

      IF (V_BC_ZMIN.EQ.0) THEN
!C Dirichlet
!C Start with zero
         DO J=JSTART-1,JEND+1
           DO I=0,NX2P
             CU2X(I,1,J)=0.d0
           END DO
           IF (RANK .EQ. 0) THEN
!C Now, set only the mean
            IF ( IC_TYPE .EQ. 5) THEN
             CU2X(0,1,J)=U2_BAR(1,J)
            ELSE
             CU2X(0,1,J)=V_BC_ZMIN_C1
            ENDIF
!C Ghost cell not used
            CU2X(0,0,J)=0.d0
           ENDIF
         END DO
      ELSE IF (V_BC_ZMIN.EQ.1) THEN
!C Neumann
         DO J=JSTART-1,JEND+1
           DO I=0,NX2P
             CU2X(I,1,J)=CU2X(I,2,J)-DZF(1)*V_BC_ZMIN_C1 
           END DO
         END DO
      ELSE IF (V_BC_ZMIN.EQ.2) THEN
!C Upstream-travelling wave proposed by Speyer/Kim
!C (initialize as zero)
       DO J=JSTART-1,JEND+1
        CU2X(0,1,J)=-CU2X(0,2,J)
       ENDDO
      ELSE
         STOP 'Error: V_BC_ZMIN must be 0, 1, or 2'
      END IF

      RETURN
      END

    
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE  APPLY_BC_VEL_RIGHT
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!C This subroutine is called after initializing the flow
!C It sets the appropriate boundary conditions including ghost cell values
!C  on the velocity field in Fourier space
!      INCLUDE 'header_duct'
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mpi_var, only : rank
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,J      

! Now, apply boundary conditions to the top of the domain
      IF (U_BC_ZMAX.EQ.0) THEN
!C Dirichlet 
!C Start with zero
         DO J=JSTART-1,JEND+1
           DO I=0,NX2P
             CU1X(I,NZ,J)=0.d0
           END DO
           IF (RANK .EQ. 0) THEN 
!C Now, set only the mean
            CU1X(0,NZ,J)=U_BC_ZMAX_C1
!C Ghost cell not used
            CU1X(0,NZ+1,J)=0.d0
           ENDIF
        ENDDO
      ELSE IF (U_BC_ZMAX.EQ.1) THEN
!C Neumann
         DO J=JSTART-1,JEND+1
           DO I=0,NX2P
             CU1X(I,NZ+1,J)=CU1X(I,NZ,J)+DZ(NZ+1)*U_BC_ZMAX_C1
           END DO
         END DO
      ELSE
         STOP 'Error: U_BC_ZMAX must be 0, or 1'
      END IF

      IF (W_BC_ZMAX.EQ.0) THEN
!C Dirichlet
!C Start with zero
!c         DO J=JSTART-1,JEND+1
!c           DO I=0,NX2P
!c             CU3X(I,NZ,J)=0.d0
!c           END DO
!C Now, set only the mean
!c         CU3X(0,NZ,J)=W_BC_ZMAX_C1
!C Ghost cell not used
!c         CU3X(0,NZ+1,J)=0.d0
!c        ENDDO 
         DO J=JSTART-1,JEND+1
           DO I=0,NX2P
             CU3X(I,NZ+1,J)=2.d0*W_BC_ZMAX_C1-CU3X(I,NZ,J)
           END DO
         END DO
      ELSE IF (W_BC_ZMAX.EQ.1) THEN
!C Neumann
         DO J=JSTART-1,JEND+1 
           DO I=0,NX2P
             CU3X(I,NZ+1,J)=CU3X(I,NZ,J)+DZ(NZ+1)*W_BC_ZMAX_C1
           END DO
         END DO
      ELSE IF (W_BC_ZMAX.EQ.4) THEN
!C     Out flow boundary 
       DO J=JSTART-1,JEND+1
           DO I=0,NX2P
             CU3X(I,NZ+1,J)=(CU3X(I,NZ,J)*(1.0/DZF(NZ)+1.0/DZF(NZ-1)) &
                   - CU3X(I,NZ-1,J)/DZF(NZ-1))*DZF(NZ)
           END DO
         END DO
      ELSE IF ( ( W_BC_ZMAX.EQ.5 ) .or. ( W_BC_ZMAX.EQ.6 ) )THEN
!C     Out flow boundary another approach 
       DO J=JSTART,JEND
           DO I=0,NX2P
             CU3X(I,NZ+1,J)=(CU3X(I,NZ,J)*(1.0/DZF(NZ)+1.0/DZF(NZ-1)) &
                   - CU3X(I,NZ-1,J)/DZF(NZ-1))*DZF(NZ)
           END DO
        END DO
 
      ELSE
        STOP 'Error: W_BC_ZMAX must be 0, or 1'
      END IF

      IF (V_BC_ZMAX.EQ.0) THEN
!C Dirichlet
!C Dirichlet
!C Start with zero
         DO J=JSTART-1,JEND+1
           DO I=0,NX2P
             CU2X(I,NZ,J)=0.d0
           END DO
           IF (RANK .EQ. 0) THEN
!C Now, set only the mean
             CU2X(0,NZ,J)=V_BC_ZMAX_C1
!C Ghost cell not used
             CU2X(0,NZ+1,J)=0.d0
           ENDIF
        ENDDO
      ELSE IF (V_BC_ZMAX.EQ.1) THEN
!C Neumann
         DO J=JSTART-1,JEND+1
           DO I=0,NX2P
             CU2X(I,NZ+1,J)=CU2X(I,NZ,J)+DZF(NZ)*V_BC_ZMAX_C1
           END DO
         END DO
      ELSE IF (V_BC_ZMAX.EQ.2) THEN
!C Upstream-travelling wave proposed by Speyer/Kim
!C (initialize as zero gradient)
       DO J=JSTART-1,JEND+1
        CU2X(0,NZ+1,J)=-CU2X(0,NZ,J)
       ENDDO
      ELSE IF (V_BC_ZMAX.EQ.4) THEN
!C     Out flow boundary 
       DO J=JSTART-1,JEND+1
           DO I=0,NX2P
             CU2X(I,NZ+1,J)=(CU2X(I,NZ,J)*(1.0/DZ(NZ)+1.0/DZ(NZ+1)) &
                   - CU2X(I,NZ-1,J)/DZ(NZ))*DZ(NZ+1)
           END DO
         END DO  

      ELSE IF ( (V_BC_ZMAX.EQ.5) .OR. (V_BC_ZMAX.EQ.6)) THEN
!C     Out flow boundary another approach
        DO J=2,NY
           DO I=0,NX2P
             CU2X(I,NZ+1,J)=(CU2X(I,NZ,J)*(1.0/DZ(NZ)+1.0/DZ(NZ+1)) &
                   - CU2X(I,NZ-1,J)/DZ(NZ))*DZ(NZ+1)
           END DO
         END DO

      ELSE
         STOP 'Error: V_BC_ZMAX must be 0, 1, or 2'
      END IF
   
      RETURN
      END
