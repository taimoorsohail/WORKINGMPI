module ntypes
 integer, parameter :: r4=4
 integer, parameter :: r8=8
 integer, parameter :: i4=4
end module ntypes

!DOMAIN
                        module Domain
                        use ntypes

 integer(i4)       :: nx, NY_TOT, NZ, ni, nj, nk, N_TH, NXM, NYM, NZM, TNKZ, TNKY
 integer(i4)       :: NKX, s11, NZ_T
 integer(i4)       :: NP,NXP,NZP,NXV,NZV, NKXV, NKXP,NX2V,NX2P,NXP_L,NX2P_L
 integer(i4)       :: n_pc,T_S_bin,NY_T,T_S_end
  
 !INCLUDE   'grid_def'

 parameter( n_pc = 24)
 parameter(T_S_bin=10672, T_S_end=10764)
 integer(i4) :: np1,np2,np3

 REAL(r8) ::   LX,LZ,LY,DX_s
 REAL(r8) ::   sum_phi_i, phi_d_sum, sum_rho_h, sum_rho_top, phi_d_check, &
               sum_phi_i_avg, temp_val, sum_phi_d_check,&
               sum_phi_d_avg, sum_phi_b2_avg, sum_phi_b2
!,DZF(0:NZ+1),DZ(0:NZ+1), GZ(0:NZ+1), GZF(0:NZ+1), GY(0:NY_TOT+1), GYF(0:NY_TOT+1)
 real(r8), allocatable, dimension(:) :: DZF, DZ, GZ, GZF, GY, GYF, DY, DYF
!PARAMETER (LX = 2.5d0, LZ=0.5d0, LY =0.2d0)

 real(r8) ::  u_0,Gravity,alpha_T, NU, kappa, PR, rho_0
 parameter (u_0=0.125, Gravity = 10.0d0,  alpha_T = 2.0d0*10.0d0**(-4.0), NU = 10.0**(-6), &
                 rho_0 = 1000.0d0, PR = 5.0d0, kappa = NU/PR)
 


 Real(r8) ::     TIME,sum_1,sum_2,dy_top, dy_bot,dx, time_sum, time_old, delta_time
 real(r8) ::     sum_Prod, sum_dissp, sum_buoyF, sum_val,    &
                 sum_dtkdt, sum_vis, sum_trans, sum_advec, dv,dv_z, vol_box
 LOGICAL         LES,xz_plane_read,write_varaview

 real(r8) ::     th_b,sum_dissp_sgs,temp_val_sgs
 real(r8),allocatable,dimension(:)      ::  H_hill_tot, DX_hill
 integer, allocatable,dimension(:)      ::  hill_ind_tot
! real buffer(0:NX-1,0:NZ-1,1:NY)
! The sizes of these strings may need to be changed
! character*11 size_str
! character*29 len_str
!! parameter     (n_pc = 32,  NY =20, NX=128, NZ=32,
!!     &                NY_TOT = NY*n_pc-(n_pc-1),
!!     &                TIME_STEP=47600,
!!     &                T_S_bin=1000)

!! REAL(r8) ::        U1(1:NZ,1:NY_TOT), U2(1:NZ,1:NY_TOT), &
!!                    U3(1:NZ,1:NY_TOT), U4(1:NZ,1:NY_TOT), &
!!                    R1(1:NX,1:NY_TOT), R2(1:NX,1:NY_TOT), &
!!                    R3(1:NX,1:NY_TOT), TH_d(1:NX,1:NY_TOT)
!
!
!
 real(r8),allocatable,dimension(:,:,:) :: U1, U2, U3, TH_WH,TH, NU_T, pr_3D
 real(r8),allocatable,dimension(:,:,:) :: U1_avg, U2_avg, U3_avg, TH_WH_avg
! real(r8),allocatable,dimension(:)     :: u1_avg_ln,u2_avg_ln,u3_avg_ln,TH_WH_avg_ln
! real(r8),allocatable,dimension(:)     :: C_DYN, C_DYN_avg

!! real(r4),allocatable,dimension(:,:,:)   :: buffer
 real(r8),allocatable,dimension(:,:)     :: z_start_bar_surf
!
 real(r8),allocatable,dimension(:)     :: prod_prof,dissip_prof,dissip_sgs_prof,buoy_prof,NU_T_prof,phi_d_prof,phi_b2_prof,g1vtk,g2vtk,g3vtk, &
                                         phi_i,z_start_bar,phi_d 

! REAL(r8) ::       TH1(1:NX,1:NZ)
! real(r4) ::       g1vtk(1:np1),g2vtk(1:np2),g3vtk(1:np3),var_1(1:3*np1*np2)
! !                  th_div_yz(1:np2,1:np3)
!
! CHARACTER*52 OutFileName, outputDIR, basename
!! CHARACTER*42 filn_th_xz_1,filn_th_xz_2,filn_th_xz_3
!
!
end module Domain

module qsort_c_module

implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A,B)
use ntypes
  real(r8), intent(in out), dimension(:) :: A
  real(r8), intent(in out), dimension(:) :: B
  integer :: iq

  if(size(A) > 1) then
     call Partition(A,B,iq)
     call QsortC(A(:iq-1),B(:iq-1))
     call QsortC(A(iq:),B(iq:))
  endif
end subroutine QsortC

subroutine Partition(A,B,marker)
use ntypes
  real(r8), intent(in out), dimension(:) :: A
  real(r8), intent(in out), dimension(:) :: B
  integer, intent(out) :: marker
  integer :: i, j
  real(r8) :: temp
  real(r8) :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp

        temp = B(j)
        B(j) = B(i)
        B(i) = temp

     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module qsort_c_module



Program add_data
use Domain
use ntypes

implicit none
integer ::      i,j,m,k,NY_min,NY_max,COUNTER, &
                T_S    

    character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename

       outputDIR='plane_data_3D/'
       basename = 'data_binary_3D'

       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR), &
       trim(basename)//"_n",T_S_bin,".pln"

       open(unit=22,file=OutFileName,access='stream',  &
       form='unformatted',status='old',             &
       convert='big_endian',iostat=s11)

        write(6,*)'Reading flow statistics', T_S_bin,' ',OutFileName
        read(22) TIME, ni, nk, nj, DX_s, LX, LY, LZ
        close(22)

        NX = ni
        NZ = nk
        NY_TOT = nj

!        allocate (g1vtk(1:ni), g2vtk(1:NY_TOT), g3vtk(1:nk))
       
         call read_grid

        COUNTER = 0;

        allocate (U1(1:ni,1:nk,1:nj), U2(1:ni,1:nk,1:nj),     &
                  U3(1:ni,1:nk,1:nj), TH_WH(1:ni,1:nk,1:nj),&
                  TH(1:ni,1:nk,1:nj), &
                  pr_3D(1:ni,1:nk,1:nj))
        allocate (NU_T(1:ni,1:NZ,1:NY_TOT)) 

        allocate (U1_avg(1:ni,1:NZ,1:NY_TOT), U2_avg(1:ni,1:NZ,1:NY_TOT),     &
                  U3_avg(1:ni,1:NZ,1:NY_TOT), TH_WH_avg(1:ni,1:NZ,1:NY_TOT))
      
        allocate (z_start_bar_surf(1:ni,1:NZ))

        allocate (prod_prof(1:NY_TOT),buoy_prof(1:NY_TOT),dissip_prof(1:NY_TOT),NU_T_prof(1:NY_TOT), &
                  dissip_sgs_prof(1:NY_TOT), phi_d_prof(1:NY_TOT),phi_b2_prof(1:NY_TOT),z_start_bar(1:NY_TOT))        

        allocate(H_hill_tot(ni), hill_ind_tot(ni), DX_hill(ni))
 
        time_sum = 0.0d0 ;
        th_wh_avg(:,:,:)= 0.0d0
        u1_avg(:,:,:)   = 0.0d0
        u2_avg(:,:,:)   = 0.0d0
        u3_avg(:,:,:)   = 0.0d0
        U1(:,:,:) = 0.0d0
        U2(:,:,:) = 0.0d0
        U3(:,:,:) = 0.0d0
        TH_WH(:,:,:) = 0.0d0
        TH(:,:,:) = 0.0d0
        pr_3D(:,:,:) = 0.0d0

        z_start_bar(:)  = 0.0d0 
        z_start_bar_surf(:,:) = 0.0d0

        phi_d_prof(:)   = 0.0d0
        phi_b2_prof(:)  = 0.0d0

        sum_phi_i_avg = 0.0d0
        sum_phi_d_avg = 0.0d0
        sum_phi_b2_avg = 0.0d0
  
        DO T_S=T_S_bin,T_S_end

        outputDIR='plane_data_3D/'
        basename = 'data_binary_3D'

        write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR), &
        trim(basename)//"_n",T_S,".pln"

        open(unit=22,file=OutFileName,access='stream',  &
        form='unformatted',status='old',             &
        convert='big_endian',iostat=s11)

        write(6,*)'Reading flow statistics', T_S,' ',OutFileName
        read(22) TIME, ni, nk, nj, DX_s, LX, LY, LZ
        read(22) (((TH(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        read(22) (((U3(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        read(22) (((U2(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        read(22) (((U1(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        read(22) (((pr_3D(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        close(22)

        th_b = 30.0d0

        do j=1,nj
        do k=1,nk
        do i=1,ni
        TH_WH(i,k,j) = TH(i,k,j)! + th_b
        end do
        end do
        end do
        
        open(unit=43,file='th_test.txt',form='formatted',status='unknown')
        do k=1,NZ
        write(43,111) GZ(K), TH(ni/2, k, NY_TOT), TH(ni, k, NY_TOT)
        end do

       COUNTER = COUNTER + 1
         open(unit=40,file='PE_statistics.txt',form='formatted',status='unknown')!,&
               ! position='append')
       IF (T_S > T_S_bin )THEN 
        delta_time = time-time_old
        time_sum = time_sum + delta_time

        call PHI_cal

        do j=1,NY_TOT
        do k=1,NZ
        do i=1,ni
         th_wh_avg(i,k,j)= th_wh_avg(i,k,j)+th_wh(i,k,j)*delta_time
         u1_avg(i,k,j)   = u1_avg(i,k,j) + u1(i,k,j)*delta_time
         u2_avg(i,k,j)   = u2_avg(i,k,j) + u2(i,k,j)*delta_time
         u3_avg(i,k,j)   = u3_avg(i,k,j) + u3(i,k,j)*delta_time
        end do
        end do
        end do
        END IF 
       time_old = time

       IF (delta_time>0) THEN
        WRITE(6,*) T_s, 'Saving PE Statistics', COUNTER
        write(40,666) COUNTER, phi_d_check, sum_phi_b2, phi_d_sum, &
                  sum_phi_i
!        write(6,*) 'COUNT =', COUNTER, 'phi_b2 =', sum_phi_b2, 'phi_d =', phi_d_sum, &
!                   'phi_i=', sum_phi_i
  
        sum_phi_i_avg = sum_phi_i_avg + sum_phi_i*delta_time
        sum_phi_d_avg = sum_phi_d_avg + phi_d_sum*delta_time
        sum_phi_b2_avg = sum_phi_b2_avg + sum_phi_b2*delta_time
        sum_phi_d_check = sum_phi_d_check + phi_d_check*delta_time
       END IF
       END DO

        do j=1,NY_TOT
        do k=1,NZ
        do i=1,ni
         th_wh_avg(i,k,j)= th_wh_avg(i,k,j)/time_sum
         u1_avg(i,k,j)   = u1_avg(i,k,j)/time_sum
         u2_avg(i,k,j)   = u2_avg(i,k,j)/time_sum
         u3_avg(i,k,j)   = u3_avg(i,k,j)/time_sum
        end do
        end do
        end do

        sum_phi_i_avg = sum_phi_i_avg/time_sum
        sum_phi_d_avg = sum_phi_d_avg/time_sum
        sum_phi_b2_avg = sum_phi_b2_avg/time_sum
        sum_phi_d_check = sum_phi_d_check/time_sum
       
       WRITE(6,*) 'Saving Mean PE Statistics' 
       WRITE(6,*) 'sum_phi_i_avg = ', sum_phi_i_avg
       WRITE(6,*) 'sum_phi_d_avg = ', sum_phi_d_avg
       WRITE(6,*) 'sum_phi_b2_avg = ',sum_phi_b2_avg
       WRITE(6,*) 'sum_phi_d_check = ', sum_phi_d_check

     open(unit=41,file='mean_PE_statistics.txt',form='formatted',status='unknown')
       write(41,*) ' sum_phi_i_avg = ' , sum_phi_i_avg
       write(41,*) ' sum_phi_d_avg = ', sum_phi_d_avg
       write(41,*) ' sum_phi_b2_avg = ', sum_phi_b2_avg
       write(41,*) ' sum_phi_d_check = ', sum_phi_d_check
!666    format(a6,E17.8,a10,E17.8,a9,E17.8,a8,E17.8)
666     format(i5, 4E17.8)
        close(40)
        close(41)
        close(42)
        close(43)
111     format( E12.5,2E18.5)
        end

subroutine read_grid
use Domain
use ntypes
implicit none

integer i,j,k
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
allocate (GY(0:nj+1))
allocate (GZ(0:NZ+1))
allocate (DY(0:nj+1))
allocate (DZ(0:NZ+1))
allocate (GYF(0:nj+1))
allocate (GZF(0:NZ+1))
allocate (DYF(0:nj+1))
allocate (DZF(0:NZ+1))
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
np1 = ni
np2 = nj
np3 = NZ


         OPEN (30,file='zgrid.txt',form='formatted',status='old')
         READ (30,*) NZ_T
!C Check to make sure that grid file is the correct dimensions
         IF (NZ_T.ne.NZ) THEN
           WRITE(6,*) 'NZ, NZ_T',NZ,NZ_T
          !STOP 'Error: zgrid.txt wrong dimensions'
         END IF
         DO K=1,NZ_T+1
           READ(30,*) GZ(k+1)
         END DO
         DO K=1,NZ_T
           READ(30,*) GZF(k+1)
         END DO
!C Define ghost cells, if needed for this grid...
         GZF(1)=GZF(2) - (GZF(3)-GZF(2))
         GZF(0)=GZF(1) - (GZF(2)-GZF(1))
         GZF(NZ_T+2)=GZF(NZ_T+1)+GZF(NZ_T+1)-GZF(NZ_T)
         GZF(NZ_T+3)=GZF(NZ_T+2)+GZF(NZ_T+2)-GZF(NZ_T+1)
         GZ(1)=GZ(2) - (GZ(3)-GZ(2))
         GZ(0)=GZ(1) - (GZ(2)-GZ(1))
         GZ(NZ_T+2)=GZ(NZ_T+1)+GZ(NZ_T+1)-GZ(NZ_T)
         GZ(NZ_T+3)=GZ(NZ_T+2)+GZ(NZ_T+2)-GZ(NZ_T+1)
!C Define grid spacing
         DO K=1,NZ_T+1
           DZ(K)=(GZF(K)-GZF(K-1))
         END DO
         DO K=1,NZ_T
           DZF(K)=(GZ(K+1)-GZ(K))
           !write(400,*) DZF(K)
         END DO
         DZ(0)=DZ(1)
         DZF(NZ_T+2)=DZF(NZ_T+1)
         DZF(NZ_T+3)=DZF(NZ_T+2) 
         CLOSE(30)

         OPEN (30,file='ygrid.txt',form='formatted',status='old')
         READ (30,*) NY_T
!C Check to make sure that grid file is the correct dimensions
         IF (NY_T.ne.NY_TOT) THEN
           WRITE(6,*) 'NY, NY_T',NY_TOT,NY_T
          !STOP 'Error: ygrid.txt wrong dimensions'
         END IF
         DO J=1,NY_T+1
           READ(30,*) GY(j+1)
         END DO
         DO J=1,NY_T
           READ(30,*) GYF(j+1)
!           write(6,*) 'Grid GYF', GYF(j)
         END DO

         GYF(1)= GYF(2) - (GYF(3)-GYF(2))
         GYF(0)= GYF(1) - (GYF(2)-GYF(1))
         GYF(NY_T+2)= GYF(NY_T+1) + GYF(NY_T+1)-GYF(NY_T)         
         GYF(NY_T+3)= GYF(NY_T+2) + GYF(NY_T+2)-GYF(NY_T+1)
         GY(1)= GY(2) - (GY(3)-GY(2))
         GY(0)= GY(1) - (GY(2)-GY(1))
         GY(NY_T+2)= GY(NY_T+1) + GY(NY_T+1)-GY(NY_T)
         GY(NY_T+3)= GY(NY_T+2) + GY(NY_T+2)-GY(NY_T+1)

         CLOSE(30)

!         do j = 1,np2
!          g2vtk(j)=GYF(j)
!         enddo

         
!         do i = 1,np1
!          g1vtk(i)=LX/dble(ni-1)*dble(i-1)
!          write(6,*) g1vtk(i)
!         enddo

!         DX_s=g1vtk(2)-g1vtk(1)

!         do k = 1,np3
!          g3vtk(k)=GZF(k)
!          write(6,*) g1vtk(i)
!         enddo

!         DZ_s=g3vtk(2)-g3vtk(1)

return
end

!---------------------------------------
      subroutine PHI_cal
!---------------------------------------

!  Note, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps
!  This subroutine should be called in SAVE_STATS_CHAN after computing
!  plane averaged statistics, with the velocity in physical space, and 
!  CRi containing the velocity in Fourier space
use ntypes
use Domain
use qsort_c_module
implicit none

INTERFACE
   FUNCTION Heavi_side (r)
     REAL(8)             ::  Heavi_side
     REAL(8), INTENT(IN) :: r
   END FUNCTION Heavi_side
END INTERFACE


      integer i,j,k,n,ii,jj,kk,NP_I, tot_ind, COUNTER
      CHARACTER*31 file_energy
!      real(r8)  ::  rho_temp, rho_trans, rho_diff, DV
      real(r8)  ::   H_domain, sum_area ,area_test, height_distributed_hill
      real(r8)  ::   rho_max, rho_min, sum_vol,sum_dv, z_star_pre,BIG_NUM
      real(r8)  ::   temp_var(1:10,1:2)    
 
      real(r8),allocatable,dimension(:,:,:) :: S1
      real(r8),allocatable,dimension(:,:)   ::  rho_1D_tot
      real(r8),allocatable,dimension(:)     ::  drho_1D_totdz_s, drho_1D_totdz_s_check
      real(r8),allocatable,dimension(:)   :: area  
      parameter (BIG_NUM = 10.d0**8.0)

        allocate(area(NY_TOT))

      ! WRITE(6,*) 'Reading Hill Profile'

        H_hill_tot(:) = 0.0d0
        hill_ind_tot(:) = 0
        DX_hill(:) = 0.0d0
        area(:) = 0.0d0
        area_test = LX*LY

        open(202,file='hill_prof.dat',form='formatted',status='old')
        DO I=1,ni
        read(202,*)  DX_hill(I), H_hill_tot(I), hill_ind_tot(I)
        END DO
        close(202)

        sum_vol =0.0d0
        do j = 1, NY_TOT
        do i=1,ni 
        IF (j<hill_ind_tot(I)) THEN
         dv = LY*(GYF(j+1)-GYF(j))*DX_s
         sum_vol = sum_vol+dv
        END IF
        end do
        end do
        WRITE(6,*) 'hill_vol', sum_vol
        WRITE(6,*) 'NY_TOT', NY_TOT 
        height_distributed_hill =  (sum_vol/area_test)

        sum_area = 0.0d0      
        do j=1, NY_TOT
        do i=1, ni
        IF(j>hill_ind_tot(I)) THEN
        sum_area = sum_area+DX_s*LY
        END IF
        end do
        area(j) = sum_area
        sum_area = 0.0d0
      end do
        sum_vol = 0.0d0
        do j=1,NY_TOT
        sum_vol = sum_vol + (GYF(j+1)-GYF(j))*area(j)
        end do

                WRITE(6,*) 'sum_vol_area = ', sum_vol
        open(203,file='area.txt',form='formatted',status='unknown')
        DO J=1,NY_TOT
        write(203,*) GYF(J), area(J)
        END DO
        close(203)
  
        write(6,*) 'GYF(NY_TOT) =', GYF(NY_TOT)
      H_domain = GYF(NY_tot) - GYF(1)

      tot_ind = NY_TOT*NZ*NX

      write(6,*) 'Calculating Potential Energy Budget'

!Started  z* calculation
       allocate (rho_1D_tot(1:tot_ind,1:2))
       allocate (drho_1D_totdz_s(1:tot_ind))
       allocate (drho_1D_totdz_s_check(1:tot_ind))

       drho_1D_totdz_s_check(:) = 0.0d0
       drho_1D_totdz_s(:)       = 0.0d0
       rho_1D_tot(:,:)          = 0.0d0


      rho_1D_tot(:,1) = 100.d0 ;
      rho_1D_tot(:,2) = -10.d0 ;

 !    H_domian = dble(rank)

      jj = 0

       do i=1,ni 
         do k=1,NZ
           do j=1,NY_TOT
           jj =  jj + 1
           ii =  NY_tot*NZ*(i-1) + NY_tot*(k-1) + j   
        ! Calculate unique grid index
        IF(j<=hill_ind_tot(I)) THEN    
         TH_WH(i,k,j) = -2.0d0
         rho_1D_tot(jj,1) = -5.0d0
        ELSE
        rho_1D_tot(jj,1) = rho_0*(-alpha_T*TH_WH(i,k,j) + 1.0d0)!+ THBAR(j,N_TH+1)          ! Storing the rho in a 1D array
        END IF 
       rho_1D_tot(jj,2) =  dble(ii)  ! dble(tot_ind*rank + jj)             ! Storing the index for each grid to a 1D array
          enddo
         enddo
       enddo



!       sorting for rho values : keep the grid index
        call QsortC(rho_1D_tot(:,1),rho_1D_tot(:,2))

!        write(6,*) 'I am here', jj, ii
!        DO ii=1,tot_ind
!         write(77,666) dble(ii), rho_1D_tot(ii,1)
!        ENDDO
!        close(77)

   
        kk = tot_ind
        DO ii = 1,tot_ind
         drho_1D_totdz_s_check(ii) = rho_1D_tot(ii,1) !!!!!!!!!!!!!!!!!!!!!!!1

!         if ( rho_1D_tot(ii,1) .gt. 90.0d0 ) then
!          kk = ii-1
!          goto 212
!         else
!        WRITE(6,*) 'rho_1D_tot(ii+1,1) is out of bounds, STOPPING'
!        stop
! steps to calculate dz*/drho
! keeping 1/dp(ii) information
 
!          sum_vol = rho_1D_tot(ii+1,1)-rho_1D_tot(ii,1)
!          if ( abs(sum_vol) .lt. 1.d0/BIG_NUM ) then
!            if (sum_vol .lt. 0.0d0 ) then
!             drho_1D_totdz_s(ii) = - BIG_NUM 
!            else
!             drho_1D_totdz_s(ii) =   BIG_NUM
!            endif 
!          else
!           drho_1D_totdz_s(ii) = 1.0d0/sum_vol
!          endif
!         endif

        ENDDO

!212      continue


  !      WRITE(6,*) 'start_rho', rho_1D_tot(tot_ind,1), rho_1D_tot(1,1)
       rho_1D_tot(:,1) = 0.0d0;
         sum_vol = 0.0d0 ;
        z_star_pre = 0.0d0 ;
        jj = 1.0d0;
        sum_dv = 0.0d0;
  
        open (unit=42,file='z_star_rho.txt',status='unknown',form='formatted')
         DO ii = kk, 1,-1 
!           NP_I = floor(rho_1D_tot(ii,2)/tot_ind)
              i = floor((rho_1D_tot(ii,2))/(NZ*NY_TOT)) 
              k = floor((rho_1D_tot(ii,2) - i*(NZ*NY_TOT))/NY_TOT) 
              j = floor(rho_1D_tot(ii,2) - i*(NZ*NY_TOT) - k*NY_TOT) + 1
                
              DV  = DX_s*abs(GZF(k+1)-GZF(k))*abs(GYF(j+1)-GYF(j))
!             calculating z* at different ii level and  mapped back to physical grid(i,j,k)
           
      !      rho_1D_tot(int(rho_1D_tot(ii,2)),1) = (sum_vol + 0.50d0*DV)/area_test
      !      temp_val = (sum_vol + 0.50d0*DV)/area_test
      !      sum_vol =  sum_vol +  DV  

             IF (sum_dv+0.50d0*DV<abs(GYF(jj+1)-GYF(jj))*area(jj)) THEN
             ! IF (sum_vol+0.50d0*DV<GYF(jj)*area(jj)) THEN
             rho_1D_tot(int(rho_1D_tot(ii,2)),1) = (sum_vol + 0.50d0*DV)/area(jj)
             !rho_1D_tot(int(rho_1D_tot(ii,2)),1) = (sum_dv + 0.50d0*DV)/area(jj)
             temp_val = (sum_vol + 0.50d0*DV)/area(jj)
             ELSE
      !       WRITE(42,*) jj, sum_dv
             sum_dv=0.0d0
             jj = jj+1
             IF (jj>NY_TOT-1) THEN
             jj = NY_TOT-1
             END IF
             rho_1D_tot(int(rho_1D_tot(ii,2)),1) = (sum_vol + 0.50d0*DV)/area(jj+1)
             !rho_1D_tot(int(rho_1D_tot(ii,2)),1) = (sum_dv + 0.50d0*DV)/area(jj+1)
             temp_val = (sum_vol + 0.50d0*DV)/area(jj+1)
             END IF
             
             sum_vol =  sum_vol +  DV  
             sum_dv = sum_dv + DV
      
             !IF (jj<6) THEN
      !  WRITE(42,*) jj, sum_dv
             !END IF
       
!             calculating dz/drho at different ii level : not yet mapped to physical grid(i,j,k) 
!              drho_1D_totdz_s(ii) =  drho_1D_totdz_s(ii)*(z_star_pre-(sum_vol + 0.50d0*DV)/area(j))  

           !   z_star_pre =   (sum_vol + 0.50d0*DV)/area(j)  !  storing the z* for (ii) level to calculate z*(ii)-z*(ii-1)
!              write(66,666)dble(j),dble(ii),drho_1D_totdz_s_check(ii),(sum_vol + 0.50d0*DV)/area, drho_1D_totdz_s(ii)
         ENDDO
!         close(66)
        WRITE(6,*) 'sum_vol_loop', sum_vol
!        interpolating dz*drhp at last point of the array (here is kk)
         drho_1D_totdz_s(kk) = 2.0d0*drho_1D_totdz_s(kk-1) - drho_1D_totdz_s(kk-2)  

        open(unit=77,file='z_star_height.txt',status='unknown', form='formatted')

         allocate (S1(1:ni,1:NZ,1:NY_tot))
          jj = 0
         do i=1,ni
         do k=1,NZ
         do j=1,NY_tot
              jj =  jj + 1             
              S1(I,K,J) = rho_1D_tot(jj,1)   
         enddo
         enddo
         enddo

 
        do k=1,NZ
            write(77,222) GZF(K), S1(256,k,NY_TOT),S1(120,k,NY_TOT),S1(512,k,NY_TOT)
            enddo
        close(77)

        do j=1,NY_tot
         z_start_bar(j) = z_start_bar(j) + sum(S1(1:ni,1:NZ,J))/dble(ni*NZ)*delta_time
        enddo 
        
        do i=1,ni
          do k=1,NZ
            z_start_bar_surf(i,k) = z_start_bar_surf(i,k) + S1(i,k,NY_tot)*delta_time           
          enddo
        enddo
         
        phi_d_prof(:) = 0.0d0
        phi_b2_prof(:) = 0.0d0
        phi_d_sum = 0.0d0
        temp_val = 0.0d0
        sum_phi_b2 = 0.0d0
        phi_d_check = 0.0d0

! Calculation of Laplacian Density term

        do i=2,ni-1
        do j=hill_ind_tot(I)+2,NY_TOT-2
        do k=2,NZ-2

        dv = DX_s*DZF(k)*(GYF(j+1)-GYF(j))
        dv_z= DX_s*DZF(k)

!        sum_val = sum_val + dv

            temp_val =  (S1(i,k,j)*(((((TH_WH(i,k,j)-TH_WH(i-1,k,j))/(DX_s))&
                        -((TH_WH(i+1,k,j)-TH_WH(i,k,j))/(DX_s)))&
                        /(2.0d0*DX_s))&         
                        +((((TH_WH(i,k,j)-TH_WH(i,k,j-1))/(GYF(j)-GYF(j-1)))&
                        -((TH_WH(i,k,j+1)-TH_WH(i,k,j))/(GYF(j+1)-GYF(j))))&
                        /(GYF(j+1)-GYF(j-1)))&
                        +((((TH_WH(i,k,j)-TH_WH(i,k-1,j))/(DZF(k-1)))&
                        -((TH_WH(i,k+1,j)-TH_WH(i,k,j))/(DZF(k))))&
                        /(2.0d0*DZF(k-1)))))

        phi_d_check = phi_d_check+rho_0*-alpha_T*kappa*Gravity*temp_val*dv
  
        end do
        end do
        end do

       !WRITE(6,*) 'phi_d_check =', phi_d_check

        do i=2,ni-1
        do j=NY_TOT-2,NY_TOT-1
        do k=2,NZ-2

        dv = DX_s*DZF(k)*(GYF(j+1)-GYF(j))
        dv_z= DX_s*DZF(k)

!        sum_val = sum_val + dv

        temp_val =  (S1(i,k,j)*(((((TH_WH(i,k,j)-TH_WH(i-1,k,j))/(DX_s))&
                    -((TH_WH(i+1,k,j)-TH_WH(i,k,j))/(DX_s)))&
                    /(2.0d0*DX_s))&
                    +((((TH_WH(i,k,j)-TH_WH(i,k,j-1))/(GYF(j)-GYF(j-1)))&
                    -((TH_WH(i,k,j+1)-TH_WH(i,k,j))/(GYF(j+1)-GYF(j))))&
                    /(GYF(j+1)-GYF(j-1)))&
                    +((((TH_WH(i,k,j)-TH_WH(i,k-1,j))/(DZF(k-1)))&
                    -((TH_WH(i,k+1,j)-TH_WH(i,k,j))/(DZF(k))))&
                    /(2.0d0*DZF(k-1)))))

        phi_d_check = phi_d_check+rho_0*-alpha_T*kappa*Gravity*temp_val*dv

        end do
        end do
        end do

        WRITE(6,*) 'phi_d_check =', phi_d_check


        !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
        !           phi_d calculation            !
        !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

        do i=2,ni-1
        do j=hill_ind_tot(I)+2,NY_TOT-2
        do k=2,NZ-2
  
!      calculation of buoyancy flux terms

        !dv = DX_s*(GZF(k+1)-GZF(k))*(GYF(j+1)-GYF(j))
        !dv_z= DX_s*(GZF(k+1)-GZF(k))

        !sum_val = sum_val + dv
!        IF (GYF(j)<(0.2d0-height_distributed_hill)) THEN
        
        dv = DX_s*(GZF(k+1)-GZF(k))*(GYF(j+1)-GYF(j))
        dv_z= DX_s*(GZF(k+1)-GZF(k)) 

         temp_val  = alpha_T*((S1(i+1,k,j)-S1(i-1,k,j))/&
                       (2.0*DX_s))*((TH_WH(i+1,k,j)-TH_WH(i-1,k,j))/&
                       (2.0*DX_s))&                             
                       +alpha_T*((S1(i,k,j+1)-S1(i,k,j-1))&
                       /(GYF(j+1)-GYF(j-1)))*((TH_WH(i,k,j+1)-TH_WH(i,k,j-1))&
                       /(GYF(j+1)-GYF(j-1)))  &                      
                       +alpha_T*((S1(i,k+1,j)-S1(i,k-1,j))/&
                       (GZF(k+1)-GZF(k-1)))*((TH_WH(i,k+1,j)-TH_WH(i,k-1,j))/&
                       (GZF(k+1)-GZF(k-1)))  
                             
           phi_d_prof(J) = phi_d_prof(J)  + rho_0*Gravity*kappa*temp_val*dv_z
           phi_d_sum = phi_d_sum + rho_0*Gravity*kappa*temp_val*dv
        !WRITE(6,*) 'TH_WH(j-1) =', TH_WH(i,k,j-1)     
 !      END IF
        end do
        end do
        end do
        
        do i=2,ni-1
        do j=NY_TOT-2, NY_TOT-1
        do k=2,NZ-2

      !calculation of buoyancy flux terms

        dv = DX_s*(GZF(k+1)-GZF(k))*(GYF(j+1)-GYF(j))
        dv_z= DX_s*(GZF(k+1)-GZF(k))

        sum_val = sum_val + dv
  !      IF (GYF(j)<(0.2d0-height_distributed_hill)) THEN
        temp_val  = alpha_T*((S1(i+1,k,j)-S1(i-1,k,j))/&
                        (2.0*DX_s))*((TH_WH(i+1,k,j)-TH_WH(i-1,k,j))/&
                        (2.0*DX_s))&
                        +alpha_T*((S1(i,k,j+1)-S1(i,k,j-1))&
                        /(GYF(j+1)-GYF(j-1)))*((TH_WH(i,k,j+1)-TH_WH(i,k,j-1))&
                        /(GYF(j+1)-GYF(j-1)))  &
                        +alpha_T*((S1(i,k+1,j)-S1(i,k-1,j))/&
                        (GZF(k+1)-GZF(k-1)))*((TH_WH(i,k+1,j)-TH_WH(i,k-1,j))/&
                        (GZF(k+1)-GZF(k-1)))

        phi_d_prof(J) = phi_d_prof(J)  + rho_0*Gravity*kappa*temp_val*dv_z
        phi_d_sum = phi_d_sum + rho_0*Gravity*kappa*temp_val*dv
        !WRITE(6,*) 'TH_WH(j-1) =', TH_WH(i,k,j-1)
  !      END IF
        end do
        end do
        end do



        WRITE(6,*) 'phi_d = ', phi_d_sum
       ! WRITE(6,*) 'TH_WH(j-1) =', TH_WH(i,k,j-1)

        !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
        !           phi_b2 calculation           !
        !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

        temp_val = 0.0d0

!        do i=2,ni-1
!        do j=hill_ind_tot(I)+3,NY_TOT-2
!        do k=2,NZ-2

!        dv_z= DX_s*DZF(k)

!            temp_val  = alpha_T*(S1(i,k,j+1)*((TH_WH(i,k,j+2)-TH_WH(i,k,j))/(GYF(j+2)-GYF(j))) &
!                               - S1(i,k,j-1)*((TH_WH(i,k,j)-TH_WH(i,k,j-2))/(GYF(j)-GYF(j-2))))/(GYF(j+1)-GYF(j-1))

!            phi_b2_prof(J) = phi_b2_prof(J)  -  rho_0*Gravity*kappa*temp_val*dv_z

!          enddo
!          enddo
!          enddo
        
!           phi_b2_prof(126) = 2.0d0*phi_b2_prof(125) - phi_b2_prof(124)
!           phi_b2_prof(127) = 2.0d0*phi_b2_prof(126) - phi_b2_prof(125)

           phi_d_prof(126) = 2.0d0*phi_d_prof(125) - phi_d_prof(124)
           phi_d_prof(127) = 2.0d0*phi_d_prof(126) - phi_d_prof(125)

           temp_val = 0.0d0

        do k=1,NZ
        do i=1,ni

        dv_z= DX_s*DZF(k-1)

           sum_phi_b2 = sum_phi_b2 + rho_0*alpha_T*Gravity*kappa*S1(i,k,NY_TOT-1)*((TH_WH(i,k,NY_TOT-1)-TH_WH(i,k,NY_TOT-2))/(GYF(NY_TOT-1)-GYF(NY_TOT-2)))*DV_z                  
        !WRITE(6,*) 'S1(i,k,NY_TOT-1)', S1(i,k,NY_TOT-1)
       
        enddo
        enddo  
        write(6,*) 'phi_b2 =', sum_phi_b2
       !write(6,*) 'phi_b2_prof((126+127)/2) =', (phi_b2_prof(126)+phi_b2_prof(127))/2
        sum_phi_i = 0.0d0
        sum_rho_h = 0.0d0
        sum_rho_top = 0.0d0
        temp_val = 0.0d0

        !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
        !           phi_i calculation            !
        !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

        do k=1,NZ
        do i=1,NX
        temp_val = -th_wh(i,k,hill_ind_tot(i)+1)*DX_s*DZF(k-1)
        sum_rho_h = temp_val+sum_rho_h
        end do
        end do
        temp_val = 0.0d0

        do k=1,NZ
        do i=1,NX
        temp_val = -th_wh(i,k,NY_TOT-1)*DX_s*DZF(k-1)
        sum_rho_top = temp_val+sum_rho_top      
        end do
        end do
        sum_phi_i = sum_phi_i+alpha_T*rho_0*Gravity*kappa*(sum_rho_top-sum_rho_h)
        WRITE(6,*) 'phi_i =', sum_phi_i
        !WRITE(6,*) 'th_wh(i,k,NY_TOT) =', th_wh(i,k,NY_TOT)

     deallocate ( area,S1, rho_1D_tot, drho_1D_totdz_s, drho_1D_totdz_s_check )

222 format(E17.8, 3E17.8)

     return
     end
