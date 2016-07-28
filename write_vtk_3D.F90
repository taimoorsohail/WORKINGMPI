         module ntypes
         integer, parameter :: r4=4
         integer, parameter :: r8=8
         integer, parameter :: i4=4
        end module ntypes

!DOMAIN
        module Domain
         use ntypes

         integer(i4)       :: NX, NY, NZ, N_TH, NXM, NYM, NZM, TNKZ, TNKY
         integer(i4)       :: NKX,s11
         integer(i4)       :: NXP,NZP,NXV,NZV, NKXV, NKXP,NX2V,NX2P,NXP_L,NX2P_L

         integer(i4)       :: ni,nj,np,mm,nk_en,nk,y_1,y_2,y_3, nk_st,ii

        parameter (nk_st=10672, nk_en=10764)

!         parameter (nk=NZ+2, nj=NY+2, ni=NX)
         !INCLUDE   'grid_def'

         integer(i4)       :: np1,np2, np3
!         parameter(np1 = ni, np3 = nj-2, np2 = nk-2)
        end module Domain
!GRID       
        module Grid
        use ntypes


        real(r8)          :: LX, LY, LZ, CSX, CSY, CSZ         !Length
        real(r8)          :: TH_3D, U_1, V_1, W_1, PR_TOT
         INTEGER           :: jstart, jend , ZSTART,  ZEND
!

        real(r8),allocatable,dimension(:) :: GX,  GY,  GZ,  DX,  DY,  DZ, &
        GXF, GYF, GZF, DXF, DYF, DZF

        real(r8),allocatable,dimension(:,:)   :: xpoint, ypoint
        end module Grid

!run_variable
        module run_variable
        use ntypes
        use Domain

        logical   :: READING_FLOW, grid_write
        real(r4),allocatable,dimension(:,:,:) :: g1vtk(:),g2vtk(:),g3vtk(:), &
                                                 var_1(:),var_2(:)

!       g1vtk(1:ni),g2vtk(1:nk),g3vtk(1:nj), var_1(1:3*np1*np2*np3)
        real(r4)  :: xL
!real(r4)  :: vrms(1:np1,1:np2),wrms(1:np1,1:np2), &
!                     u1(1:np1,1:np2),u2(1:np1,1:np2), u3(1:np1,1:np2),
!                     &
!                     th_div(1:np1,1:np2,1:N_TH),th_wh(1:np1,1:np2,1:N_TH+1),
!                     &
!                     cp(1:np1,1:np2)


      Real(r8) DX_s,DZ_s,time_old, time_sum, delta_time,time,dt,ubulk
      real(r8),allocatable,dimension(:,:,:) :: x,y,z,u,v,w,u1,u2,u3
      real(r8),allocatable,dimension(:,:,:) :: th,pr_3D
      real(r8),allocatable,dimension(:,:,:) :: x_c,y_c,z_c,th_wh_c, &
                                             u2_c,u3_c,f1_c
      real(r8),allocatable,dimension(:,:,:) :: u1_avg,u2_avg,u3_avg,s1

      real(r8),allocatable,dimension(:,:) :: u2_avg_2D, th_avg_2D,thv_avg_2D
      real(r8),allocatable,dimension(:,:,:) :: varSP_zstar
      real(r8),allocatable,dimension(:,:) ::  rho_1D_tot

      real(r4),allocatable,dimension(:,:,:) :: th_wh,th_wh_avg, omega_x
      Real(r8) ::  th_b
      real(r8),allocatable,dimension(:)      ::  H_hill_tot, DX_hill
      integer, allocatable,dimension(:)      ::  hill_ind_tot
      real(r8) ::  u_0,Gravity,alpha_T, NU, rho_0
      parameter (u_0=0.125, Gravity = 10.0d0,  alpha_T = 2.0d0*10.0d0**(-4.0),&
         NU = 10.0**(-6), rho_0 = 1000.0d0)
      real(r8) ::   sum_Prod, sum_dissp, sum_buoyF, sum_Prod_avg, sum_dissp_avg, &
                        sum_buoyF_avg, sum_wstress,  &
      sum_dtkdt, sum_vis, sum_trans, sum_advec, dv,ds, vol_box

      integer  :: i_1, j_1, i_2, i_3,j_3
      parameter (i_1 = 152, j_1 = 45, i_2 = i_1, i_3 = i_1, j_3 =118 )

      integer ::  imax,jmax,kmax,pp,jj,p_p,p_m
      integer ::  debug,ier,itot ! local_Y(3,ni)
      integer ::  tecini,tecdat,teczne,tecnod,tecfil,tecend
      integer ::  visdouble,disdouble,i_strt,i_end,j_strt,j_end
      character*1 nulchar
!      real(r8) :: var_1_spectrum(2*np),var_2_spectrum(2*np)


!      CHARACTER*34 file_name_1,file_name_3
!      CHARACTER*18 file_name_2
!      CHARACTER*31 file_name_4

        LOGICAL   ::  TKE_FIELD, WRITE_TECP, X_PLANE,Y_PLANE,INTER_pol_y, &
         INTER_pol_z, vorticity_cal, PE_cal, WRITE_PARA_View,CAL_STAT_AVG, IBM

        end module run_variable

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

!write(6,*) 'MADE IT INTO QSORTC'
!write(6,*) 'iq = ', iq

! write(6,*) 'INSIDE QsortC: rho_1D_tot(5000,1) = ', rho_1D_tot(5000,1)

        if (size(A) > 1) then
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
         i = 0
         j = size(A) + 1

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

         Program main
        use ntypes
        use Domain
        use Grid
        use run_variable
      implicit none

      integer :: i,j,k,kk,p

!     integer np1,np2,np3
       !FILENAMES
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss

!      real(r4),allocatable,dimension(:) :: g1vtk,g2vtk,g3vtk
!      real(4) :: g1vtk(1:np1),g3vtk(1:np2),g2vtk(1:np3), var_1(1:3*np1*np2*np3)
!      real(4) :: th_wh(np1,np3,np2)
       real(8) :: sum_1, sum_2
              CHARACTER*5 file_num
              CHARACTER*39 file_name_1,file_time
              CHARACTER*16 base_name
          character(len=18) :: title

!        X_PLANE    = .FALSE.
!        Y_PLANE    = .FALSE.

!        WRITE_TECP = .FALSE.
        WRITE_PARA_View = .FALSE.
        IBM = .FALSE.
        PE_cal = .FALSE.
        CAL_STAT_AVG = .TRUE.
        vorticity_cal = .FALSE.
!        Inter_pol_y  = .FALSE.
!        Inter_pol_z  = .FALSE.

       outputDIR='plane_data_3D/'
       basename = 'data_binary_3D'

       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR), &
       trim(basename)//"_n",nk_st,".pln"

       open(unit=22,file=OutFileName,access='stream',  &
       form='unformatted',status='old',             &
       convert='big_endian',iostat=s11)

         write(6,*)'Reading time and grid data', nk_st,' ', OutFileName

        read(22) time, ni, nk, nj, DX_s, LX, LY, LZ
        close(22)
        
        NZ = nk-2
        NX = ni
        NY = nj-2
        np = ni+2
        
        WRITE(6,*) 'NX',NX,'NY',NY,'NZ',NZ,'np',np
        WRITE(6,*) 'TIME = ',TIME
        
        IF (CAL_STAT_AVG) THEN
        
       allocate(th_wh_avg(nk,nj,np),u1_avg(nk,nj,np),u2_avg(nk,nj,np),u3_avg(nk,nj,np))
       allocate(s1(nk,nj,np))

       th_wh_avg(:,:,:)= 0.0d0
       u1_avg(:,:,:)   = 0.0d0
       u2_avg(:,:,:)   = 0.0d0
       u3_avg(:,:,:)   = 0.0d0
       s1(:,:,:)       = 0.0d0

       time_sum = 0.0d0
      ENDIF

        IF (IBM) THEN
        WRITE(6,*) 'Reading Hill Profile'

        allocate(H_hill_tot(ni), hill_ind_tot(ni), DX_hill(ni))
      
        H_hill_tot(:) = 0.0d0
        hill_ind_tot(:) = 0
        DX_hill(:) = 0.0d0

        open(202,file='hill_prof.dat',form='formatted',status='old') 
        DO I=1,ni
        read(202,*)  DX_hill(I), H_hill_tot(I), hill_ind_tot(I)   
        END DO
        close(202)
        
!        open(333,file='hill_prof_test.txt',form='formatted',status='unknown')
!        DO I=1,ni
!         write(333,111) DX_hill(I), H_hill_tot(I), hill_ind_tot(I)
!        END DO
!111     format(2f12.6,i5)        
!        close(333)

        END IF
      
        allocate (u(ni,nk,nj),v(ni,nk,nj),w(ni,nk,nj), &
                 th(ni,nk,nj), th_wh(nk,nj,np), pr_3D(ni,nk,nj),u1(nk,nj,np), &
                 u2(nk,nj,np),u3(nk,nj,np))
        allocate(u2_avg_2D(nk,nj), th_avg_2D(nk,nj), thv_avg_2D(nk,nj))
 
        u(:,:,:) = 0.0d0
        v(:,:,:) = 0.0d0
        w(:,:,:) = 0.0d0
        th(:,:,:) = 0.0d0
        th_wh(:,:,:) = 0.0d0
        pr_3D(:,:,:) = 0.0d0
        u1(:,:,:) = 0.0d0
        u2(:,:,:) = 0.0d0
        u3(:,:,:) = 0.0d0
        

     write(6,*) 'Grid is reading ', ni

      call read_grid

      write(6,*) 'Grid reading done', ni

       th_b = 30.0d0

       allocate(omega_x(nk,nj,np))
       omega_x(:,:,:) = 0.d0 

       DO kk=nk_st,nk_en      
       outputDIR='plane_data_3D/'
       basename = 'data_binary_3D'

       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR), &
       trim(basename)//"_n",kk,".pln"

       open(unit=22,file=OutFileName,access='stream',  &
       form='unformatted',status='old',             &
       convert='big_endian',iostat=s11)

         write(6,*)'Reading flow statistics', kk,' ',OutFileName
        read(22) TIME, ni, nk, nj, DX_s, LX, LY, LZ
        read(22) (((th(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        read(22) (((u(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        read(22) (((w(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        read(22) (((v(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        read(22) (((pr_3D(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        close(22)

     !   open(66,file='plane_data_3D/time_bulk_3D_test.txt',form='formatted', &
     !  status='unknown')        
        
     !   DO I = 1,NX
     !   TH_3D = th(I, 124,113)
     !   U_1 = u(I, 124,113)
     !   W_1 =  w(I, 124,113)
     !   V_1 = v(I,124,113)
     !   PR_TOT = pr_3D(I,124,113)
     !   write(66,565) DX_s*(I-1) , TH_3D, U_1, W_1, V_1, PR_TOT
     !   END DO
!565 !   format(f12.5,5f13.8)
     !   close(66)

        WRITE(6,*) 'Transforming Co-ordinate System'

      DO p=1,np-2
        DO j=1,nj
          DO k=1,nk
             th_wh(k,j,p)=th(p,k,j) + th_b !(-y(i,j,np)+y(1,nj,np))
             u1(k,j,p)   = v(p,k,j)
             u2(k,j,p)   = w(p,k,j)
             u3(k,j,p)   = u(p,k,j)
          END DO
        END DO
       END DO

!       DO p=np-1,np
        p=np-1
        DO j=1,nj
          DO k=1,nk
             u1(k,j,p)   = u1(k,j,1)
             u2(k,j,p)   = u2(k,j,1)
             u3(k,j,p)   = u3(k,j,1)
             th_wh(k,j,p)=th_wh(k,j,1)
          END DO
        END DO

        p=np
        DO j=1,nj
          DO k=1,nk
             u1(k,j,p)   = u1(k,j,2)
             u2(k,j,p)   = u2(k,j,2)
             u3(k,j,p)   = u3(k,j,2)
             th_wh(k,j,p)=th_wh(k,j,2)
          END DO
        END DO

      if (vorticity_cal) then

      WRITE(6,*) 'Calculating Vorticity'

      do p=1,np
       p_p = p+1
       p_m = p
       if (p==np) THEN
        p_p = 1
       endif
        do j=1,nj-1
        do k=1,nk
         omega_x(k,j,p)= (u2(k,j+1,p_m)-u2(k,j,p_m))/(GYF(j)-GYF(j-1)) &
                   -(u3(k,j,p_p)-u3(k,j,p_m))/DX_s

      end do
      end do
      end do

      do jj = 1,2
      do p=1,np
      do j=1,nj
      do k=2,nk
        omega_x(k,j,p)= 0.5*(omega_x(k,j,p) + omega_x(k-1,j,p))
      end do
      end do
      end do
      end do

!      open(333,file='vorticity_test.txt',form='formatted', &
!                                                     status='unknown')
!      do p=1,np
!      do j=nj/2
!      do k=nk/2
!                write(333,111) omega_x(k,j,p)
!      end do
!      end do
!      end do
!111   format(2f12.6,i5)
!      close(333)

      WRITE(6,*) 'Done calculating Vorticity'

      endif

      IF (WRITE_PARA_View) THEN

      WRITE(6,*) 'Writing Instantaneous Paraview'

      call plane_parav_vel(kk)

      ENDIF

        IF( (CAL_STAT_AVG).AND. (kk .gt. nk_st) ) THEN

        WRITE(6,*) 'kk', kk, '> nk_st', nk_st

       delta_time = time-time_old
       time_sum = time_sum + delta_time
        
       do p=1,np
       do j=1,nj
       do k=1,nk
        th_wh_avg(k,j,p)= th_wh_avg(k,j,p)+th_wh(k,j,p)*delta_time
        u1_avg(k,j,p)= u1_avg(k,j,p) + u1(k,j,p)*delta_time
        u2_avg(k,j,p)= u2_avg(k,j,p) + u2(k,j,p)*delta_time
        u3_avg(k,j,p)= u3_avg(k,j,p) + u3(k,j,p)*delta_time
       end do
       end do
       end do     
         ENDIF
        
      time_old = time
      END DO
      
         IF (CAL_STAT_AVG) THEN
       do p=1,np
       do j=1,nj
       do k=1,nk
        th_wh_avg(k,j,p)=th_wh_avg(k,j,p)/time_sum
        u1_avg(k,j,p)= u1_avg(k,j,p)/time_sum
        u2_avg(k,j,p)= u2_avg(k,j,p)/time_sum
        u3_avg(k,j,p)= u3_avg(k,j,p)/time_sum

        s1(k,j,p) = 0.5d0*((u1_avg(k,j,p))**2+(u2_avg(k,j,p))**2+&
                                                (u3_avg(k,j,p))**2)
       end do
       end do
       end do

        WRITE(6,*) 'Writing average paraview'

       call plane_parav_vel(1)
   
       write(6,*) 'Done with averaging'

       write(6,*) 'Opening tke_budget.dat, mean_tke_budget.dat & 
                                                and mean_ke_budget.dat'

       open(70,file='tke_budget.dat',form='formatted', &
       status='unknown')

       open(69,file='mean_tke_budget.dat',form='formatted', &
       status='unknown')

       open(67,file='mean_ke_budget.dat',form='formatted', &
       status='unknown')

! Calculating mean KE budget

         call call_stat_cal_avg

        write(6,*) 'time = ', time, &
                   'sum_dissp = ', sum_dissp, &
                    'sum_buoyF = ', sum_buoyF, &
                    'sum_advec = ', sum_advec, &
                    'sum_wstress = ', sum_wstress
        write(67,112) time, sum_dissp, sum_buoyF, sum_advec, sum_wstress
                !sum_dissp_1,sum_dissp_2,sum_dissp_3,sum_dissp_4, &
                !sum_buoyF_1, sum_buoyF_2, sum_buoyF_3, sum_buoyF_4

                  !,sum_advec_1,sum_advec_2,sum_advec_3,sum_advec_4, &
                  !sum_advec_1+sum_advec_2 +sum_advec_3 +sum_advec_4
                
! k starts at nk_st+1 because it is subtracted from the mean, which only has a
! value at nk_st+1
       
       do kk=nk_st+1,nk_en

       outputDIR='plane_data_3D/'
       basename = 'data_binary_3D'

       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR), &
       trim(basename)//"_n",nk_st+1,".pln"

       open(unit=22,file=OutFileName,access='stream',  &
       form='unformatted',status='old',             &
       convert='big_endian',iostat=s11)

!        write(6,'(i6,2a)') nk_st+1,' ', OutFileName

        read(22) time, ni, nk, nj, DX_s, LX, LY, LZ
        close(22)

        NZ = nk-2
        NX = ni
        NY = nj-2
        np = ni+2

        deallocate(u,v,w,th)

       allocate (u(ni,nk,nj))
       allocate (v(ni,nk,nj))
       allocate (w(ni,nk,nj))
       allocate (th(ni,nk,nj))

       outputDIR='plane_data_3D/'
       basename = 'data_binary_3D'

       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR), &
       trim(basename)//"_n",kk,".pln"

       open(unit=22,file=OutFileName,access='stream',  &
       form='unformatted',status='old',             &
       convert='big_endian',iostat=s11)

!         write(6,'(i6,2a)') kk,' ',OutFileName
        read(22) TIME, ni, nk, nj, DX_s, LX, LY, LZ
        read(22) (((th(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        read(22) (((u(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        read(22) (((w(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        read(22) (((v(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        read(22) (((pr_3D(I,K,J),I=1,ni),K=1,nk),J=1,nj)
        close(22)

      DO p=1,np-2
        DO j=1,nj
          DO k=1,nk
             th_wh(k,j,p)=th(p,k,j) + th_b !(-y(i,j,np)+y(1,nj,np))
             u1(k,j,p)   = v(p,k,j)
             u2(k,j,p)   = w(p,k,j)
             u3(k,j,p)   = u(p,k,j)
          END DO
        END DO
       END DO

!       DO p=np-1,np
        p=np-1
        DO j=1,nj
          DO k=1,nk
             u1(k,j,p)   = u1(k,j,1)
             u2(k,j,p)   = u2(k,j,1)
             u3(k,j,p)   = u3(k,j,1)
             th_wh(k,j,p)=th_wh(k,j,1)
          END DO
        END DO

        p=np
        DO j=1,nj
          DO k=1,nk
             u1(k,j,p)   = u1(k,j,2)
             u2(k,j,p)   = u2(k,j,2)
             u3(k,j,p)   = u3(k,j,2)
             th_wh(k,j,p)=th_wh(k,j,2)
          END DO
        END DO

        deallocate(u,v,w,th)

       allocate (u(nk,nj,np))
       allocate (v(nk,nj,np))
       allocate (w(nk,nj,np))
       allocate (th(nk,nj,np))

!      Calculating tke_budget
       call call_stat_cal

        IF (kk.eq.nk_st+1) THEN
        delta_time = 2.00000000004366
        ELSE       
        delta_time = time-time_old
        END IF
      
        !  delta_time = time-time_old
       
     WRITE(6,*) time, delta_time, sum_dtkdt, sum_advec, sum_Prod, sum_dissp, sum_buoyF

     write(70,113) time, delta_time, sum_dtkdt, sum_advec, sum_Prod, sum_dissp, sum_buoyF

       time_sum = time_sum + delta_time

        sum_Prod_avg=  sum_Prod_avg+sum_Prod*delta_time
        sum_dissp_avg= sum_dissp_avg+sum_dissp*delta_time
        sum_buoyF_avg= sum_buoyF_avg+sum_buoyF*delta_time

        time_old = time
       
!       write(68,113)TIME,sum_dtkdt,sum_advec, sum_Prod_1, sum_dissp_1,
!sum_buoyF_1
!       write(69,113)TIME,sum_dtkdt,sum_advec, sum_Prod_2, sum_dissp_2,
!sum_buoyF_2
!       write(70,113)TIME,sum_dtkdt,sum_advec, sum_Prod_3, sum_dissp_3,
!sum_buoyF_3
!       write(71,113)TIME,sum_dtkdt,sum_advec, sum_Prod_4, sum_dissp_4,
!sum_buoyF_4

! Calculating Phi_d and Phi_b2
       
         IF (PE_cal) THEN
        write(6,'(a)') 'Before call z_star', kk

        call z_star_cal

         write(6,'(a)') 'After call z_star', kk
        END IF
        end do

        sum_Prod_avg = sum_Prod_avg/time_sum
        sum_dissp_avg = sum_dissp_avg/time_sum
        sum_buoyF_avg = sum_buoyF_avg/time_sum


       WRITE(6,*)  sum_Prod_avg, sum_dissp_avg, sum_buoyF_avg

       write(69,114) sum_Prod_avg, sum_dissp_avg, sum_buoyF_avg

112    format(f12.6,3E18.9)
113    format(f12.6,9E18.9)
114    format(3f12.6)       
       close(70)
       close(67)
       close(69)
         END IF

!        stop
        end

         subroutine read_grid


        use ntypes
        use Grid
        use domain

        implicit none

         integer :: i,j,k, NZ_T, NY_T


!Local Variables
         integer           :: s_1
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! grid allocation
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

allocate (GX(0:NX+1), stat=s_1)
allocate (GY(0:NY+1), stat=s_1)
allocate (GZ(0:NZ+1), stat=s_1)
allocate (DX(0:NX+1), stat=s_1)
allocate (DY(0:NY+1), stat=s_1)
allocate (DZ(0:NZ+1), stat=s_1)
allocate (GXF(0:NX), stat=s_1)
allocate (GYF(0:NY+1), stat=s_1)
allocate (GZF(0:NZ+1), stat=s_1)
allocate (DXF(0:NX), stat=s_1)
allocate (DYF(0:NY+1), stat=s_1)
allocate (DZF(0:NZ+1), stat=s_1)

        NZM = NZ-1
        NYM = NY-1

        if (s_1.NE.0) then
          write(6,*) "Error Allocating Grid Variables"
         endif

         OPEN (30,file='zgrid.txt',form='formatted',status='old')
         READ (30,*) NZ_T
!C Check to make sure that grid file is the correct dimensions
         IF (NZ_T.ne.NZ) THEN
           WRITE(6,*) 'NZ, NZ_T',NZ,NZ_T
           STOP 'Error: zgrid.txt wrong dimensions'
         END IF
         DO K=1,NZ+1
           READ(30,*) GZ(k)
         END DO
         DO K=1,NZ
           READ(30,*) GZF(k)
         END DO
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
           write(400,*) DZF(K)
         END DO
         DZ(0)=DZ(1)
         DZF(NZ+1)=DZF(NZ)

         OPEN (30,file='./ygrid.txt',form='formatted',status='old')
         READ (30,*) NY_T
!C Check to make sure that grid file is the correct dimensions
         IF (NY_T.ne.NY) THEN
           WRITE(6,*) 'NY, NY_T',NY,NY_T
           STOP 'Error: ygrid.txt wrong dimensions'
         END IF
         DO J=1,NY+1
           READ(30,*) GY(j)
         END DO
         DO J=1,NY
           READ(30,*) GYF(j)
         END DO
         CLOSE(30)



!C Define ghost cells
           GYF(0)=2.d0*GYF(1)-GYF(2)
           GYF(NY+1)=2.d0*GYF(NY)-GYF(NYM)
           GY(0)=2.d0*GY(1)-GY(2)
!C Define grid spacing
         DO J=1,NY+1
           DY(J)=(GYF(J)-GYF(J-1))
         END DO
         DO J=1,NY
           DYF(J)=(GY(J+1)-GY(J))
         END DO
         DY(0)=DY(1)
         DYF(NY+1)=DYF(NY)

        return
        end subroutine         

         subroutine plane_parav_vel(kk)

         use ntypes
         use Grid
         use Domain
         use run_variable

        implicit none
             
       integer i,j,k,kk,N
       !FILENAMES
!       parameter(np1 = ni, np2 = nj-2, np3 = nk-2)
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss
        np1 = ni
        np2 = nj-2
        np3 = nk-1
 
       allocate(g1vtk(1:np1),g2vtk(1:np2),g3vtk(1:np3),     &
                var_1(1:3*np1*np2*np3),var_2(1:3*np1*np2*np3))

        g1vtk(:) = 0.0
        g2vtk(:) = 0.0
        g3vtk(:) = 0.0
        var_1(:) = 0.0
        var_2(:) = 0.0

       do k=1,np3
         g3vtk(k) = GZF(k)
       enddo

       do j=1,np2       
         g2vtk(j) = GYF(j)
       enddo

       do i=1,np1      
         g1vtk(i) = real(DX_s)*(i-1)
       enddo

!         DO N=1,N_TH
!         th_wh(k,j,N)=dble(CTHX(0,k,j,N)) +  TH_BAR(K,j,N) !
!         ENDDO


 
!      IF (Non_linear_ST) THEN
!          call density_TC(.true.,.true.)
!       ELSE 
!       do j=0,NY+1
!       do k=0,NZ+1
!         th_wh(k,j,N_TH+1)= -alpha_w*dble(CTHX(0,K,J,1)) +
!         gamma_w*dble(CTHX(0,K,J,2)) +  TH_BAR(K,j,N_TH+1) !
!       enddo
!       enddo
!       ENDIF

!       outputDIR='plane_3D_para/'
!       basename = 'data_parav'

       IF((CAL_STAT_AVG) .AND. (kk .eq. 1) ) THEN
        outputDIR='plane_data_3D/'

        basename = 'data_parav_avg'

        IF (IBM) THEN
         DO I=1,ni
         DO J=1,hill_ind_tot(I)+1
         u1_avg(:,J,I)=1.0d0
         u2_avg(:,J,I)=1.0d0
         u3_avg(:,J,I)=1.0d0
         th_wh_avg(:,J,I)=100.0d0
         ENDDO
         ENDDO
         ENDIF

        write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR),  &
             trim(basename)//".vtk"
       ELSE
        outputDIR='plane_data_3D/'
        basename = 'data_parav'

         IF (IBM) THEN
         DO I=1,ni
         DO J=1,hill_ind_tot(I)+1
         u(I,:,J)=1.0d0
         v(I,:,J)=1.0d0
         w(I,:,J)=1.0d0
         th_wh(:,J,I)=100.0d0
         ENDDO
         ENDDO
         ENDIF


        write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR),  &
             trim(basename)//"_n",kk,".vtk"
       ENDIF

       open(unit=13,file=OutFileName,access='stream', &
      form='unformatted',status='unknown',            &
      convert='big_endian',iostat=s11)




        WRITE(6,*) 'Wrting in ', OutFileName,' for paraview'
!HEADER: note termination with char(10)
        
        write(13) "# vtk DataFile Version 3.0"//char(10)
        write(13) trim(BaseName)//char(10)
        write(13) "BINARY"//char(10)
        write(13) "DATASET RECTILINEAR_GRID"//char(10)

        write(ss,fmt='(A10,3I5)') "DIMENSIONS",np1,np3,np2
        write(13) ss//char(10)
!Xgrid
        write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",np1," float"
        write(13) char(10)//ss//char(10)
        do i = 1,np1
        write(13) g1vtk(i)
        enddo
!Ygrid
        write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np3," float"
        write(13) char(10)//ss//char(10)
        do k = 1,np3
        write(13) g3vtk(k)
        enddo
!Zgrid
        write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np2," float"
        write(13) char(10)//ss//char(10)
        do j = 1,np2
        write(13) g2vtk(j)
        enddo
  
        
!Field
      IF((CAL_STAT_AVG) .AND. (kk .eq. 1) ) THEN 
        jj = 1
        do j = 1,np2
         do k = 1, np3
          do i = 1, np1
           var_1(jj) =  real(u1_avg(k,j,i))
           jj = jj+1
           var_1(jj) =  real(u2_avg(k,j,i))
           jj = jj+1
           var_1(jj) =  real(u3_avg(k,j,i))
           jj = jj+1
        end do 
          end do
            end do
         jj = 1
        do j = 1,np2
         do k = 1, np3
          do i = 1, np1
           var_2(jj) =  real(th_wh_avg(k,j,i))
           jj = jj+1 
         enddo
         enddo
        enddo
        ELSE
       jj = 1 
        do j = 1,np2
         do k = 1, np3
          do i = 1, np1
           var_1(jj) =  real(u(i,k,j))
           jj = jj+1
           var_1(jj) =  real(v(i,k,j))   
           jj = jj+1
           var_1(jj) =  real(w(i,k,j))
           jj = jj+1
         enddo 
        enddo
        enddo
      jj = 1
        do j = 1,np2
         do k = 1, np3
          do i = 1, np1
           var_2(jj) =  real(th_wh(k,j,i))
           jj = jj+1
         enddo
         enddo
        enddo  
         ENDIF
 
!Field
        IF((CAL_STAT_AVG) .AND. (kk .eq. 1) ) THEN
        write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2*np3
        write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors_avg float 1"//char(10)
        write(13) var_1
         write(13) "SCALARS TH_wh_avg float 1"//char(10)
         write(13) "LOOKUP_TABLE default"//char(10)
         write(13) var_2
        ELSE
        write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2*np3
        write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var_1
         write(13) "SCALARS TH_wh float 1"//char(10)
         write(13) "LOOKUP_TABLE default"//char(10)
         write(13) var_2
!         write(13) "SCALARS omega_x float 1"//char(10)
!         write(13) "LOOKUP_TABLE default"//char(10)
!         write(13) th_wh
        END IF

!Close VTK File
        close(13)


        deallocate(g1vtk,g2vtk,g3vtk,var_1,var_2)
!if (allocated(SPplane) ) deallocate(SPplane)
!if (allocated(DPplane) ) deallocate(DPplane)

        return
        end subroutine plane_parav_vel

        subroutine call_stat_cal_avg
         use ntypes
         use Grid
         use Domain
         use run_variable

         implicit none

         integer  :: k,j,p
         real(r8) :: temp_val

!       s1 =   Gravity*alpha_T*th*v

      sum_Prod  = 0.0d0 ;
      sum_dissp = 0.0d0 ;
      sum_buoyF = 0.0d0 ;       
      sum_dtkdt = 0.0d0 ;
      sum_vis   = 0.0d0 ;
      sum_trans = 0.0d0 ;
      sum_advec = 0.0d0 ;
      vol_box   = 0.0d0 ;
      sum_wstress = 0.0d0;

!      sum_advec_1 = 0.0d0 ;
!      sum_advec_2 = 0.0d0 ;
!      sum_advec_3 = 0.0d0 ;
!      sum_advec_4 = 0.0d0 ;

!      sum_Prod_1  = 0.0d0 ;
!      sum_Prod_2  = 0.0d0 ;
!      sum_Prod_3  = 0.0d0 ;
!      sum_Prod_4  = 0.0d0 ;

!      sum_buoyF_1  = 0.0d0 ;
!      sum_buoyF_2  = 0.0d0 ;
!      sum_buoyF_3  = 0.0d0 ;
!      sum_buoyF_4  = 0.0d0 ;

!      sum_dissp_1  = 0.0d0 ;
!      sum_dissp_2  = 0.0d0 ;
!      sum_dissp_3  = 0.0d0 ;
!      sum_dissp_4  = 0.0d0 ;

      do p=2,np-1
      do j=2,nj-1
      do k=2,nk-1
        dv        =  DZF(k-1)*DYF(j-1)*DX_s
        vol_box   = vol_box  + dv

!      calculation of advection terms
 
         temp_val = - u1_avg(k,j,p)*(s1(k+1,j,p)-s1(k-1,j,p))/&
                        (2.0d0*DZF(k-1))&
!     u*d<E>/dx
                    - u2_avg(k,j,p)*(s1(k,j+1,p)-s1(k,j-1,p))/&
                        (2.0d0*DYF(j-1))&
!     v*d<E>/dy
                    - u3_avg(k,j,p)*(s1(k,j,p+1)-s1(k,j,p-1))/&
                        (2.0d0*DX_s)
!     w*d<E>/dz

       sum_advec =sum_advec + temp_val*dv 

!         if ( (i > i_1) .and. (j <= j_1) ) then
!         sum_advec_1 = sum_advec_1 +  temp_val*dv
!        endif
!        if ( (i <= i_2)  ) then
!         sum_advec_2 = sum_advec_2 +  temp_val*dv
!        endif
!        if ( (i > i_3) .and. (j > j_3) ) then
!         sum_advec_3 = sum_advec_3 +  temp_val*dv
!        endif
!        if ( (i > i_1) .and. ( (j <= j_3) .and. (j > j_1)) ) then
!         sum_advec_4 = sum_advec_4 +  temp_val*dv
!        endif

!      calculation of buoyancy flux terms

        temp_val = Gravity*alpha_T*th_wh_avg(k,j,p)*(u2_avg(k,j+1,p)&
                   +u2_avg(k,j,p))/2.0d0
        sum_buoyF = sum_buoyF +  temp_val*dv

!        if ( (i > i_1) .and. (j <= j_1) ) then
!         sum_buoyF_1 = sum_buoyF_1 +  temp_val*dv
!        endif
!        if ( (i <= i_2)  ) then
!         sum_buoyF_2 = sum_buoyF_2 +  temp_val*dv
!        endif
!        if ( (i > i_3) .and. (j > j_3) ) then
!         sum_buoyF_3 = sum_buoyF_3 +  temp_val*dv
!        endif
!        if ( (i > i_1) .and. ( (j <= j_3) .and. (j > j_1)) ) then
!         sum_buoyF_4 = sum_buoyF_4 +  temp_val*dv
!        endif

!      calculation of dissipation terms

        temp_val  =   ((u1_avg(k+1,j,p)-u1_avg(k,j,p))/DZF(k-1))**2.0+ &
                     (0.50*(u1_avg(k+1,j+1,p)+u1_avg(k,j+1,p) &
                        -u1_avg(k+1,j-1,p)-&
                     u1_avg(k,j-1,p))/(2.0d0*DYF(j-1)))**2.0 + &
                    ((u1_avg(k,j,p+1)-u1_avg(k,j,p-1))/(2.0*DX_s))**2.0&
                        + &
!                    Done with (dudx_i)^2.0
                    (0.50*(u2_avg(k+1,j+1,p)+u2_avg(k+1,j,p)&
                        -u2_avg(k-1,j+1,p)-&
                        u2_avg(k-1,j,p))/(2.0d0*DZF(k-1)))**2.0+ &
                     ((u2_avg(k,j+1,p)-u2_avg(k,j,p))/DYF(j-1))**2.0 + &
                     ((u2_avg(k,j,p+1)-u2_avg(k,j,p-1))/&
                        (2.0*DX_s))**2.0 &
!                    Done with (dvdx_i)^2.0 
                    + ((u3_avg(k+1,j,p)-u3_avg(k-1,j,p))/&
                        (2.0*DZF(k-1)))**2.0 &
                    + ((u3_avg(k,j+1,p)-u3_avg(k,j-1,p))/&
                        (2.0*DYF(j-1)))**2.0 &
                    + ((u3_avg(k,j,p+1)-u3_avg(k,j,p-1))/(2.0*DX_s))**2.0 
!                    Done with (dwdx_i)^2.0

        sum_dissp  = sum_dissp - NU*temp_val*dv
      
!        if ( (i > i_1) .and. (j <= j_1) ) then
!         sum_dissp_1  = sum_dissp_1 - NU*temp_val*dv
!       endif
!       if ( i <= i_2)  then
!         sum_dissp_2  = sum_dissp_2 - NU*temp_val*dv
!       endif
!       if ( (i > i_3) .and. (j > j_3) ) then
!         sum_dissp_3  = sum_dissp_3 - NU*temp_val*dv
!       endif
!       if ( (i > i_1) .and. ( (j <= j_3) .and. (j > j_1)) ) then
!         sum_dissp_4  = sum_dissp_4 - NU*temp_val*dv
!       endif
     
      enddo
      enddo
      enddo 

      do p = 2,np-1
        do k = 2,nk-1

        dS = DZF(k-1)*DX_s

! calculation of KE from wind stress

      temp_val =((u3_avg(k,nj-1,p)+u3_avg(k,nj-2,p))/2.0)*&
                ((u3_avg(k,nj-1,p)-u3_avg(k,nj-2,p))/(2.0*DYF(nj-1)))
!                    Done with (u*du/dx_j)
        sum_wstress  = sum_wstress + NU*temp_val*(dS)

        end do
        end do

      sum_advec  = sum_advec*rho_0
      sum_buoyF  = sum_buoyF*rho_0
      sum_dissp  = sum_dissp*rho_0
      sum_wstress = sum_wstress*rho_0

!      sum_advec_1  = sum_advec_1*rho_0 ;
!      sum_advec_2  = sum_advec_2*rho_0 ;
!      sum_advec_3  = sum_advec_3*rho_0 ;
!      sum_advec_4  = sum_advec_4*rho_0 ;

!      sum_buoyF_1  = sum_buoyF_1*rho_0 ;
!      sum_buoyF_2  = sum_buoyF_2*rho_0 ;
!      sum_buoyF_3  = sum_buoyF_3*rho_0 ;
!      sum_buoyF_4  = sum_buoyF_4*rho_0 ;

!      sum_dissp_1  = sum_dissp_1*rho_0 ;
!      sum_dissp_2  = sum_dissp_2*rho_0 ;
!      sum_dissp_3  = sum_dissp_3*rho_0 ;
!      sum_dissp_4  = sum_dissp_4*rho_0 ;

!      sum_buoyF  = sum_buoyF/vol_box
!      sum_dissp  = sum_dissp/vol_box

return
        end subroutine


        subroutine call_stat_cal
         use ntypes
         use Grid
         use Domain
         use run_variable

         implicit none

         integer  :: k,j,p
         real(r8) :: temp_val

 

          u(:,:,:)  = 0.0d0 
          w(:,:,:)  = 0.0d0
          v(:,:,:)  = 0.0d0
          th(:,:,:) = 0.0d0 

         u = u1-u1_avg
         v = u2-u2_avg
         w = u3-u3_avg
         th= th_wh-th_wh_avg

         !s1 =   Gravity*alpha_T*th*v

      sum_Prod  = 0.0d0 ;
      sum_dissp = 0.0d0 ;
      sum_buoyF = 0.0d0 ;       
      sum_dtkdt = 0.0d0 ;
      sum_vis   = 0.0d0 ;
      sum_trans = 0.0d0 ;
      sum_advec = 0.0d0 ;
      vol_box   = 0.0d0 ;

!      sum_Prod_1  = 0.0d0 ;
!      sum_Prod_2  = 0.0d0 ;
!      sum_Prod_3  = 0.0d0 ;
!      sum_Prod_4  = 0.0d0 ;

!      sum_buoyF_1  = 0.0d0 ;
!      sum_buoyF_2  = 0.0d0 ;
!      sum_buoyF_3  = 0.0d0 ;
!      sum_buoyF_4  = 0.0d0 ;

!      sum_dissp_1  = 0.0d0 ;
!      sum_dissp_2  = 0.0d0 ;
!      sum_dissp_3  = 0.0d0 ;
!      sum_dissp_4  = 0.0d0 ;

      do p=2,np-1
      do j=2,nj-1
      do k=2,nk-1
        dv        =  DZF(k-1)*DYF(j-1)*DX_s
        vol_box   = vol_box  + dv
        s1(k,j,p) = Gravity*alpha_T*th(k,j,p)*(v(k,j+1,p)+v(k,j,p))/2.0d0

        
        sum_buoyF = sum_buoyF +  s1(k,j,p)*dv

        !if ( (i > i_1) .and. (j <= j_1) ) then
        ! sum_buoyF_1 = sum_buoyF_1 +  s1(i,j,p)*dv
        !endif
        !if ( (i <= i_2)  ) then
        ! sum_buoyF_2 = sum_buoyF_2 +  s1(i,j,p)*dv
        !endif
        !if ( (i > i_3) .and. (j > j_3) ) then
        ! sum_buoyF_3 = sum_buoyF_3 +  s1(i,j,p)*dv
        !endif
        !if ( (i > i_1) .and. ( (j <= j_3) .and. (j > j_1)) ) then
        ! sum_buoyF_4 = sum_buoyF_4 +  s1(i,j,p)*dv
        !endif




        temp_val = u(k,j,p)*u(k,j,p)*(u1_avg(k+1,j,p)-u1_avg(k,j,p))/&
                        DZF(k-1) &
                 + 0.50d0*0.50d0*0.50d0*(v(k,j+1,p)+v(k,j,p))* &
        (u(k+1,j,p)+u(k,j,p))*(u2_avg(k+1,j+1,p)+u2_avg(k+1,j,p)-&
        u2_avg(k-1,j+1,p)-u2_avg(k-1,j,p))/(2.0d0*DZF(k-1)) &
                 + &
        0.50d0*w(k,j,p)*(u(k+1,j,p)+u(k,j,p))*(u3_avg(k+1,j,p) &
        -u3_avg(k-1,j,p))/(2.0*DZF(k-1)) &                 
!           Done with (.)ud()/dx
                 + &
        0.50d0*0.50d0*0.50d0*(v(k,j+1,p)+v(k,j,p))*(u(k+1,j,p) &
        +u(k,j,p))*(u1_avg(k+1,j+1,p)+u1_avg(k,j+1,p)-u1_avg(k+1,j-1,p)&
        -u1_avg(k,j-1,p))/(2.0d0*DYF(j-1)) &
                 + &
        0.50d0*0.50d0*(v(k,j+1,p)+v(k,j,p))*(v(k,j+1,p)+v(k,j,p))&
        *(u2_avg(k,j+1,p)-u2_avg(k,j,p))/DYF(j-1) &
                 + &
        0.50d0*w(k,j,p)*(v(k,j+1,p)+v(k,j,p))*(u3_avg(k,j+1,p)&
        -u3_avg(k,j-1,p))/(2.0*DYF(j-1))
!           Done with (.)vd()/dy
        sum_Prod  = sum_Prod - temp_val*dv

       ! if ( (i > i_1) .and. (j <= j_1) ) then
       !  sum_Prod_1  = sum_Prod_1 - temp_val*dv
       ! endif
       ! if ( i <= i_2 ) then
       !  sum_Prod_2  = sum_Prod_2 - temp_val*dv
       ! endif
       ! if ( (i > i_3) .and. (j > j_3) ) then
       !  sum_Prod_3  = sum_Prod_3 - temp_val*dv
       ! endif
       ! if ( (i > i_1) .and. ( (j <= j_3) .and. (j > j_1)) ) then
       !  sum_Prod_4  = sum_Prod_4 - temp_val*dv
       ! endif



        temp_val  =   ((u(k+1,j,p)-u(k,j,p))/DZF(k-1))**2.0 &
                    + &
        (0.50*(u(k+1,j+1,p)+u(k,j+1,p)-u(k+1,j-1,p)-u(k,j-1,p))/&
                (2.0d0*DYF(j-1)))**2.0 &
                    + ((u(k,j,p+1)-u(k,j,p-1))/(2.0*DX_s))**2.0 &   
!                    Done with (dudx_i)^2.0
                    + &
        (0.50*(v(k+1,j+1,p)+v(k+1,j,p)-v(k-1,j+1,p)-v(k-1,j,p))/&
        (2.0d0*DZF(k-1)))**2.0 &
                    + ((v(k,j+1,p)-v(k,j,p))/DYF(j-1))**2.0 &
                    + ((v(k,j,p+1)-v(k,j,p-1))/(2.0*DX_s))**2.0 &
!                    Done with (dvdx_i)^2.0 
                    + ((w(k+1,j,p)-w(k-1,j,p))/(2.0*DZF(k-1)))**2.0 &
                    + ((w(k,j+1,p)-w(k,j-1,p))/(2.0*DYF(j-1)))**2.0 &
                    + ((w(k,j,p+1)-w(k,j,p-1))/(2.0*DX_s))**2.0
!                    Done with (dwdx_i)^2.0

       sum_dissp  = sum_dissp - NU*temp_val*dv

      ! if ( (i > i_1) .and. (j <= j_1) ) then
      !   sum_dissp_1  = sum_dissp_1 - NU*temp_val*dv
      ! endif
      ! if ( i .lt. i_2)  then
      !   sum_dissp_2  = sum_dissp_2 - NU*temp_val*dv
      ! endif
      ! if ( (i > i_3) .and. (j > j_3) ) then
      !   sum_dissp_3  = sum_dissp_3 - NU*temp_val*dv
      ! endif
      ! if ( (i > i_1) .and. ( (j <= j_3) .and. (j > j_1)) ) then
      !   sum_dissp_4  = sum_dissp_4 - NU*temp_val*dv
      ! endif

      enddo
      enddo
      enddo 

      sum_buoyF  = sum_buoyF*rho_0
      sum_dissp  = sum_dissp*rho_0
      sum_Prod   = sum_Prod*rho_0

      !sum_Prod_1  = sum_Prod_1*rho_0 ;
      !sum_Prod_2  = sum_Prod_2*rho_0 ;
      !sum_Prod_3  = sum_Prod_3*rho_0 ;
      !sum_Prod_4  = sum_Prod_4*rho_0 ;

      !sum_buoyF_1  = sum_buoyF_1*rho_0 ;
      !sum_buoyF_2  = sum_buoyF_2*rho_0 ;
      !sum_buoyF_3  = sum_buoyF_3*rho_0 ;
      !sum_buoyF_4  = sum_buoyF_4*rho_0 ;

      !sum_dissp_1  = sum_dissp_1*rho_0 ;
      !sum_dissp_2  = sum_dissp_2*rho_0 ;
      !sum_dissp_3  = sum_dissp_3*rho_0 ;
      !sum_dissp_4  = sum_dissp_4*rho_0 ;


!     sum_buoyF  = sum_buoyF !/vol_box
!     sum_Prod   = sum_Prod  !/vol_box
!     sum_dissp  = sum_dissp !/vol_box
    
        return
        end subroutine
 
        subroutine z_star_cal
        use Grid
        use Domain
        use run_variable
        use qsort_c_module
        use ntypes
        implicit none

      integer   ::   i, j, k, kk, tot_ind
      real(r8)  ::   area, rho_min, rho_max, sum_vol, z_star_pre, BIG_NUM
!     real(r8)  ::   DX_s,DY_s,GZF(0:nz+1)
!     real(r8)  ::   u_0,Gravity,alpha_T, NU, rho_0
      real(r8)  ::   sum_val, temp_val, dv_z
      real(r8),allocatable,dimension(:,:,:) ::  S1S
!     real(r8),allocatable,dimension(:,:)   ::  rho_1D_tot
      real(r8),allocatable,dimension(:,:)   ::  z_start_bar_surf
      real(r8),allocatable,dimension(:)     ::  drho_1D_totdz_s
      real(r8),allocatable,dimension(:)     ::  drho_1D_totdz_s_check
      real(r8),allocatable,dimension(:)     ::  z_start_bar,phi_d_prof
      real(r8),allocatable,dimension(:)     ::  phi_b2_prof,phi_d_total

      parameter (BIG_NUM = 10.d0**8.0)

      area = DX_s*(i-1)*GZF(nk-1)
      tot_ind = ni*nj*nk

!!Start  z* calculation
       allocate (S1S(1:ni,1:nk,1:nj))
       allocate (rho_1D_tot(1:tot_ind,1:2))
       allocate (varSP_zstar(1:ni,1:nk,1:nj))
       allocate (z_start_bar_surf(1:ni,1:nk))
        allocate (z_start_bar(1:nj))

       drho_1D_totdz_s_check(:) = 0.0d0
       drho_1D_totdz_s(:)       = 0.0d0
       rho_1D_tot(:,:)          = 0.0d0
       varSP_zstar(:,:,:)       = 0.0d0
       S1S(:,:,:)               = 0.0d0
       z_start_bar_surf(:,:)    = 0.0d0
       z_start_bar(:)           = 0.0d0
!      rho_1D_tot(:,1) = 100.d0 ;
!      rho_1D_tot(:,2) = -10.d0 ;

!write(6,*) 'shape varSP_zstar', shape(varSP_zstar)
!write(6,*) 'nxp2, nyp2, nzp2', nxp2, nyp2, nzp2
!write(6,*) 'size rho_1D_tot', shape(rho_1D_tot)


      jj = 0
       do i=1,ni
         do k=1,nk
           do j=1,nj
           jj =  jj + 1
           ii =  nj*nk*(i-1) + nj*(k-1) + j           ! Calculate unique
                                                                !grid index
           rho_1D_tot(jj,1) = -alpha_T*th_wh(i,k,j)+1.0d0 !where T_0=10C
           ! rho_1D_tot(jj,1) = -alpha_T*varSP_zstar(i,j,k) !+ 1.0d0+
           ! THBAR(j,N_TH+1)
           ! Storing the rho in a 1D array
           rho_1D_tot(jj,2) =  dble(ii)  ! dble(tot_ind*rank + jj)             
           ! Storing the index for each grid to a 1D array
           enddo
         enddo
       enddo

!write(6,*) 'shape rho_1D_tot = ', shape(rho_1D_tot)

     !      sorting for rho values : keep the grid index
        call QsortC(rho_1D_tot(:,1),rho_1D_tot(:,2))

! open(unit=56,file="rho_1D_tot.dat",form='formatted',status='unknown')
!do ii = 1,100
! write(56,'(2es20.10)') rho_1D_tot(ii,1), rho_1D_tot(ii,2)
!enddo

!do ii = 1,100
! write(56,'(2es20.10)') rho_1D_tot(jj-ii,1), rho_1D_tot(jj-ii,2)
!enddo

!open(unit=56,file="rho_1D_tot.dat",form='formatted',status='unknown')

        kk = tot_ind
        DO ii = 1,tot_ind
         drho_1D_totdz_s_check(ii) = rho_1D_tot(ii,1)
!!!!!!!!!!!!!!!!!!!!!!!1

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

         DO ii = kk, 1,-1 
!           NP_I = floor(rho_1D_tot(ii,2)/tot_ind)
              i = floor((rho_1D_tot(ii,2))/(nk*nj)) 
              k = floor((rho_1D_tot(ii,2) - i*(nk*nj))/nj) 
              j = floor(rho_1D_tot(ii,2) - i*(nk*nj) - k*nj) + 1

              DV  = DX_s*DZ_s*abs(GYF(j+1)-GYF(j))
!             calculating z* at different ii level and  mapped back to
!             physical grid(i,j,k)
              rho_1D_tot(int(rho_1D_tot(ii,2)),1) = (sum_vol + &
                                                         0.50d0*DV)/area
              sum_vol =  sum_vol +  DV  
 
!             calculating dz*drho at different ii level : not yet mapped
!             to physical grid(i,j,k) 
              drho_1D_totdz_s(ii) = drho_1D_totdz_s(ii)*&
                                (z_star_pre-(sum_vol + 0.50d0*DV)/area)

              z_star_pre =   (sum_vol + 0.50d0*DV)/area  !  storing the z* for (ii) level to calculate z*(ii)-z*(ii-1)
!              write(66,666)dble(j),dble(ii),drho_1D_totdz_s_check(ii),(sum_vol
!              + 0.50d0*DV)/area, drho_1D_totdz_s(ii)
         ENDDO
!         close(66)

!        interpolating dz*drhp at last point of the array (here is kk)
         drho_1D_totdz_s(kk) = 2.0d0*drho_1D_totdz_s(kk-1) &
                                                - drho_1D_totdz_s(kk-2)

         jj = 0
         do i=1,ni
          do k=1,nk
           do j=1,nj
              jj =  jj + 1              
              S1S(I,K,J) = rho_1D_tot(jj,1)   
           enddo
          enddo
         enddo


        do j=1,nj
         z_start_bar(j) = z_start_bar(j) + &
                           sum(S1S(1:ni,1:nk,nj))/dble(ni*nk)*delta_time
        enddo 
        
        do i=1,ni
          do k=1,nk
            z_start_bar_surf(i,k) = z_start_bar_surf(i,k) + &
                                         S1S(i,k,nj)*delta_time
          enddo
        enddo
         
!        phi_d_prof(:) = 0.0d0
!        phi_b2_prof(:) = 0.0d0

        do j=3,nj-2
        do k=2,nk-1
        do i=2,ni-1
!      calculation of buoyancy flux terms

        dv = DX_s*DZF(k-1)*(GYF(j+1)-GYF(j))
        dv_z= DX_s*DZF(k-1)
        sum_val = sum_val + dv

            temp_val = alpha_T*((S1S(i+1,k,j)-S1S(i-1,k,j))/(2.0*DX_s))&
                      *((TH_WH(i+1,k,j)-TH_WH(i-1,k,j))/(2.0*DX_s)) &
                      +alpha_T*((S1S(i,k,j+1)-S1S(i,k,j-1))/(GYF(j+1)&
                      -GYF(j-1)))*((TH_WH(i,k,j+1)-TH_WH(i,k,j-1)) &
                      /(GYF(j+1)-GYF(j-1)))  &                      
                      +alpha_T*((S1S(i,k+1,j)-S1S(i,k-1,j))/(2.0*DZF(k-1)))&
                      *((TH_WH(i,k+1,j)-TH_WH(i,k-1,j))/(2.0*DZF(k-1)))

            phi_d_prof(J) = phi_d_prof(J) &
                            + Gravity*NU*temp_val*dv_z*delta_time 

            temp_val  = alpha_T*(S1S(i,k,j+1)*((TH_WH(i,k,j+2)&
                                -TH_WH(i,k,j))/(GYF(j+2)-GYF(j))) &
                        - S1S(i,k,j-1)*((TH_WH(i,k,j)-TH_WH(i,k,j-2))&
                                /(GYF(j)-GYF(j-2))))/(GYF(j+1)-GYF(j-1))

            phi_b2_prof(J) = phi_b2_prof(J)&
                                -Gravity*NU*temp_val*dv_z*delta_time

          enddo
          enddo
          enddo
        
           phi_b2_prof(2) = 2.0d0*phi_b2_prof(3) - phi_b2_prof(4)
           phi_b2_prof(1) = 2.0d0*phi_b2_prof(2) - phi_b2_prof(3)

           phi_d_prof(2) = 2.0d0*phi_d_prof(3) - phi_d_prof(4)
           phi_d_prof(1) = 2.0d0*phi_d_prof(2) - phi_d_prof(3) 

           temp_val = 0.0d0

        do k=1,nk
        do i=1,ni
           temp_val = temp_val + alpha_T*Gravity*NU*S1(i,k,1)*&
                      ((TH_WH(i,k,2)-TH_WH(i,k,1))/(GYF(2)-GYF(1)))*dv_z
        enddo
        enddo  

        write(6,*) temp_val 
!           phi_d = temp_val
!666    format(5E17.8)


     deallocate(S1S,rho_1D_tot, drho_1D_totdz_s, drho_1D_totdz_s_check)
     deallocate(varSP_zstar, z_start_bar_surf, z_start_bar)

     return
             end

