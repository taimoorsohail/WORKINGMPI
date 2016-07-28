      Program main
!      INCLUDE 'header_fft'
      implicit none
      INTEGER   NX, NY, NZ, N_TH

      PARAMETER (NX=512)
      PARAMETER (NY=129)
      PARAMETER (NZ=257)
      PARAMETER (N_TH=1)

      integer ni,nj,nk,nk_st,n_t, st,en,tstep,N,ii
      integer nip,njp,NXP,NXV,NXP_L,NX_PAR,rank,NP
      parameter (ni=NX, nj=NZ+2,nk_st=1102,nk=2000,NP=257)
      parameter (n_t = nk-nk_st-1)
      PARAMETER (NXV=ceiling(real(NX+2)/real(NP))*NP)
      PARAMETER (NXP=NXV/NP-1)
      Real*8 x(ni,nj),y(ni,nj)
      Real*8 u(ni,nj),v(ni,nj),w(ni,nj),p(ni,nj),th(ni,nj)
!     &       th_wh(ni,nj)

      real*8   RI,U_0,omega,rho_0,nu,dt,dx,xloc,TIME,PI
      PARAMETER( RI = 14.93d0, U_0 = 0.0395d0, omega = 1.d0, &
                rho_0 = 1.d0, nu=1.0d-7 ) 
      integer  i,j,imax,jmax,kmax,k,kk,i_md,ind_lt,ind_rt
      integer  zindex,zindex2,xindex,xindex2,index_Mfac
      parameter (zindex=201,zindex2=240,xindex=341,xindex2=107)
   
       integer np1,np2
       !FILENAMES
       parameter(np1 = ni, np2 = nj)
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss

 !      real(r4),allocatable,dimension(:) :: g1vtk,g2vtk
       real(4) :: g1vtk(1:np1),g2vtk(1:np2), var_1(1:3*np1*np2)
       real(4) :: th_wh(np1,np2)
!       th_wh(1:np1,1:np2,N_TH+1)
!       real(r4),allocatable,dimension(:,:) :: var1 
!       real(r8),allocatable,dimension(:,:,:) :: th_wh
       !STATUS VARIABLES
       integer :: s11
   
!      REAL*8     U_data(0:n_t+1,0:NZ+1,0:NY+1), t(0:n_t+1)	!(nk_st:nk)
!      COMPLEX*16 CU_data(0:n_t/2,0:NZ+1,0:NY+1)
!      EQUIVALENCE (U_data,CU_data)

      integer  debug,ier,itot
      integer  tecini,tecdat,teczne,tecnod,tecfil,tecend
      integer  visdouble,disdouble
      character*1 nulchar

      CHARACTER*4 PID
      CHARACTER*5 file_num
      CHARACTER*35 file_name_1,file_time
      CHARACTER*8 base_name
 	  character(len=18) :: title

      LOGICAL VEL_FIELD, VEL_FIELD_RAW
!       PARAMETER (NKXV=ceiling(real(NKX+1)/real(NP))*NP) 
!       PARAMETER (NKXP=NKXV/NP-1)
!       PARAMETER (NX2V=ceiling(real(NXV/2)/real(NP))*NP)
!       PARAMETER (NX2P=NX2V/NP-1)    

      base_name = 'data_ZX_'
      print *, 'folder name=', base_name
      
      PI  = ATAN(1.0)*4.0

      write(6,*) 'U_0', U_0, 'RI', RI

      VEL_FIELD = .false.
      VEL_FIELD_RAW = .true.
      

      nulchar = char(0)
      debug   = 0
      visdouble = 0
      disdouble = 1
      imax = ni
      jmax = nj
      kmax = 1        
      kk = 0

       
      kk = -1
       
	DO k=nk_st,nk
	  kk= k !kk+1
		
      file_num = CHAR(MOD(k,100000)/10000+48) &
                //CHAR(MOD(k,10000)/1000+48)  &
     		//CHAR(MOD(k,1000)/100+48)    & 
     		//CHAR(MOD(k,100)/10+48)      & 
     		//CHAR(MOD(k,10)+48)
       
     	st=1
     
	 DO rank=0,NP-1
	   I=rank+1
       PID = '_'//CHAR(MOD(I,1000)/100+48) &  
     		//CHAR(MOD(I,100)/10+48)   & 
     		//CHAR(MOD(I,10)+48)
       
       file_name_1 = 'plane_data_ZX/'//base_name//file_num//PID//'.pln'

	   NXP_L = NX-(NXP+1)*rank-1
	   NX_PAR = min(NXP,NXP_L)
	   if(rank.EQ.0)  write(6,'(2i6,2a)') k,kk,' ',file_name_1

       open(22,file=file_name_1,status='old',form='unformatted')
!        print *, file_name_1
       read(22) time, dt, tstep, nip, njp
!        if(rank.EQ.0) then
!        	write(*,'(2f8.4,i6,f8.4,2i6)') time, dt, tstep, xloc, nip, njp
!        endif
       if(njp.NE.nj)	then 
       	print *, '!!! WARNING !!! Check NJ', nj,'NJP', njp, tstep
       	STOP
       endif
       en=st+nip-1
       read(22) y(1,:),dx,u(st:en,:),         &
                       w(st:en,:),v(st:en,:), &
                       th(st:en,:)
       close(22)
       st=en+1
	 ENDDO
	 x(1,:) = 0.d0
	 DO i=2,ni
	 	y(i,:) = y(1,:)
	 	x(i,:) = x(i-1,:)+dx
	 ENDDO
! 	 print *, 'dx ', dx
! 	 print *, 'End index:', en, 'nx:', ni
 
!       read(20,*) time, dt,u_bulk,v_bulk,u_bc,v_bc 
!       write(6,*)k,kk,dt,time
      
!       t(kk) = time


 !     do i=1,ni       
 !      do j=nj,1,-1
 !       U_data(kk-1,i-1,j-1)   = u(i,j) 
 !      end do
 !     end do	! end of loop over 'i' index
!_____________________________________________________________________       
       
      IF(VEL_FIELD_RAW) THEN

       outputDIR='plane_data_zx_raw/'
       basename = 'data_parav_zx'


       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR), &
        trim(basename)//"_n",kk,".pln"

       open(unit=33,file=OutFileName,access='stream',  &
       form='unformatted',status='unknown',            &
       convert='big_endian',iostat=s11)


       WRITE(6,*) 'Wrting in ', OutFileName
       write(33) y(1,:),dx,u(:,:),         &
                       w(:,:),v(:,:), &
                       th(:,:)
       close(33)

      ENDIF




      IF(VEL_FIELD) THEN

       do i=1,np1
         g1vtk(i) = x(i,1)
       enddo

       do j=1,np2
         g2vtk(j) = y(1,j)
       enddo



       outputDIR='plane_data_zx/'
       basename = 'data_parav_zx'


       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR), &
        trim(basename)//"_n",kk,".vtk"

       open(unit=13,file=OutFileName,access='stream',  &
       form='unformatted',status='unknown',            & 
       convert='big_endian',iostat=s11)


        WRITE(6,*) 'Wrting in ', OutFileName,' for paraview'


      !
      !
!HEADER: note termination with char(10)
       
        write(13) "# vtk DataFile Version 3.0"//char(10)
        write(13) trim(BaseName)//char(10)
        write(13) "BINARY"//char(10)
        write(13) "DATASET RECTILINEAR_GRID"//char(10)

        write(ss,fmt='(A10,3I5)') "DIMENSIONS",np1,np2,1
        write(13) ss//char(10)
!Xgrid
        write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",np1," float"
        write(13) char(10)//ss//char(10)
        do i = 1,np1
        write(13) g1vtk(i)
        enddo
!Ygrid
        write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np2," float"
        write(13) char(10)//ss//char(10)
        do j = 1,np2
          write(13) g2vtk(j)
        enddo
!Zgrid
        write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",1," float"
        write(13) char(10)//ss//char(10)
        write(13) 0.2

!Field
        ii = 1
        do j = 1, np2
        do i = 1, np1
        var_1(ii) =  real(v(i,j))
        ii = ii+1
        var_1(ii) =  real(u(i,j))
        ii = ii+1
        var_1(ii) =  real(w(i,j))
        ii = ii+1
         DO N=1,N_TH
         th_wh(i,j)= real(TH(i,j)) + 30.0 !+  THBAR(j,N) !
         ENDDO
        enddo
        enddo

!        do j = 1, np2
!        do i = 1, np1
!         th_wh(i,j,N_TH+1)=-alpha_T*THX(int(NXP/2),i,j,1)+THBAR(j,N_TH+1)
!        enddo
!        enddo


!Field
        write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
        write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var_1
!        write(13) "SCALARS Temp_div float 1"//char(10)
!        write(13) "LOOKUP_TABLE default"//char(10)
!        write(13) real(dble(THX(int(NXP/2),1:np1,1:np2,1)))
        write(13) "SCALARS Temp_wh float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(th_wh(1:np1,1:np2))
!        write(13) "SCALARS Density_wh float 1"//char(10)
!        write(13) "LOOKUP_TABLE default"//char(10)
!        write(13) real(th_wh(1:np1,1:np2,2))


!Close VTK File
        close(13)
      ENDIF
        
	ENDDO		! END OF LOOPING OVER FILES (i.e. end of TIME loop)
!_____________________________________________________________________       
!_____________________________________________________________________       

      write(6,*)'Total steps',  kk

        
        end

