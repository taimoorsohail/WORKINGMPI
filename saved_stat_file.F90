     subroutine plane_parav_plane_yz(kk)

!      INCLUDE 'header_duct'
 use ntypes
 use Domain
 use Grid, only : dx,GZF, GYF, GXF
 use run_variable
 use variable_stat
 implicit none
             
       integer i,j,k,np1,np2,kk,N
       !FILENAMES
       parameter(np1 = NZ, np2 = NY)
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss
       
 !      real(r4),allocatable,dimension(:) :: g1vtk,g2vtk
       real(r4) :: g1vtk(1:np1),g2vtk(1:np2), var_1(1:3*np1*np2)

!       th_wh(1:np1,1:np2,N_TH+1)
!       real(r4),allocatable,dimension(:,:) :: var1 
!       real(r8),allocatable,dimension(:,:,:) :: th_wh
       !STATUS VARIABLES
       integer :: s11



       do k=1,np1
         g1vtk(k) = GZF(k)
       enddo

       do j=1,np2       
         g2vtk(j) = GYF(J)
       enddo



       outputDIR='plane_data_yz/'
       basename = 'data_parav_yz'


       write(OutFileName,'(a,a,i4.4,a2,i5.5,a4)') trim(outputDIR), &
       trim(basename)//"_p",ipn,"_n",kk,".vtk"

       open(unit=13,file=OutFileName,access='stream',  &
      form='unformatted',status='unknown',             &
      convert='big_endian',iostat=s11)


        WRITE(6,*) 'Wrting in ', OutFileName,' for paraview for plane ipn=', ipn
!HEADER: note termination with char(10)
	
	write(13) "# vtk DataFile Version 3.0"//char(10)
	write(13) trim(BaseName)//char(10)
	write(13) "BINARY"//char(10)
	write(13) "DATASET RECTILINEAR_GRID"//char(10)

	write(ss,fmt='(A10,3I5)') "DIMENSIONS",1,np1,np2
	write(13) ss//char(10)
!Xgrid
	write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",1," float"
	write(13) char(10)//ss//char(10)
	write(13) real(GXF(ipn))
!Ygrid
	write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np1," float"
	write(13) char(10)//ss//char(10)
	do k = 1,np1
	write(13) g1vtk(k)
	enddo
!Zgrid
	write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np2," float"
	write(13) char(10)//ss//char(10)
	do j = 1,np2
	  write(13) g2vtk(j)
	enddo

!Field
        k = 1
        do j = 1, np2
        do i = 1, np1
        var_1(k) =  real(U1X(ipn,i,j))
        k = k+1
        var_1(k) =  real(U3X(ipn,i,j))
        k = k+1
        var_1(k) =  real(U2X(ipn,i,j))
        k = k+1
         DO N=1,N_TH
         th_wh(i,j,N)= THX(ipn,i,j,N) +  THBAR(j,N) !
         ENDDO
        enddo
        enddo
 
        do j = 1, np2
        do i = 1, np1
         th_wh(i,j,N_TH+1)=-alpha_T*THX(ipn,i,j,1)+THBAR(j,N_TH+1) 
        enddo
        enddo
     
        
        
!Field
	write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
	write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var_1
        write(13) "SCALARS Temp_div float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(dble(THX(ipn,1:np1,1:np2,1)))
        write(13) "SCALARS Temp_wh float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(th_wh(1:np1,1:np2,1))
        write(13) "SCALARS Density_wh float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(th_wh(1:np1,1:np2,2))

        WRITE(6,*) 'Finished Wrting in ', OutFileName,' for paraview'
        
!Close VTK File
	close(13)


        return
        end subroutine plane_parav_plane_yz
      
          subroutine plane_parav_plane_zx(kk)

!      INCLUDE 'header_duct'
 use ntypes
 use Domain
 use Grid, only : dx,GZF, GYF
 use run_variable
 use variable_stat
 use mpi_var, only: rank
 implicit none

       integer i,j,k,np1,np3,kk,N
       !FILENAMES
       parameter(np1 = NX, np3 = NZ)
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss

       real(r8),allocatable,dimension(:,:) :: U1_ZX,U2_ZX,U3_ZX,TH_ZX,K_E_ZX
       real(r4) :: g1vtk(1:np1),g3vtk(1:np3), var_2(1:3*np1*np3)

!       th_wh(1:np1,1:np2,N_TH+1)
!       real(r4),allocatable,dimension(:,:) :: var1
!       real(r8),allocatable,dimension(:,:,:) :: th_wh
       !STATUS VARIABLES
       integer :: s11


       allocate (U1_ZX(1:(NXP+1)*NP,1:NZ+2))
       allocate (U2_ZX(1:(NXP+1)*NP,1:NZ+2))
       allocate (U3_ZX(1:(NXP+1)*NP,1:NZ+2))
       allocate (TH_ZX(1:(NXP+1)*NP,1:NZ+2))
       allocate (K_E_ZX(1:(NXP+1)*NP,1:NZ+2))

       U1_ZX(:,:)=0.0d0
       U2_ZX(:,:)=0.0d0
       U3_ZX(:,:)=0.0d0
       TH_ZX(:,:)=0.0d0
      K_E_ZX(:,:)=0.0d0

!      jpn1 = int(NY)
       CALL  MPI_COMBINE (U1X(0:NXP,0:NZ+1,jpn1),U1_ZX,NXP+1,NZ+2)
       CALL  MPI_COMBINE (U2X(0:NXP,0:NZ+1,jpn1),U2_ZX,NXP+1,NZ+2)
       CALL  MPI_COMBINE (U3X(0:NXP,0:NZ+1,jpn1),U3_ZX,NXP+1,NZ+2)
       CALL  MPI_COMBINE (THX(0:NXP,0:NZ+1,jpn1,1),TH_ZX,NXP+1,NZ+2)
        

        do k=1,np1
         g1vtk(k) = dx(1)*(k-1)
        enddo
       
       do j=1,np3
         g3vtk(j) = GZF(J)
       enddo


       if (rank.eq.0)then
       outputDIR='plane_data_zx/'
       basename = 'data_parav_zx'
      

        write(OutFileName,'(a,a,i4.4,a2,i5.5,a4)') trim(outputDIR), &
        trim(basename)//"_p",jpn1,"_n",kk,".vtk"

!       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR), &
!       trim(basename)//"_n",kk,".vtk"  

       open(unit=13,file=OutFileName,access='stream',  &
       form='unformatted',status='unknown',             &
       convert='big_endian',iostat=s11)


        WRITE(6,*) 'Wrting in ', OutFileName,' for paraview for plane jpn1=', jpn1
!HEADER: note termination with char(10)

        write(13) "# vtk DataFile Version 3.0"//char(10)
        write(13) trim(BaseName)//char(10)
        write(13) "BINARY"//char(10)
        write(13) "DATASET RECTILINEAR_GRID"//char(10)

        write(ss,fmt='(A10,3I5)') "DIMENSIONS",np1,np3,1
        write(13) ss//char(10)
!Xgrid
        write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",np1," float"
        write(13) char(10)//ss//char(10)
        do k = 1,np1
        write(13) g1vtk(k)
        enddo
!Ygrid
        write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np3," float"
        write(13) char(10)//ss//char(10)
        do j = 1,np3
          write(13) g3vtk(j)
        enddo
    
!Zgrid
        write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",1," float"
        write(13) char(10)//ss//char(10)
        write(13) real(GYF(jpn1))
        

!Field
        j = 1
        do k = 1, np3
        do i = 1, np1
        var_2(j) =  real(U1_ZX(i,k))
        j = j+1
        var_2(j) =  real(U3_ZX(i,k))
        j = j+1
        var_2(j) =  real(U2_ZX(i,k))
        j = j+1
        K_E_ZX(i,k)=0.5d0*((real(U1_ZX(i,k)))**2+(real(U2_ZX(i,k)))**2+(real(U3_ZX(i,k)))**2)
        enddo
        enddo

!Field
        write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np3
        write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var_2
        write(13) "SCALARS Temp_div float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(TH_ZX(1:np1,1:np3))
        write(13) "SCALARS Kinetic_Energy float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(K_E_ZX(1:np1,1:np3))

!Close VTK File
        close(13)
        endif

        deallocate(U1_ZX,U2_ZX,U3_ZX,TH_ZX,K_E_ZX)

        return
        end subroutine plane_parav_plane_zx

        subroutine plane_parav_plane_zx_2(kk)

!      INCLUDE 'header_duct'
 use ntypes
 use Domain
 use Grid, only : dx,GZF, GYF
 use run_variable
 use variable_stat
 use mpi_var, only: rank
 implicit none

       integer i,j,k,np1,np3,kk,N
       !FILENAMES
       parameter(np1 = NX, np3 = NZ)
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss

       real(r8),allocatable,dimension(:,:) :: U1_ZX_2,U2_ZX_2,U3_ZX_2,TH_ZX_2,K_E_ZX_2
       real(r4) :: g1vtk(1:np1),g3vtk(1:np3), var_2(1:3*np1*np3)

!       th_wh(1:np1,1:np2,N_TH+1)
!       real(r4),allocatable,dimension(:,:) :: var1
!       real(r8),allocatable,dimension(:,:,:) :: th_wh
       !STATUS VARIABLES
       integer :: s11

       allocate (U1_ZX_2(1:(NXP+1)*NP,1:NZ+2))
       allocate (U2_ZX_2(1:(NXP+1)*NP,1:NZ+2))
       allocate (U3_ZX_2(1:(NXP+1)*NP,1:NZ+2))
       allocate (TH_ZX_2(1:(NXP+1)*NP,1:NZ+2))
       allocate (K_E_ZX_2(1:(NXP+1)*NP,1:NZ+2))

       U1_ZX_2(:,:)=0.0d0
       U2_ZX_2(:,:)=0.0d0
       U3_ZX_2(:,:)=0.0d0
       TH_ZX_2(:,:)=0.0d0
      K_E_ZX_2(:,:)=0.0d0

!      jpn2 = int(NY-16)
       CALL  MPI_COMBINE (U1X(0:NXP,0:NZ+1,jpn2),U1_ZX_2,NXP+1,NZ+2)
       CALL  MPI_COMBINE (U2X(0:NXP,0:NZ+1,jpn2),U2_ZX_2,NXP+1,NZ+2)
       CALL  MPI_COMBINE (U3X(0:NXP,0:NZ+1,jpn2),U3_ZX_2,NXP+1,NZ+2)
       CALL  MPI_COMBINE (THX(0:NXP,0:NZ+1,jpn2,1),TH_ZX_2,NXP+1,NZ+2)

       do k=1,np1
         g1vtk(k) = dx(1)*(k-1)
       enddo

       do j=1,np3
         g3vtk(j) = GZF(J)
       enddo

       if (rank.eq.0)then
       outputDIR='plane_data_zx/'
       basename = 'data_parav_zx'


       write(OutFileName,'(a,a,i4.4,a2,i5.5,a4)') trim(outputDIR), &
       trim(basename)//"_p",jpn2,"_n",kk,".vtk"

       open(unit=13,file=OutFileName,access='stream',  &
      form='unformatted',status='unknown',             &
      convert='big_endian',iostat=s11)

        WRITE(6,*) 'Wrting in ', OutFileName,' for paraview for plane jpn2=', jpn2
!HEADER: note termination with char(10)

        write(13) "# vtk DataFile Version 3.0"//char(10)
        write(13) trim(BaseName)//char(10)
        write(13) "BINARY"//char(10)
        write(13) "DATASET RECTILINEAR_GRID"//char(10)

        write(ss,fmt='(A10,3I5)') "DIMENSIONS",np1,np3,1
        write(13) ss//char(10)
!Xgrid
        write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",np1," float"
        write(13) char(10)//ss//char(10)
        do k = 1,np1
        write(13) g1vtk(k)
        enddo
!Ygrid
        write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np3," float"
        write(13) char(10)//ss//char(10)
        do j = 1,np3
          write(13) g3vtk(j)
        enddo

!Zgrid
        write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",1," float"
        write(13) char(10)//ss//char(10)
        write(13) real(GYF(jpn2))


!Field
        j = 1
        do k = 1, np3
        do i = 1, np1
        var_2(j) =  real(U1_ZX_2(i,k))
        j = j+1
        var_2(j) =  real(U3_ZX_2(i,k))
        j = j+1
        var_2(j) =  real(U2_ZX_2(i,k))
        j = j+1
        K_E_ZX_2(i,k)=0.5d0*((real(U1_ZX_2(i,k)))**2+(real(U2_ZX_2(i,k)))**2+(real(U3_ZX_2(i,k)))**2)
        enddo
        enddo

!Field
        write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np3
        write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var_2
        write(13) "SCALARS Temp_div float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(TH_ZX_2(1:np1,1:np3))
        write(13) "SCALARS Kinetic_Energy float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(K_E_ZX_2(1:np1,1:np3))

!Close VTK File
        close(13)
        endif

        deallocate(U1_ZX_2,U2_ZX_2,U3_ZX_2,TH_ZX_2,K_E_ZX_2)

        return
        end subroutine plane_parav_plane_zx_2

     subroutine plane_parav_plane_xy(kk)

!      INCLUDE 'header_duct'
 use ntypes
 use Domain
 use Grid, only : dx,GZF, GYF
 use run_variable
 use variable_stat
 use mpi_var, only: rank
 implicit none
             
       integer i,j,k,np1,np2,kk,N
       !FILENAMES
       parameter(np1 = NX, np2 = NY)
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss
       
       real(r8),allocatable,dimension(:,:) :: U1_XY,U2_XY,U3_XY,TH_XY,K_E
       real(r4) :: g1vtk(1:np1),g2vtk(1:np2), var_1(1:3*np1*np2)

!       th_wh(1:np1,1:np2,N_TH+1)
!       real(r4),allocatable,dimension(:,:) :: var1 
!       real(r8),allocatable,dimension(:,:,:) :: th_wh
       !STATUS VARIABLES
       integer :: s11


       allocate (U1_XY(1:(NXP+1)*NP,1:NY+2))
       allocate (U2_XY(1:(NXP+1)*NP,1:NY+2))
       allocate (U3_XY(1:(NXP+1)*NP,1:NY+2))
       allocate (TH_XY(1:(NXP+1)*NP,1:NY+2))
	allocate (K_E(1:(NXP+1)*NP,1:NY+2))

       U1_XY(:,:)=0.0d0
       U2_XY(:,:)=0.0d0
       U3_XY(:,:)=0.0d0 
       TH_XY(:,:)=0.0d0
	K_E(:,:)=0.0d0

!      kpn = int((NZ+1)/2)
       CALL  MPI_COMBINE (U1X(0:NXP,kpn1,0:NY+1),U1_XY,NXP+1,NY+2)
       CALL  MPI_COMBINE (U2X(0:NXP,kpn1,0:NY+1),U2_XY,NXP+1,NY+2)
       CALL  MPI_COMBINE (U3X(0:NXP,kpn1,0:NY+1),U3_XY,NXP+1,NY+2)
       CALL  MPI_COMBINE (THX(0:NXP,kpn1,0:NY+1,1),TH_XY,NXP+1,NY+2)

       do k=1,np1
         g1vtk(k) = dx(1)*(k-1)
       enddo

       do j=1,np2       
         g2vtk(j) = GYF(J)
       enddo


       if (rank.eq.0)then
       outputDIR='plane_data_xy/'
       basename = 'data_parav_xy'

         IF (IBM) THEN
         DO I=0,(NXP+1)*NP-1
         DO J=1,hill_ind_tot(I)+1
         U1_XY(I+1,J)=1.0d0
         U2_XY(I+1,J)=1.0d0
         U3_XY(I+1,J)=1.0d0
         TH_XY(I+1,J)=100.0d0        
         ENDDO
         ENDDO
         ENDIF

       write(OutFileName,'(a,a,i4.4,a2,i5.5,a4)') trim(outputDIR), &
       trim(basename)//"_p",kpn1,"_n",kk,".vtk"

       open(unit=13,file=OutFileName,access='stream',  &
      form='unformatted',status='unknown',             &
      convert='big_endian',iostat=s11)


        WRITE(6,*) 'Wrting in ', OutFileName,' for paraview for plane kpn1=', kpn1

!HEADER: note termination with char(10)
	
	write(13) "# vtk DataFile Version 3.0"//char(10)
	write(13) trim(BaseName)//char(10)
	write(13) "BINARY"//char(10)
	write(13) "DATASET RECTILINEAR_GRID"//char(10)

	write(ss,fmt='(A10,3I5)') "DIMENSIONS",np1,1,np2
	write(13) ss//char(10)
!Xgrid
	write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",np1," float"
	write(13) char(10)//ss//char(10)
        do k = 1,np1
        write(13) g1vtk(k)
        enddo
!Ygrid
	write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",1," float"
	write(13) char(10)//ss//char(10)
	write(13) real(GZF(kpn1))
!Zgrid
	write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np2," float"
	write(13) char(10)//ss//char(10)
	do j = 1,np2
	  write(13) g2vtk(j)
	enddo

!Field
        k = 1
        do j = 1, np2
        do i = 1, np1
        var_1(k) =  real(U1_XY(i,j))
        k = k+1
        var_1(k) =  real(U3_XY(i,j))
        k = k+1
        var_1(k) =  real(U2_XY(i,j))
        k = k+1
        K_E(i,j) =0.5d0*((real(U1_XY(i,j)))**2+(real(U2_XY(i,j)))**2+(real(U3_XY(i,j)))**2)
        enddo
        enddo
       
!Field
	write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
	write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var_1
        write(13) "SCALARS Temp_div float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(TH_XY(1:np1,1:np2))
        write(13) "SCALARS Kinetic_Energy float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(K_E(1:np1,1:np2))
        
!Close VTK File
	close(13)
        endif

        deallocate(U1_XY,U2_XY,U3_XY,TH_XY,K_E)

        return
        end subroutine plane_parav_plane_xy

      subroutine plane_parav_plane_xy_2(kk)

!      INCLUDE 'header_duct'
 use ntypes
 use Domain
 use Grid, only : dx,GZF, GYF
 use run_variable
 use variable_stat
 use mpi_var, only: rank
 implicit none
             
       integer i,j,k,np1,np2,kk,N
       !FILENAMES
       parameter(np1 = NX, np2 = NY)
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss
       
       real(r8),allocatable,dimension(:,:) :: U1_XY,U2_XY,U3_XY,TH_XY,K_E
       real(r4) :: g1vtk(1:np1),g2vtk(1:np2), var_1(1:3*np1*np2)

!       th_wh(1:np1,1:np2,N_TH+1)
!       real(r4),allocatable,dimension(:,:) :: var1 
!       real(r8),allocatable,dimension(:,:,:) :: th_wh
       !STATUS VARIABLES
       integer :: s11


       allocate (U1_XY(1:(NXP+1)*NP,1:NY+2))
       allocate (U2_XY(1:(NXP+1)*NP,1:NY+2))
       allocate (U3_XY(1:(NXP+1)*NP,1:NY+2))
       allocate (TH_XY(1:(NXP+1)*NP,1:NY+2))
        allocate (K_E(1:(NXP+1)*NP,1:NY+2))

       U1_XY(:,:)=0.0d0
       U2_XY(:,:)=0.0d0
       U3_XY(:,:)=0.0d0 
       TH_XY(:,:)=0.0d0
        K_E(:,:)=0.0d0

!      kpn = int((NZ+1)/2)
       CALL  MPI_COMBINE (U1X(0:NXP,kpn2,0:NY+1),U1_XY,NXP+1,NY+2)
       CALL  MPI_COMBINE (U2X(0:NXP,kpn2,0:NY+1),U2_XY,NXP+1,NY+2)
       CALL  MPI_COMBINE (U3X(0:NXP,kpn2,0:NY+1),U3_XY,NXP+1,NY+2)
       CALL  MPI_COMBINE (THX(0:NXP,kpn2,0:NY+1,1),TH_XY,NXP+1,NY+2)

       do k=1,np1
         g1vtk(k) = dx(1)*(k-1)
       enddo

       do j=1,np2       
         g2vtk(j) = GYF(J)
       enddo


       if (rank.eq.0)then
       outputDIR='plane_data_xy/'
       basename = 'data_parav_xy'

         IF (IBM) THEN
         DO I=0,(NXP+1)*NP-1
         DO J=1,hill_ind_tot(I)+1
         U1_XY(I+1,J)=1.0d0
         U2_XY(I+1,J)=1.0d0
         U3_XY(I+1,J)=1.0d0
         TH_XY(I+1,J)=100.0d0        
         ENDDO
         ENDDO
         ENDIF

       write(OutFileName,'(a,a,i4.4,a2,i5.5,a4)') trim(outputDIR), &
       trim(basename)//"_p",kpn2,"_n",kk,".vtk"

       open(unit=13,file=OutFileName,access='stream',  &
      form='unformatted',status='unknown',             &
      convert='big_endian',iostat=s11)


        WRITE(6,*) 'Wrting in ', OutFileName,' for paraview for plane kpn2=', kpn2

!HEADER: note termination with char(10)
        
        write(13) "# vtk DataFile Version 3.0"//char(10)
        write(13) trim(BaseName)//char(10)
        write(13) "BINARY"//char(10)
        write(13) "DATASET RECTILINEAR_GRID"//char(10)

        write(ss,fmt='(A10,3I5)') "DIMENSIONS",np1,1,np2
        write(13) ss//char(10)
!Xgrid
        write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",np1," float"
        write(13) char(10)//ss//char(10)
        do k = 1,np1
        write(13) g1vtk(k)
        enddo
!Ygrid
        write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",1," float"
        write(13) char(10)//ss//char(10)
        write(13) real(GZF(kpn2))
!Zgrid
        write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np2," float"
        write(13) char(10)//ss//char(10)
        do j = 1,np2
          write(13) g2vtk(j)
        enddo

!Field
        k = 1
        do j = 1, np2
        do i = 1, np1
        var_1(k) =  real(U1_XY(i,j))
        k = k+1
        var_1(k) =  real(U3_XY(i,j))
        k = k+1
        var_1(k) =  real(U2_XY(i,j))
        k = k+1
        K_E(i,j)=0.5d0*((real(U1_XY(i,j)))**2+(real(U2_XY(i,j)))**2+(real(U3_XY(i,j)))**2)
        enddo
        enddo
       
!Field
        write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
        write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var_1
        write(13) "SCALARS Temp_div float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(TH_XY(1:np1,1:np2))
        write(13) "SCALARS Kinetic_Energy float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(K_E(1:np1,1:np2))
        
!Close VTK File
        close(13)
        endif

        deallocate(U1_XY,U2_XY,U3_XY,TH_XY,K_E)

        return
        end subroutine plane_parav_plane_xy_2

      subroutine plane_parav_vel(kk)

!      INCLUDE 'header_duct'
 use ntypes
 use Domain
 use Grid, only : dx,GZF, GYF
 use run_variable
 use variable_stat
 implicit none
             
       integer i,j,k,np1,np2,kk,N
       !FILENAMES
       parameter(np1 = NZ, np2 = NY)
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss
       
 !      real(r4),allocatable,dimension(:) :: g1vtk,g2vtk
       real(r4)  g1vtk(1:np1),g2vtk(1:np2), var_1(1:3*np1*np2)
!      th_wh(1:np1,1:np2,N_TH+1)
!       real(r4),allocatable,dimension(:,:) :: var1 
!       real(r8),allocatable,dimension(:,:,:) :: th_wh
       !STATUS VARIABLES
       integer :: s11

!       allocate (th_wh(0:NZ+1,0:NY+1,1:N_TH+1))
!       allocate (var1(0:NZ+1,0:NY+1))
!       allocate(g1vtk(0:NZ+1), g2vtk(0:NY+1))



!       do j=0,NY+1
!         g1vtk(j) = xpoint(k,j)
!       enddo


       do k=1,np1
         g1vtk(k) = GZF(k)
       enddo

       do j=1,np2       
         g2vtk(j) = GYF(J)
       enddo

!         DO N=1,N_TH
!         th_wh(k,j,N)=dble(CTHX(0,k,j,N)) +  TH_BAR(K,j,N) !
!         ENDDO


 
!      IF (Non_linear_ST) THEN
!          call density_TC(.true.,.true.)
!       ELSE 
!       do j=0,NY+1
!       do k=0,NZ+1
!         th_wh(k,j,N_TH+1)= -alpha_w*dble(CTHX(0,K,J,1)) + gamma_w*dble(CTHX(0,K,J,2)) +  TH_BAR(K,j,N_TH+1) !
!       enddo
!       enddo
!       ENDIF

       outputDIR='plane_data_yz/'
       basename = 'data_parav'


       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR),  &
      trim(basename)//"_n",kk,".vtk"

       open(unit=13,file=OutFileName,access='stream', &
      form='unformatted',status='unknown',            &
      convert='big_endian',iostat=s11)


        WRITE(6,*) 'Wrting in ', OutFileName,' for paraview'
!HEADER: note termination with char(10)
	
	write(13) "# vtk DataFile Version 3.0"//char(10)
	write(13) trim(BaseName)//char(10)
	write(13) "BINARY"//char(10)
	write(13) "DATASET RECTILINEAR_GRID"//char(10)

	write(ss,fmt='(A10,3I5)') "DIMENSIONS",1,np1,np2
	write(13) ss//char(10)
!Xgrid
	write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",1," float"
	write(13) char(10)//ss//char(10)
	write(13) 0.01
!Ygrid
	write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np1," float"
	write(13) char(10)//ss//char(10)
	do k = 1,np1
	write(13) g1vtk(k)
	enddo
!Zgrid
	write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np2," float"
	write(13) char(10)//ss//char(10)
	do j = 1,np2
	  write(13) g2vtk(j)
	enddo

!Field
        k = 1
        do j = 1, np2
        do i = 1, np1
        var_1(k) =  real(dble(CU1X(0,i,j)))
        k = k+1
        var_1(k) =  real(dble(CU3X(0,i,j)))
        k = k+1
        var_1(k) =  real(dble(CU2X(0,i,j)))
        k = k+1
         DO N=1,N_TH
         th_wh(i,j,N)=dble(CTHX(0,i,j,N)) +  THBAR(j,N) !
         ENDDO
        enddo
        enddo
 
        do j = 1, np2
        do i = 1, np1
         th_wh(i,j,N_TH+1)=-alpha_T*dble(CTHX(0,i,j,1))+THBAR(j,N_TH+1) 
        enddo
        enddo
     
        
!Field
	write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
	write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var_1
        write(13) "SCALARS Temp_div float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(dble(CTHX(0,1:np1,1:np2,1)))
        write(13) "SCALARS Temp_wh float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(th_wh(1:np1,1:np2,1))
        write(13) "SCALARS Density_wh float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(th_wh(1:np1,1:np2,2))
!         write(13) "SCALARS Salinity float 1"//char(10)
!         write(13) "LOOKUP_TABLE default"//char(10)
!         write(13) real(dble(CTHX(0,0:NZ+1,0:NY+1,2)))
!         write(13) "SCALARS Temp_w float 1"//char(10)
!         write(13) "LOOKUP_TABLE default"//char(10)
!         write(13) real(th_wh(0:NZ+1,0:NY+1,1))
!         write(13) "SCALARS Salinity_w float 1"//char(10)
!         write(13) "LOOKUP_TABLE default"//char(10)
!         write(13) real(th_wh(0:NZ+1,0:NY+1,2))
!         write(13) "SCALARS Density_w float 1"//char(10)
!         write(13) "LOOKUP_TABLE default"//char(10)
!         write(13) real(th_wh(0:NZ+1,0:NY+1,3))
         write(13) "SCALARS Pressure float 1"//char(10)
         write(13) "LOOKUP_TABLE default"//char(10)
         write(13) real(dble(CPX(0,1:np1,1:np2)))
        write(13) "SCALARS Urms float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(wrms(1:np1,1:np2))
        write(13) "SCALARS Wrms float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(vrms(1:np1,1:np2))
        write(13) "SCALARS Z_star float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(Z_star_bar(1:np1,1:np2))
        write(13) "SCALARS back_pot float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(back_potential(1:np1,1:np2))

!Close VTK File
	close(13)

!deallocate(g1vtk,g2vtk,var1 )
!if (allocated(SPplane) ) deallocate(SPplane)
!if (allocated(DPplane) ) deallocate(DPplane)

        return
        end subroutine plane_parav_vel 






subroutine plot_para_tke(kk)

!      INCLUDE 'header_duct'
 use ntypes
 use Domain
 use Grid, only : dx,GZF, GYF
 use run_variable
 use variable_stat
 use les_chan_var
 implicit none
             
       integer i,j,k,np1,np2,kk,N
       !FILENAMES
       parameter(np1 = NZ, np2 = NY)
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss


      real(r4)  g1vtk(1:np1),g2vtk(1:np2), var_1(1:3*np1*np2)
      real(r8) :: var_C_dyn(1:np1,1:np2)
      integer :: s11




!       do j=0,NY+1
!         g1vtk(j) = xpoint(k,j)
!       enddo


       do k=1,np1
         g1vtk(k) = GZF(k)
       enddo

       do j=1,np2       
         g2vtk(j) = GYF(J)
       enddo

!        do j=1,np2
!        do k=1,np1
!         if(denominator_sum(k,j) .ne. 0.) then
!          var_C_dyn(k,j) = -0.50d0*numerator_sum(k,j)/denominator_sum(k,j)
!         else
!          var_C_dyn(k,j) = 0.0
!         endif 
!          if(var_C_dyn(k,j) .lt. 0.d0 ) then
!            var_C_dyn(k,j) = 0.0d0
!          endif
!          var_C_dyn(k,j) = DELTA_Y(k,j)**2.d0*Sbar_MEAN(k,j)
!        enddo
!        enddo




!         DO N=1,N_TH
!         th_wh(k,j,N)=dble(CTHX(0,k,j,N)) +  TH_BAR(K,j,N) !
!         ENDDO


 
!      IF (Non_linear_ST) THEN
!          call density_TC(.true.,.true.)
!       ELSE 
!       do j=0,NY+1
!       do k=0,NZ+1
!         th_wh(k,j,N_TH+1)= -alpha_w*dble(CTHX(0,K,J,1)) + gamma_w*dble(CTHX(0,K,J,2)) +  TH_BAR(K,j,N_TH+1) !
!       enddo
!       enddo
!       ENDIF

       outputDIR='plane_tke/'
       basename = 'data_tke_parav'


       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR),  &
      trim(basename)//"_n",kk,".vtk"

       open(unit=13,file=OutFileName,access='stream', &
      form='unformatted',status='unknown',            &
      convert='big_endian',iostat=s11)


        WRITE(6,*) 'Wrting in ', OutFileName,' for paraview'
!HEADER: note termination with char(10)
	
	write(13) "# vtk DataFile Version 3.0"//char(10)
	write(13) trim(BaseName)//char(10)
	write(13) "BINARY"//char(10)
	write(13) "DATASET RECTILINEAR_GRID"//char(10)

!        write(ss,fmt='(A4,2I4,A6)') "TIME",1,1, " double"
!        write(13) char(10)//ss//char(10)
!        write(13) real(TIME/(60.00*60.00))

	write(ss,fmt='(A10,3I5)') "DIMENSIONS",1,np1,np2
	write(13) ss//char(10)
        
!Xgrid
	write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",1," float"
	write(13) char(10)//ss//char(10)
	write(13) 1.0
!Ygrid
	write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np1," float"
	write(13) char(10)//ss//char(10)
	do k = 1,np1
	write(13) g1vtk(k)
	enddo
!Zgrid
	write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np2," float"
	write(13) char(10)//ss//char(10)
	do j = 1,np2
	  write(13) g2vtk(j)
	enddo

!Field
        k = 1
        do j = 1, np2
        do i = 1, np1
        var_1(k) =  real(dble(CR1X(0,i,j)))
        k = k+1
        var_1(k) =  real(dble(CR3X(0,i,j)))
        k = k+1
        var_1(k) =  real(dble(CR2X(0,i,j)))
        k = k+1
!         DO N=1,N_TH
!         th_wh(i,j,N)=dble(CTHX(0,i,j,N)) +  THBAR(j,N) !
!         ENDDO
        enddo
        enddo
 
!        do j = 1, np2
!        do i = 1, np1
!        th_wh(i,j,N_TH+1)=-alpha_T*dble(CTHX(0,i,j,1))+THBAR(j,N_TH+1) 
!        enddo
!        enddo
     
        
!Field
	write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
	write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var_1
        write(13) "SCALARS tke float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_mean(1:np1,1:np2))
        write(13) "SCALARS dtkedt float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_1(1:np1,1:np2))
        write(13) "SCALARS Tur_advection float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_2(1:np1,1:np2))
        write(13) "SCALARS tur_production float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_3(1:np1,1:np2))
        write(13) "SCALARS viscous_dissipation float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_4(1:np1,1:np2))
        write(13) "SCALARS Buoyancy_flux float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_5(1:np1,1:np2,1))
        write(13) "SCALARS pressure_transport float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_6_1(1:np1,1:np2))
        write(13) "SCALARS trubulent_transport float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_6_2(1:np1,1:np2))
        write(13) "SCALARS dissipation float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_7(1:np1,1:np2))
        IF (LES) THEN
        write(13) "SCALARS Kappa_T_mean float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(Kappa_T_MEAN(1:np1,1:np2,1))
        write(13) "SCALARS NU_T_mean float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(NU_T_mean(1:np1,1:np2))
        write(13) "SCALARS Prod_sgs float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_sgs_p(1:np1,1:np2))
        write(13) "SCALARS Dissip_sgs float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_sgs_diss(1:np1,1:np2)) 
        ENDIF    

!Close VTK File
	close(13)




	return
        end subroutine plot_para_tke


       subroutine plot_para_eng(kk)

!      INCLUDE 'header_duct'
 use ntypes
 use Domain
 use Grid, only : dx,GZF, GYF
 use run_variable
 use variable_stat
 implicit none
             
       integer i,j,k,np1,np2,kk,N
       !FILENAMES
       parameter(np1 = NZ, np2 = NY)
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss


      real(r4)  g1vtk(1:np1),g2vtk(1:np2), var_1(1:3*np1*np2)
      integer :: s11

      return
end subroutine plot_para_eng


    subroutine plane_parav_3D(kk)

!      INCLUDE 'header_duct'
 use ntypes
 use Domain
 use Grid, only : dx,LX,LZ,LY,GZF, GYF, GX
 use run_variable
 use variable_stat
 use mpi_var, only: rank, IERROR
 use TIME_STEP_VAR

implicit none
   integer :: I, J, K, tot_n, send_rank, en, st, kk
   real(r8),allocatable,dimension(:,:,:,:) :: array_temp_3D
   real(r8),allocatable,dimension(:,:,:)      :: array_3D
 


      !FILENAMES
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss
        
       real(r8) :: TH, UX_1, UX_2, UX_3, PR_3D

       !STATUS VARIABLES
       integer :: s11
        
        allocate (array_temp_3D(0:NXP,0:NZ+1,0:NY+1,1:NP))
        allocate (array_3D(0:((NXP+1)*NP-1), 0:NZ+1, 0:NY+1))

        array_temp_3D(:,:,:,:) = 0.0d0
        array_3D(:,:,:) = 0.0d0

        send_rank = 0
        tot_n = (NXP+1)*(NY+2)*(NZ+2)
       CALL MPI_GATHER(THX(0:NXP,0:NZ+1,0:NY+1,1),tot_n,MPI_DOUBLE_PRECISION, & 
       array_temp_3D,tot_n,MPI_DOUBLE_PRECISION,send_rank,MPI_COMM_WORLD,IERROR)

 IF ( RANK .EQ. 0) THEN
        array_3D(:,:,:) = 0.0d0
        st=1
        DO I = 1,NP
          en=st+NXP
       IF (en<=NX) THEN
         DO J =0,NY+1
         DO K=0,NZ+1
           array_3D(st:en,K,J) =  array_temp_3D(:,K,J,I)
         ENDDO
         ENDDO
        st=en+1
         ELSE
         END IF
        ENDDO
      
         outputDIR='plane_data_3D/'
       basename = 'data_binary_3D'

       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR), &
       trim(basename)//"_n",kk,".pln"

       open(unit=13,file=OutFileName,access='stream',  &
       form='unformatted',status='unknown',             &
       convert='big_endian',iostat=s11)

        WRITE(6,*) 'Wrting TIME, Grid,  TH in binary ', OutFileName
        write(13) time,NX,NZ+2,NY+2,DX(1),LX,LZ,LY    
        write(13) (((array_3D(I,K,J),I=1,NX),K=0,NZ+1),J=0,NY+1)
        open(66,file='plane_data_3D/time_bulk_3D.txt',form='formatted', &
       status='unknown')

      DO I = 1,NX
        TH = array_3D(I, 125,113)

        write(66,565) GX(I) , TH
        END DO
565    format(f12.5,f13.8)        

         END IF
        
        array_temp_3D(:,:,:,:) = 0.0d0
        array_3D(:,:,:) = 0.0d0

     CALL MPI_GATHER(U1X(0:NXP,0:NZ+1,0:NY+1),tot_n,MPI_DOUBLE_PRECISION, & 
     array_temp_3D,tot_n,MPI_DOUBLE_PRECISION,send_rank,MPI_COMM_WORLD,IERROR)

IF ( RANK .EQ. 0) THEN
        array_3D(:,:,:) = 0.0d0
        st=1
        DO I = 1,NP
          en=st+NXP
       IF (en<=NX) THEN
         DO J =0,NY+1
         DO K=0,NZ+1
           array_3D(st:en,K,J) =  array_temp_3D(:,K,J,I)
         ENDDO
         ENDDO
        st=en+1
         ELSE
         END IF
        ENDDO
    
         WRITE(6,*) 'Wrting U1X in binary ', OutFileName

        write(13) (((array_3D(I,K,J),I=1,NX),K=0,NZ+1),J=0,NY+1)

      DO I = 1,NX
        UX_1 = array_3D(I, 125,113)

        write(66,566) GX(I), UX_1
        END DO
566    format(f12.8,f13.8)

         END IF

        array_temp_3D(:,:,:,:) = 0.0d0
        array_3D(:,:,:) = 0.0d0
         
       CALL MPI_GATHER(U2X(0:NXP,0:NZ+1,0:NY+1),tot_n,MPI_DOUBLE_PRECISION, &
       array_temp_3D,tot_n,MPI_DOUBLE_PRECISION,send_rank,MPI_COMM_WORLD,IERROR)

IF ( RANK .EQ. 0) THEN
        array_3D(:,:,:) = 0.0d0
        st=1
        DO I = 1,NP
          en=st+NXP
       IF (en<=NX) THEN
         DO J =0,NY+1
         DO K=0,NZ+1
           array_3D(st:en,K,J) =  array_temp_3D(:,K,J,I)
         ENDDO
         ENDDO
        st=en+1
         ELSE
         END IF
        ENDDO

        WRITE(6,*) 'Wrting U2X in binary ', OutFileName

        write(13) (((array_3D(I,K,J),I=1,NX),K=0,NZ+1),J=0,NY+1)

         DO I = 1,NX
        UX_2 = array_3D(I, 125,113)

        write(66,567) GX(I),UX_2
        END DO
567    format(f12.8,f13.8)

         END IF

        array_temp_3D(:,:,:,:) = 0.0d0
        array_3D(:,:,:) = 0.0d0

        CALL MPI_GATHER(U3X(0:NXP,0:NZ+1,0:NY+1),tot_n,MPI_DOUBLE_PRECISION, &
      array_temp_3D,tot_n,MPI_DOUBLE_PRECISION,send_rank,MPI_COMM_WORLD,IERROR)

IF ( RANK .EQ. 0) THEN
        array_3D(:,:,:) = 0.0d0
        st=1
        DO I = 1,NP
          en=st+NXP
       IF (en<=NX) THEN
         DO J =0,NY+1
         DO K=0,NZ+1
           array_3D(st:en,K,J) =  array_temp_3D(:,K,J,I)
         ENDDO
         ENDDO
        st=en+1
         ELSE
         END IF
        ENDDO

        WRITE(6,*) 'Wrting U3X in binary ', OutFileName

        write(13) (((array_3D(I,K,J),I=1,NX),K=0,NZ+1),J=0,NY+1)

         DO I = 1,NX
        UX_3 = array_3D(I, 125,113)

        write(66,568) GX(I), UX_3
        END DO
568    format(f12.8,f13.8)
        END IF
         
        array_temp_3D(:,:,:,:) = 0.0d0
         array_3D(:,:,:) = 0.0d0

      CALL MPI_GATHER(PX(0:NXP,0:NZ+1,0:NY+1),tot_n,MPI_DOUBLE_PRECISION, &
      array_temp_3D,tot_n,MPI_DOUBLE_PRECISION,send_rank,MPI_COMM_WORLD,IERROR)

IF ( RANK .EQ. 0) THEN
        array_3D(:,:,:) = 0.0d0
        st=1
        DO I = 1,NP
          en=st+NXP
       IF (en<=NX) THEN
         DO J =0,NY+1
         DO K=0,NZ+1
           array_3D(st:en,K,J) =  array_temp_3D(:,K,J,I)
         ENDDO
         ENDDO
        st=en+1
         ELSE
         END IF
        ENDDO

        WRITE(6,*) 'Wrting PRX in binary ', OutFileName

        write(13) (((array_3D(I,K,J),I=1,NX),K=0,NZ+1),J=0,NY+1)

       DO I = 1,NX
        PR_3D = array_3D(I, 125,113)

        write(66,569) GX(I), PR_3D
        END DO
569    format(f12.8,f13.8)

       close(13)
       close(66)
       END IF

        deallocate(array_3D, array_temp_3D)
        
        return
        end subroutine plane_parav_3D

      subroutine plane_3D_binary(kk)

use ntypes
use Domain
use Grid, only : dx
use run_variable, only : U1X,U2X,U3X,THX,TIME,TIME_STEP
use mg_vari, only : INIT_FLAG
use TIME_STEP_VAR, only: DELTA_T
use mpi_var, only: rank
implicit none

      integer  ind
      integer  i,j,k, NI,kmax,kk
!      character(len=100) :: InFileName, OutFileName
!      character(len=100) :: dataDIR, outputDIR
!      character(len=80) :: basename
      CHARACTER*35 file_name
      CHARACTER*6 file_num
      CHARACTER*3 PID

       I=1+RANK       

       PID = CHAR(MOD(I,1000)/100+48)     &
             //CHAR(MOD(I,100)/10+48)    &
             //CHAR(MOD(I,10)+48) 

       file_num = CHAR(MOD(kk,100000)/10000+48) &
             //CHAR(MOD(kk,10000)/1000+48) &
             //CHAR(MOD(kk,1000)/100+48)     &
             //CHAR(MOD(kk,100)/10+48)       &
             //CHAR(MOD(kk,10)+48)//'_'

       file_name = 'plane_3D/data_'//file_num//PID//'.pln'

!       outputDIR='plane_3D/'
!       basename = 'data_3D'


!       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR),  &
!      trim(basename)//trim(PID)//"_n",kk,".pln"

        NI = min(NXP,NXP_L)

!        write(6,*) 'Saving span-wise plane information'
       open(22,file=file_name,form='unformatted',status='unknown')
       write(22) TIME,DELTA_T,TIME_STEP,NI+1,NZ+2,NY+2
       write(22) dx(1),U3X(0:NI,0:NZ+1,0:NY+1),               &
            U2X(0:NI,0:NZ+1,0:NY+1), U1X(0:NI,0:NZ+1,0:NY+1), &
            THX(0:NI,0:NZ+1,0:NY+1,1)
       close(22)

!      deallocate (xpoint, ypoint)

      return
      end subroutine plane_3D_binary


subroutine plane_ZX_binary(kk,ind)

use ntypes
use Domain
use Grid, only : dx, GZF
use run_variable, only : U1X,U2X,U3X,THX,PX,TIME,TIME_STEP
use mpi_var, only : RANK
use TIME_STEP_VAR, only: DELTA_T
implicit none

      integer  ind
      integer  i,j, NI,kk
!      CHARACTER*29 file_name
!      real(r8),allocatable,dimension(:,:)   :: xpoint, ypoint
!      allocate (xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1))

      CHARACTER*35 file_name
      CHARACTER*6 file_num
      CHARACTER*3 PID

       I=1+RANK

       PID = CHAR(MOD(I,1000)/100+48)     &
             //CHAR(MOD(I,100)/10+48)    &
             //CHAR(MOD(I,10)+48)

       file_num = CHAR(MOD(kk,100000)/10000+48) &
             //CHAR(MOD(kk,10000)/1000+48) &
             //CHAR(MOD(kk,1000)/100+48)     &
             //CHAR(MOD(kk,100)/10+48)       &
             //CHAR(MOD(kk,10)+48)//'_'

       file_name = 'plane_data_ZX/data_ZX_'//file_num//PID//'.pln'

!      

       J=ind
       NI = min(NXP,NXP_L)

!        write(6,*) 'Saving span-wise plane information'
       open(22,file=file_name,form='unformatted',status='unknown')
       write(22) TIME,DELTA_T,TIME_STEP, NI+1,NZ+2
       write(22) GZF(0:NZ+1),dx(1),U3X(0:NI,0:NZ+1,J), &
                        U2X(0:NI,0:NZ+1,J), U1X(0:NI,0:NZ+1,J), &
                        THX(0:NI,0:NZ+1,J,1)
       close(22)


      return
      end subroutine plane_ZX_binary

subroutine plane_ZX_binary_2(kk,ind)

use ntypes
use Domain
use Grid, only : dx, GZF
use run_variable, only : U1X,U2X,U3X,THX,PX,TIME,TIME_STEP
use mpi_var, only : RANK
use TIME_STEP_VAR, only: DELTA_T
implicit none

      integer  ind
      integer  i,j, NI,kk
!      CHARACTER*29 file_name
!      real(r8),allocatable,dimension(:,:)   :: xpoint, ypoint
!      allocate (xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1))

      CHARACTER*37 file_name
      CHARACTER*6 file_num
      CHARACTER*3 PID

       I=1+RANK

       PID = CHAR(MOD(I,1000)/100+48)     &
             //CHAR(MOD(I,100)/10+48)    &
             //CHAR(MOD(I,10)+48)

       file_num = CHAR(MOD(kk,100000)/10000+48) &
             //CHAR(MOD(kk,10000)/1000+48) &
             //CHAR(MOD(kk,1000)/100+48)     &
             //CHAR(MOD(kk,100)/10+48)       &
             //CHAR(MOD(kk,10)+48)//'_'

       file_name = 'plane_data_ZX_2/data_ZX_'//file_num//PID//'.pln'

!      

       J=ind
       NI = min(NXP,NXP_L)

!        write(6,*) 'Saving span-wise plane information'
       open(22,file=file_name,form='unformatted',status='unknown')
       write(22) TIME,DELTA_T,TIME_STEP, NI+1,NZ+2
       write(22) GZF(0:NZ+1),dx(1),U3X(0:NI,0:NZ+1,J), &
                        U2X(0:NI,0:NZ+1,J), U1X(0:NI,0:NZ+1,J), &
                        THX(0:NI,0:NZ+1,J,1)
       close(22)


      return
      end subroutine plane_ZX_binary_2




      subroutine plane_XY_binary(kk,ind)

use ntypes
use Domain
use Grid, only : dx, GYF
use run_variable, only : U1X,U2X,U3X,THX,PX,TIME,TIME_STEP
use mpi_var, only : RANK
use TIME_STEP_VAR, only: DELTA_T
implicit none

      integer  ind
      integer  i,j,k, NI,kk
!      CHARACTER*29 file_name
!      real(r8),allocatable,dimension(:,:)   :: xpoint, ypoint
!      allocate (xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1))

      CHARACTER*35 file_name
      CHARACTER*6 file_num
      CHARACTER*3 PID

       I=1+RANK

       PID = CHAR(MOD(I,1000)/100+48)     &
             //CHAR(MOD(I,100)/10+48)    &
             //CHAR(MOD(I,10)+48)

       file_num = CHAR(MOD(kk,100000)/10000+48) &
             //CHAR(MOD(kk,10000)/1000+48) &
             //CHAR(MOD(kk,1000)/100+48)     &
             //CHAR(MOD(kk,100)/10+48)       &
             //CHAR(MOD(kk,10)+48)//'_'

       file_name = 'plane_data_XY/data_XY_'//file_num//PID//'.pln'

!      

       k=ind
       NI = min(NXP,NXP_L)

!        write(6,*) 'Saving span-wise plane information'
       open(22,file=file_name,form='unformatted',status='unknown')
       write(22) TIME,DELTA_T,TIME_STEP, NI+1,NY+2
       write(22) GYF(0:NY+1),dx(1),U3X(0:NI,k,0:NY+1), &
                        U2X(0:NI,k,0:NY+1), U1X(0:NI,k,0:NY+1), &
                        THX(0:NI,k,0:NY+1,1)
                                                                                                                         
       close(22)


      return
      end subroutine plane_XY_binary

