program main
implicit none

integer NP, NX, NXV, NX2V,NX2P, NXP, NZV,NZP,NKXV,NKX, NZ, NKXP
real(8) :: NXV_R

PARAMETER (NP = 49, NX=96, NKX = NX/3, NZ=513)

PARAMETER (NZV=anint(real(NZ+2)/real(NP))*NP)
PARAMETER (NXV=anint(real(NX+2)/real(NP))*NP)
PARAMETER (NXP=NXV/NP-1)
PARAMETER (NZP=NZV/NP-1)
PARAMETER (NKXV=anint(real(NKX+1)/real(NP))*NP)
PARAMETER (NKXP=NKXV/NP-1)
PARAMETER (NX2V=anint(real(NXV/2)/real(NP))*NP)
PARAMETER (NX2P=NX2V/NP-1)



write(6,*) 'Hello',NXV_R
WRITE(6,*) 'NX= ', NX, ' NXV= ', NXV, 'NX2P = ', NX2P, 'NKXP = ', NKXP, 'NKXV= ', NKXV, 'NX2V', NX2V, 'NXP',NXP
WRITE(6,*) 'NZ= ', NZ, ' NZV= ', NZV, 'NZP = ', NZP

stop
end
