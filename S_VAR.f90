MODULE SOLVER_VARIABLES
IMPLICIT NONE
!----------------------'PAR' VARIABLES--------------------------
INTEGER,PARAMETER:: NX=2048, NY=500, NZ=1
INTEGER,PARAMETER:: NMAX=2048, MST=22
REAL(8),PARAMETER:: ZL=1.D0
INTEGER,PARAMETER:: NFXX=7, NFXZ=0, NVXE0=3333, N1=33, N2=20
REAL(8)          :: XL,YL
REAL(8)          :: PI,KC,FR,AMP,UC,UC_OLD,UM
!---------------------NUMERICAL PARAMETER-----------------------
REAL(8)                                :: TID, DT
INTEGER                                :: NNN,INS,NST,IIO,IIC
REAL(8)                                :: AOA,THETA,STD
!------------------------FLOW PARAMETER-------------------------
REAL(8)                                :: REB,PRL,DTDY
REAL(8)                                :: UXBI,UXBO
!-----------------------------GRID------------------------------
REAL(8),PARAMETER                      :: TV=2.10D0, GR=1.018D0
REAL(8)                                :: GS
REAL(8),DIMENSION(0:NX+2)              :: X,X12
REAL(8),DIMENSION(0:NY+2)              :: Y,Y12,DY,DY12
REAL(8),DIMENSION(0:NZ+2)              :: Z,Z12
REAL(8)                                :: DX,DZ
!----------------------------FFT--------------------------------
REAL(8),DIMENSION(0:NX)                :: KX
REAL(8),DIMENSION(0:NZ)                :: KZ
REAL(8),DIMENSION(0:NX,0:NZ)           :: KXZ
!---------------------------TDMA--------------------------------
REAL(8),DIMENSION(1:NX,1:NY,1:NZ,3)    :: AG,ACT
REAL(8),DIMENSION(1:NY,3)              :: AC,AE!,ACHI
!----------------------VELOCITY AND B.C.------------------------
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: UX,UUX,UY,UUY,UZ
REAL(8),DIMENSION(0:NY+2,0:NZ+2)       :: VXI,VXO,VYI,VYO,VZI,VZO
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: TT
REAL(8),DIMENSION(0:NY+2,0:NZ+2)       :: TTI, TTO
REAL(8)                                :: VXM
!
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: UX_OLD,UY_OLD,TT_OLD
!-------------------------PRESSURE------------------------------
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: P
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: GXP,GYP,GZP
REAL(8),DIMENSION(0:NY+2,0:NZ+2)       :: PIN,POUT
!--------------------------DELTA U------------------------------
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: DUX,DUY,DUZ,DTT
REAL(8),DIMENSION(0:NX+2,0:NY+2):: DUDY
!-----------------------ADVECTION TERM--------------------------
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: FX,FX1,FY,FY1,FZ,FZ1
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: FT,FT1
!-------------------------STATISTICS----------------------------
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: DUXC,DUYC
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: DUZC,RX
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: RY,RZ
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: XMM,YMM
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: ZMM,UXC
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: UYC,UZC
REAL(8),DIMENSION(0:NY+2,0:MST)        :: ST
!---------------------------------------------------------------
REAL(8),DIMENSION(0:NX+2,0:NY+2,0:NZ+2):: OMEGA_Z
END MODULE
