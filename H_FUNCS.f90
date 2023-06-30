FUNCTION ARC(AS,X)
IMPLICIT NONE
REAL(8):: ARC,X
REAL(8),DIMENSION(0:3):: AS
ARC= AS(3)*X*X*X+AS(2)*X*X+AS(1)*X+AS(0)
END
!------------------------------------------------------------------------------
FUNCTION RTBIS(ARC,AS,XS1,XE2,XACC)
IMPLICIT NONE
INTEGER,PARAMETER     :: JMAX=40
INTEGER               :: J
REAL(8),DIMENSION(0:3):: AS
REAL(8)               :: FMID,F,RTBIS,XACC
REAL(8)               :: XS1,XE2,DX,XMID,ARC
!------------------------------------------------
FMID   = ARC(AS,XE2)
F      = ARC(AS,XS1)
IF(F*FMID.GE.0.) THEN
PRINT*, 'Root must be bracketed for bisection.'
END IF
IF(F.LT.0.)THEN
   RTBIS= XS1
   DX   = XE2-XS1
ELSE
RTBIS= XE2
DX   = XS1-XE2
ENDIF
DO J=1,JMAX
   DX   = DX*0.5D0
   XMID = RTBIS+DX
   FMID = ARC(AS,XMID)
   IF(FMID.LE.0.)RTBIS=XMID
   IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
END DO
PAUSE 'too many bisections'
END
!------------------------------------------------------------------------------
SUBROUTINE zbrak(fx,AS,x11,x22,n,xb1,xb2,nb)
IMPLICIT NONE
INTEGER,PARAMETER         :: NBMAX=20
INTEGER                   :: n,nb
REAL(8),EXTERNAL          :: fx
INTEGER                   :: i,nbb1
REAL(8)                   :: x11,x22
REAL(8),DIMENSION(1:NBMAX):: xb1,xb2
REAL(8),DIMENSION(0:3)    :: AS
REAL(8)                   :: dx,fc,xz
REAL(8)                   :: fp
nbb1=0
xz=x11
dx=(x22-x11)/n
fp=fx(AS,xz)
do i=1,n
   xz=xz+dx
   fc=fx(AS,xz)
   if(fc*fp.le.0.) then
      nbb1    = nbb1+1
      xb1(nbb1)= xz-dx
      xb2(nbb1)= xz
      if(nbb1.eq.nb)goto 1
    endif
    fp=fc
END DO
1     continue
nb=nbb1
return
END


