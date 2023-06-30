!In this subroutine, the interception between the normal line from the points belonging -----------
!to the 3*3 packet and the physical boundary ------------------------------------------------------
SUBROUTINE H_INTERSECTION (H_NX,H_NY,X,Y,H_S_INTER)
!
USE IB_VARIABLES, ONLY: H_M_START,H_N_START,H_M_END,H_N_END,H_A_X,H_A_Y,H_DELTA_X
IMPLICIT NONE
INTEGER                             :: I,J
INTEGER                             :: H_NX,H_NY,NBB,K1
REAL(8)                             :: XACC,ROOT,RTBIS
REAL(8)                             :: X_INT1,X_INT2
REAL(8),EXTERNAL                    :: ARC
INTEGER,PARAMETER                   :: N_N=10,NBMAX=20
REAL(8),DIMENSION(0:3)              :: A_S
REAL(8),DIMENSION(1:NBMAX)          :: X_ST,X_EN
REAL(8),DIMENSION(0:H_NX+2)         :: X
REAL(8),DIMENSION(0:H_NY+2)         :: Y
REAL(8),DIMENSION(0:H_NX+2,0:H_NY+2):: H_S_INTER
!
ROOT=1000.0D0
DO J=H_N_START,H_N_END
   DO I=H_M_START,H_M_END
      A_S(3)=(2.0D0*H_A_Y(0)*H_A_Y(0)+2.0D0*H_A_X(0)*H_A_X(0))
      A_S(2)=(3.0D0*H_A_Y(0)*H_A_Y(1)+3.0D0*H_A_X(0)*H_A_X(1))
      A_S(1)=(H_A_Y(1)*H_A_Y(1)+2.0D0*H_A_Y(2)*H_A_Y(0)+H_A_X(1)*                                  &
              H_A_X(1)+2.0D0*H_A_X(2)*H_A_X(0)-2.0D0*                                              &
              H_A_Y(0)*(Y(J))-2.0D0*(X(I))*H_A_X(0))
      A_S(0)=(H_A_Y(2)*H_A_Y(1)+H_A_X(2)*H_A_X(1)-H_A_Y(1)*(Y(J))-H_A_X(1)*(X(I)))
!
      X_INT1= 0.0D0+H_DELTA_X/DBLE(1000000)
      X_INT2= 2.0D0*H_DELTA_X
      CALL ZBRAK(ARC,A_S,X_INT1,X_INT2,N_N,X_ST,X_EN,NBB)
      DO K1=1,NBB
         XACC= (1.0E-10)*(X_ST(K1)+X_EN(K1))/2.0D0
         ROOT= RTBIS(ARC,A_S,X_ST(K1),X_EN(K1),XACC)
      END DO
      IF (NBB==0) THEN
          ROOT=1000.0D0
      END IF
      H_S_INTER(I,J)= ROOT
   END DO
END DO
END SUBROUTINE
