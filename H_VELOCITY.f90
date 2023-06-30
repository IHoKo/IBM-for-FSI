!---------------------- IN THIS SUBROUTINE, THE VELOCITY OF THE GHOST POINTS ARE INTERPOLATED USING
!THE VALUE OF DATA POINTS THAT ARE INDENTIFIED IN 'INTERPOLATION_DATA' ----------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE H_VELOCITY(H_NX,H_NY,H_N_NB,I_IP,J_IP,X,Y,H_I_NB,H_J_NB,                                &
                      H_X_INT,H_Y_INT,H_DIS_INPOL,H_BODY,U)
USE IB_VARIABLES, ONLY: N_ORDER
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER                             :: I,J
INTEGER                             :: H_NX,H_NY
INTEGER                             :: H_N_NB
INTEGER,PARAMETER                   :: N_S= N_ORDER+1
REAL(8)                             :: DEL_0,W_N,ER
REAL(8),DIMENSION(1:2,1:H_NX)       :: U_INTER
REAL(8),DIMENSION(1:N_S)            :: DEL,VEL
REAL(8),DIMENSION(1:2,1:4)          :: A_NUM
REAL(8),DIMENSION(1:2)              :: A_DEN
INTEGER,DIMENSION(1:H_NX)           :: H_I_NB,H_J_NB
REAL(8),DIMENSION(1:H_NX)           :: H_BODY
REAL(8),DIMENSION(0:H_NX+2)         :: X
REAL(8),DIMENSION(0:H_NY+2)         :: Y
REAL(8),DIMENSION(1:3,1:H_NX)       :: H_DIS_INPOL
REAL(8),DIMENSION(1:2,1:H_NX)       :: H_X_INT,H_Y_INT
INTEGER,DIMENSION(1:4,1:2,1:H_NX)   :: I_IP,J_IP
REAL(8),DIMENSION(0:H_NX+2,0:H_NY+2):: U
!
DO I=1,H_N_NB
   DO J=1,2
      A_NUM(J,3)= DABS(X(I_IP(4,J,I))-H_X_INT(J,I))*                                               &
                  DABS(Y(J_IP(1,J,I))-H_Y_INT(J,I))
!
      A_NUM(J,4)= DABS(H_X_INT(J,I)-X(I_IP(1,J,I)))*                                               &
                  DABS(Y(J_IP(1,J,I))-H_Y_INT(J,I))
!
      A_NUM(J,1)= DABS(X(I_IP(4,J,I))-H_X_INT(J,I))*                                               &
                  DABS(H_Y_INT(J,I)-Y(J_IP(4,J,I)))
!
      A_NUM(J,2)= DABS(H_X_INT(J,I)-X(I_IP(1,J,I)))*                                               &
                  DABS(H_Y_INT(J,I)-Y(J_IP(4,J,I)))
!
      A_DEN(J)  = DABS(X(I_IP(4,J,I))-X(I_IP(1,J,I)))*                                             &
                  DABS(Y(J_IP(1,J,I))-Y(J_IP(4,J,I)))
!
      U_INTER(J,I)= A_NUM(J,3)/A_DEN(J)*U(I_IP(3,J,I),J_IP(3,J,I))+                                &
                    A_NUM(J,4)/A_DEN(J)*U(I_IP(4,J,I),J_IP(4,J,I))+                                &
                    A_NUM(J,1)/A_DEN(J)*U(I_IP(1,J,I),J_IP(1,J,I))+                                &
                    A_NUM(J,2)/A_DEN(J)*U(I_IP(2,J,I),J_IP(2,J,I))
   END DO
!------------------------------------------ INTERPOLATION -----------------------------------------
   DEL(1)= 0.0D0
   DEL(2)= H_DIS_INPOL(2,I)+H_DIS_INPOL(1,I)
   DEL(3)= H_DIS_INPOL(3,I)+H_DIS_INPOL(1,I)
   VEL(1)= H_BODY(I)
   VEL(2)= U_INTER(1,I)
   VEL(3)= U_INTER(2,I)
   DEL_0 = H_DIS_INPOL(1,I)
!
   CALL POLINT(DEL,VEL,N_S,DEL_0,W_N,ER)
   U(H_I_NB(I),H_J_NB(I))= W_N
END DO
END SUBROUTINE
!--------------------------------------------------------------------------------------------------
SUBROUTINE POLINT(XA,YA,NN,XP,YP,DY)
IMPLICIT NONE
INTEGER,PARAMETER         :: MMAX=40
REAL(8),DIMENSION(1:NN)   :: XA,YA
REAL(8),DIMENSION(0:MMAX) :: C,D
INTEGER                   :: NS,MM,NN,IV
REAL(8)                   :: DIF,DIFT,YP,XP,HO,HP,WP,DEN,DY
!
NS = 1
DIF= DABS(XP-XA(1))
DO IV=1,NN
   DIFT=DABS(XP-XA(IV))
   IF (DIFT.LT.DIF) THEN
       NS = IV
       DIF= DIFT
   ENDIF
   C(IV)= YA(IV)
   D(IV)= YA(IV)
END DO
YP= YA(NS)
NS= NS-1
DO MM=1,NN-1
   DO IV=1,NN-MM
      HO = XA(IV)-XP
      HP = XA(IV+MM)-XP
      WP = C(IV+1)-D(IV)
      DEN= HO-HP
      IF(DEN.EQ.0.)PAUSE
      DEN  = WP/DEN
      D(IV)= HP*DEN
      C(IV)= HO*DEN
   END DO
   IF (2*NS.LT.NN-MM)THEN
       DY= C(NS+1)
   ELSE
       DY= D(NS)
       NS= NS-1
   ENDIF
   YP= YP+DY
END DO
RETURN
END
