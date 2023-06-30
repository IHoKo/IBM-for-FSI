!-------------- In this subroutine, the calculated velocity for forcinf points is replaced with the
!------------------------------------------------------------ value obtained by solution of N-S eq.
!--------------------------------------------------------------------------------------------------
SUBROUTINE BC_INTERPOLATION
USE IB_VARIABLES, ONLY: HU_N_NB,HV_N_NB,HT_N_NB,H_U,H_V,H_T,H_U_BODY,H_V_BODY,H_T_BODY,            &
                        H_X_CENTER,H_Y_CENTER,H_RADIUS,X_GRID,Y_GRID,X12_GRID,Y12_GRID,            &
                        HU_I_NB,HU_J_NB,HV_I_NB,HV_J_NB,HT_I_NB,HT_J_NB,HU_X_INT,HU_Y_INT,         &
                        HV_X_INT,HV_Y_INT,HT_X_INT,HT_Y_INT,HU_DIS_INPOL,HV_DIS_INPOL,            &
                        HT_DIS_INPOL,IU_IP,JU_IP,IV_IP,JV_IP,IT_IP,JT_IP

USE SOLVER_VARIABLES, ONLY: UX,UY,NX,NY,NZ,Y12,X12,Y,X,TT,UC
IMPLICIT NONE
INTEGER:: I,J,K
REAL(8):: R1,R2,R
!--------------------------------------------------------------------------------------------------
DO K=0,NZ+1
   DO J=0,NY+1
      DO I=0,NX+1
         R1= DSQRT((X12(I)-H_X_CENTER)*(X12(I)-H_X_CENTER)+                                        &
                   (Y(J)-H_Y_CENTER)*(Y(J)-H_Y_CENTER))
         R2= DSQRT((X(I)-H_X_CENTER)*(X(I)-H_X_CENTER)+                                            &
                   (Y12(J)-H_Y_CENTER)*(Y12(J)-H_Y_CENTER))
         IF (R1<=H_RADIUS) THEN
             UX(I,J,K)= UC
         END IF
         IF (R2<=H_RADIUS) THEN
             UY(I,J,K)= 0.0D0
         END IF
      END DO
   END DO
END DO
!------------------------------------------ TEMPRATURE --------------------------------------------
DO K=0,NZ+1
   DO J=0,NY+1
      DO I=0,NX+1
         R= DSQRT((X(I)-H_X_CENTER)*(X(I)-H_X_CENTER)+                                             &
                  (Y(J)-H_Y_CENTER)*(Y(J)-H_Y_CENTER))
         IF (R<=H_RADIUS) THEN
             TT(I,J,K)= 1.0D0
         END IF
      END DO
   END DO
END DO
!--------------------------------------------------------------------------------------------------
DO K=0,NZ+1
   DO J=0,NY+1
      DO I=0,NX+1
         H_U(I,J)= UX(I,J,K)
         H_V(I,J)= UY(I,J,K)
         H_T(I,J)= TT(I,J,K)
       END DO
    END DO
!------------------------------------------- CYLINDER ---------------------------------------------
    DO I=1,HU_N_NB
       H_U_BODY(I)= UC
    END DO
    DO I=1,HV_N_NB
       H_V_BODY(I)= 0.0D0
    END DO
!
    DO I=1,HT_N_NB
       H_T_BODY(I)= 1.0D0
    END DO
!--------------------------------------------------------------------------------------------------
    CALL H_VELOCITY(NX,NY,HU_N_NB,IU_IP,JU_IP,X12_GRID,Y_GRID,HU_I_NB,HU_J_NB,HU_X_INT,HU_Y_INT,HU_DIS_INPOL,H_U_BODY,H_U)
    CALL H_VELOCITY(NX,NY,HV_N_NB,IV_IP,JV_IP,X_GRID,Y12_GRID,HV_I_NB,HV_J_NB,HV_X_INT,HV_Y_INT,HV_DIS_INPOL,H_V_BODY,H_V)
    CALL H_VELOCITY(NX,NY,HT_N_NB,IT_IP,JT_IP,X_GRID,Y_GRID,HT_I_NB,HT_J_NB,HT_X_INT,HT_Y_INT,HT_DIS_INPOL,H_T_BODY,H_T)
!--------------------------------------------------------------------------------------------------
    DO J=0,NY+1
       DO I=0,NX+1
          UX(I,J,K)= H_U(I,J)
          UY(I,J,K)= H_V(I,J)
          TT(I,J,K)= H_T(I,J)
       END DO
    END DO
 END DO
END SUBROUTINE
