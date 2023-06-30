! The forcing points, which are used to implement the boundary conditions, are identified in this subroutine
!--------------------------------------------------------------------------------------------------
SUBROUTINE H_FNB(H_NX,H_NY,X,Y,H_X_INTER,H_Y_INTER,H_S_INTER,H_LANDA_X,H_LANDA_Y,                  &
                 H_SLOPE,H_DISTANT,H_DELTA,H_B,FLAG,                                               &
                 H_N_NB,H_I_NB,H_J_NB,H_X_NB,H_Y_NB,H_LAN_X,H_LAN_Y,H_B_NB,                        &
                 H_SLOPE_NB,DIS_1,X_INT,Y_INT)
!
USE IB_VARIABLES, ONLY: H_M_START,H_N_START,H_M_END,H_N_END
IMPLICIT NONE
INTEGER                             :: H_NX,H_NY
INTEGER                             :: I,J
INTEGER                             :: H_N_NB
REAL(8),DIMENSION(0:H_NX+2)         :: X
REAL(8),DIMENSION(0:H_NY+2)         :: Y
INTEGER,DIMENSION(0:H_NX+2,0:H_NY+2):: FLAG
REAL(8),DIMENSION(0:H_NX+2,0:H_NY+2):: H_S_INTER,H_X_INTER,H_Y_INTER,H_LANDA_X,H_LANDA_Y,          &
                                       H_SLOPE,H_DISTANT,H_DELTA,H_B
REAL(8),DIMENSION(1:H_NX)           :: DIS_1,H_LAN_X,H_LAN_Y,H_X_NB,H_Y_NB,H_B_NB,                 &
                                       H_SLOPE_NB,X_INT,Y_INT
INTEGER,DIMENSION(1:H_NX)           :: H_I_NB,H_J_NB
!------------------------------------------------
DO I= H_M_START,H_M_END
   DO J= H_N_START,H_N_END
      IF (DABS(H_DELTA(I,J)-1.0D0)<0.0001   .AND.                                                 &
          DABS(H_DELTA(I+1,J)+1.0D0)<0.0001 .AND.  I+1<=H_M_END.OR.                               &
          DABS(H_DELTA(I,J)-1.0D0)<0.0001   .AND.                                                 &
          DABS(H_DELTA(I-1,J)+1.0D0)<0.0001 .AND.I-1>=H_M_START.OR.                               &
          DABS(H_DELTA(I,J)-1.0D0)<0.0001   .AND.                                                 &
          DABS(H_DELTA(I,J+1)+1.0D0)<0.0001 .AND.  J+1<=H_N_END.OR.                               &
          DABS(H_DELTA(I,J)-1.0D0)<0.0001   .AND.                                                 &
          DABS(H_DELTA(I,J-1)+1.0D0)<0.0001 .AND.J-1>=H_N_START) THEN
          IF (FLAG(I,J)==100 .AND.H_S_INTER(I,J)/=1000.0D0) THEN                 !
              H_N_NB            = H_N_NB+1
              H_I_NB(H_N_NB)    = I
              H_J_NB(H_N_NB)    = J
              FLAG(I,J)         = 0
              H_X_NB(H_N_NB)    = X(I)
              H_Y_NB(H_N_NB)    = Y(J)
              H_LAN_X(H_N_NB)   =-H_LANDA_X(I,J)
              H_LAN_Y(H_N_NB)   =-H_LANDA_Y(I,J)
              H_B_NB(H_N_NB)    = H_B(I,J)
              H_SLOPE_NB(H_N_NB)= H_SLOPE(I,J)
              DIS_1(H_N_NB)     = H_DISTANT(I,J)
              X_INT(H_N_NB)     = H_X_INTER(I,J)
              Y_INT(H_N_NB)     = H_Y_INTER(I,J)
          END IF
      END IF
   END DO
END DO
END SUBROUTINE
