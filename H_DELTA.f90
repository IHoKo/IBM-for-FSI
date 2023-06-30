!-------------- In this subroutine the inner product of normal vector of the surface and the vector
!from every point belonging to the 3*3 packet to the surface is calculated ------------------------
SUBROUTINE H_INNER_PRODUCT(H_NX,H_NY,X,Y,H_X_INTER,H_Y_INTER,H_S_INTER,H_LANDA_X,                  &
                           H_LANDA_Y,H_SLOPE,H_DISTANT,H_DELTA,H_B)
USE IB_VARIABLES, ONLY: H_M_START,H_N_START,H_M_END,H_N_END,H_A_X,H_A_Y
IMPLICIT NONE
INTEGER                             :: I,J
INTEGER                             :: H_NX,H_NY
REAL(8)                             :: PI
REAL(8)                             :: NORMAL_X,NORMAL_Y
REAL(8),DIMENSION(0:H_NX+2)         :: X
REAL(8),DIMENSION(0:H_NY+2)         :: Y
REAL(8),DIMENSION(0:H_NX+2,0:H_NY+2):: H_S_INTER,H_X_INTER,H_Y_INTER,H_LANDA_X,H_LANDA_Y,          &
                                       H_SLOPE,H_DISTANT,H_DELTA,H_B
!--------------------------------------------------------------------------------------------------
PI= ACOS(-1.0D0)
DO J=H_N_START,H_N_END
   DO I=H_M_START,H_M_END
      H_X_INTER(I,J)= H_A_X(0)*H_S_INTER(I,J)*H_S_INTER(I,J)+                                      &
                      H_A_X(1)*H_S_INTER(I,J)+H_A_X(2)
      H_Y_INTER(I,J)= H_A_Y(0)*H_S_INTER(I,J)*H_S_INTER(I,J)+                                      &
                      H_A_Y(1)*H_S_INTER(I,J)+H_A_Y(2)
!--------------------------------------------------------------------------------------------------
!    NORMAL VECTOR FROM EACH POINT BELONGING TO THE 3*3PACKAGE TO THE PHISICAL BOUNDARY SURFACE
!--------------------------------------------------------------------------------------------------
      H_LANDA_X(I,J)= (H_X_INTER(I,J)-X(I))/DSQRT((H_X_INTER(I,J)-X(I))*(H_X_INTER(I,J)-           &
                       X(I))+(H_Y_INTER(I,J)-Y(J))*(H_Y_INTER(I,J)-Y(J)))
      H_LANDA_Y(I,J)= (H_Y_INTER(I,J)-Y(J))/DSQRT((H_X_INTER(I,J)-X(I))*(H_X_INTER(I,J)-           &
                       X(I))+(H_Y_INTER(I,J)-Y(J))*(H_Y_INTER(I,J)-Y(J)))
!--------------------------------------------------------------------------------------------------
!                                  THE SLOPE OF THE NORMAL VECTOR
!--------------------------------------------------------------------------------------------------
      H_SLOPE(I,J)= -(2.0D0*H_A_X(0)*H_S_INTER(I,J)+H_A_X(1))/                                     &
                     (2.0D0*H_A_Y(0)*H_S_INTER(I,J)+H_A_Y(1))
      IF ((2.0D0*H_A_Y(0)*H_S_INTER(I,J)+H_A_Y(1))==0) THEN
           H_SLOPE(I,J)= TAN(PI/2.0D0)
      END IF
!--------------------------------------------------------------------------------------------------
!                         NORMAL VECTOR OF THE PHYSICAL BOUNDARY SURFACE
!--------------------------------------------------------------------------------------------------
      NORMAL_X= -(2.0D0*H_A_Y(0)*H_S_INTER(I,J)+H_A_Y(1))/DSQRT                                    &
                ((2.0D0*H_A_Y(0)*H_S_INTER(I,J)+H_A_Y(1))                                          &
                *(2.0D0*H_A_Y(0)*H_S_INTER(I,J)+H_A_Y(1))+                                         &
                 (2.0D0*H_A_X(0)*H_S_INTER(I,J)+H_A_X(1))*                                         &
                 (2.0D0*H_A_X(0)*H_S_INTER(I,J)+H_A_X(1)))
!
      NORMAL_Y=  (2.0D0*H_A_X(0)*H_S_INTER(I,J)+H_A_X(1))/DSQRT                                    &
                ((2.0D0*H_A_Y(0)*H_S_INTER(I,J)+H_A_Y(1))                                          &
                *(2.0D0*H_A_Y(0)*H_S_INTER(I,J)+H_A_Y(1))+                                         &
                 (2.0D0*H_A_X(0)*H_S_INTER(I,J)+H_A_X(1))*                                         &
                 (2.0D0*H_A_X(0)*H_S_INTER(I,J)+H_A_X(1)))
!
      H_DISTANT(I,J)= DSQRT((H_X_INTER(I,J)-X(I))*(H_X_INTER(I,J)-X(I))+                           &
                            (H_Y_INTER(I,J)-Y(J))*(H_Y_INTER(I,J)-Y(J)))
      H_DELTA(I,J)  = NORMAL_X*H_LANDA_X(I,J)+NORMAL_Y*H_LANDA_Y(I,J)
      H_B(I,J)      = Y(J)+(2.0D0*H_A_X(0)*H_S_INTER(I,J)+H_A_X(1))/                               &
                      (2.0D0*H_A_Y(0)*H_S_INTER(I,J)+H_A_Y(1))*X(I)
   END DO
END DO
END SUBROUTINE
