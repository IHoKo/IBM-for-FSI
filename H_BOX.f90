!-------------- In order to identify the forcing points, first we consider the packet of 3*3 points
!aroud every physical boundary points.-------------------------------------------------------------
!-------------- This subroutine finds the points belonging to the mentioned packet ----------------
SUBROUTINE BOX(X_P,Y_P)
USE IB_VARIABLES, ONLY: H_M_START,H_N_START,H_M_END,H_DEL,                                         &
                        H_N_END,H_DELTA_X,Y_GRID
USE SOLVER_VARIABLES, ONLY: TV,NY
IMPLICIT NONE
REAL(8):: X_P,Y_P
INTEGER:: N
!--------------------------------------------------------------------------------------------------
H_M_START= NINT((X_P)/H_DELTA_X)-1
N        = INT(TV/H_DEL)
H_N_START= NINT((Y_P-Y_GRID(NY/2-N))/H_DEL)+NY/2-N-1
H_M_END  = H_M_START+2
H_N_END  = H_N_START+2
!
END SUBROUTINE

