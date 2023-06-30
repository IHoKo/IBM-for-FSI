!------- In this subroutine, the Lagrangian points that construct the geometry is defined ---------
!-------------------------- It has to be considered that the distance between every two neighboring
!points has to be equal to the local grid size!----------------------------------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE GEOMETRY_CYLINDER
USE IB_VARIABLES, ONLY: H_RADIUS,H_X_CENTER,H_Y_CENTER,H_X_PHYSICAL,                               &
                        H_Y_PHYSICAL,H_N_PHYSICAL,H_DELTA_X,H_RADIUS
USE SOLVER_VARIABLES, ONLY: NX
IMPLICIT NONE
REAL(8),DIMENSION(1:NX):: TETA,DELTA_S
INTEGER                :: I
REAL(8)                :: PI
!------------------------ DEFINING THE GEOMETRY OF A CIRCULAR CYLINDER ----------------------------
PI             = ACOS(-1.0D0)
DELTA_S(1)     = 0.0D0
H_X_PHYSICAL(1)= H_RADIUS+H_X_CENTER
H_Y_PHYSICAL(1)= H_Y_CENTER
H_N_PHYSICAL   = 1
!
DO I=2,NX
   DELTA_S(I)     = DELTA_S(I-1)+H_DELTA_X
   TETA(I)        = DELTA_S(I)/(H_RADIUS)
   IF (TETA(I)<=(PI/2.0D0)) THEN
       H_N_PHYSICAL      = H_N_PHYSICAL+1
       H_X_PHYSICAL(I)= H_RADIUS*COS(TETA(I))+H_X_PHYSICAL(1)-H_RADIUS
       H_Y_PHYSICAL(I)= H_RADIUS*SIN(TETA(I))+H_Y_PHYSICAL(1)
   END IF
END DO
!
DO I=H_N_PHYSICAL,1,-1
   H_X_PHYSICAL(2*H_N_PHYSICAL-I+1)= H_X_CENTER-(H_X_PHYSICAL(I)-H_X_CENTER)
   H_Y_PHYSICAL(2*H_N_PHYSICAL-I+1)= H_Y_PHYSICAL(I)
END DO
H_N_PHYSICAL= 2*H_N_PHYSICAL
!
DO I=H_N_PHYSICAL,3,-1
   H_X_PHYSICAL(2*H_N_PHYSICAL-I+1)= H_X_PHYSICAL(I-1)
   H_Y_PHYSICAL(2*H_N_PHYSICAL-I+1)= -H_Y_PHYSICAL(I-1)
END DO
H_N_PHYSICAL= 2*H_N_PHYSICAL-2
WRITE(*,*)H_N_PHYSICAL
END SUBROUTINE
