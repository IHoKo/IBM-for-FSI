!------------------------------------- INTERFACE DESCRIPTION --------------------------------------
!In order to define the boundary, the interface tracking algorithm employes piecewise polynomial
!fits to compute the necessary shape derivatives at each location. This subroutine is used to calculate
!the coefficients of the polynomial that is constructed for every three neighboring points.
SUBROUTINE COEFFICIENT (K_P)
USE IB_VARIABLES, ONLY:H_X_POINT,H_Y_POINT,H_A_X,H_A_Y
IMPLICIT NONE
REAL(8) :: ARC_L1,ARC_L2
REAL(8) :: ALPHA,BETA,MU,ETHA
INTEGER :: K_P
!--------------------------------------------------------------------------------------------------
ARC_L1  = DSQRT((H_X_POINT(1,K_P)-H_X_POINT(0,K_P))*                                               &
                (H_X_POINT(1,K_P)-H_X_POINT(0,K_P))+                                               &
                (H_Y_POINT(1,K_P)-H_Y_POINT(0,K_P))*                                               &
                (H_Y_POINT(1,K_P)-H_Y_POINT(0,K_P)))
ARC_L2  = DSQRT((H_X_POINT(2,K_P)-H_X_POINT(0,K_P))*                                               &
                (H_X_POINT(2,K_P)-H_X_POINT(0,K_P))+                                               &
                (H_Y_POINT(2,K_P)-H_Y_POINT(0,K_P))*                                               &
                       (H_Y_POINT(2,K_P)-H_Y_POINT(0,K_P)))
ALPHA   = -ARC_L2*ARC_L2/(ARC_L1*ARC_L1)
BETA    = -ARC_L1/(ALPHA*ARC_L1+ARC_L2)
MU      = ALPHA*(H_X_POINT(1,K_P)-H_X_POINT(0,K_P))+                                               &
                (H_X_POINT(2,K_P)-H_X_POINT(0,K_P))
ETHA    = ALPHA*(H_Y_POINT(1,K_P)-H_Y_POINT(0,K_P))+                                               &
                (H_Y_POINT(2,K_P)-H_Y_POINT(0,K_P))
H_A_X(0)=((H_X_POINT(1,K_P)-H_X_POINT(0,K_P))+BETA*MU)/(ARC_L1*ARC_L1)
H_A_X(1)=  MU/(ALPHA*ARC_L1+ARC_L2)
H_A_X(2)=  H_X_POINT(0,K_P)
H_A_Y(0)=((H_Y_POINT(1,K_P)-H_Y_POINT(0,K_P))+BETA*ETHA)/(ARC_L1*ARC_L1)
H_A_Y(1)=  ETHA/(ALPHA*ARC_L1+ARC_L2)
H_A_Y(2)=  H_Y_POINT(0,K_P)
!
END SUBROUTINE
