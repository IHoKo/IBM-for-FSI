!This subroutine basically calls other subroutines that help identifying the numerical boundary points.
!Also, the subroutines that find th datapoints in interpolation stencil are called here.
!The last step of implementing boundary conditions which is interpolation of the velocity of
!the numerical boundary points is performed by 'h_velocity' later.
SUBROUTINE IMPBC
USE IB_VARIABLES
USE SOLVER_VARIABLES, ONLY: NX,NY
IMPLICIT NONE
REAL(8):: X_BOX,Y_BOX,PI
INTEGER:: K
!----------------------------------------- INITIALIZATION -----------------------------------------
H_N_PHYSICAL= 0
HU_N_NB     = 0
HV_N_NB     = 0
HT_N_NB     = 0
PI          = ACOS(-1.0D0)
!
CALL H_GRID
!
UFLAG = 100
VFLAG = 100
TFLAG = 100
!---------------------------------------- CYLINDER GEOMETRY ---------------------------------------
CALL GEOMETRY_CYLINDER
!
DO K=1,H_N_PHYSICAL
   IF (K==1) THEN
       H_X_POINT(0,K)= H_X_PHYSICAL(H_N_PHYSICAL)
       H_X_POINT(1,K)= H_X_PHYSICAL(K)
       H_X_POINT(2,K)= H_X_PHYSICAL(K+1)
       H_Y_POINT(0,K)= H_Y_PHYSICAL(H_N_PHYSICAL)
       H_Y_POINT(1,K)= H_Y_PHYSICAL(K)
       H_Y_POINT(2,K)= H_Y_PHYSICAL(K+1)
    ELSEIF (K==H_N_PHYSICAL) THEN
            H_X_POINT(0,K)= H_X_PHYSICAL(H_N_PHYSICAL-1)
            H_X_POINT(1,K)= H_X_PHYSICAL(H_N_PHYSICAL)
            H_X_POINT(2,K)= H_X_PHYSICAL(1)
            H_Y_POINT(0,K)= H_Y_PHYSICAL(H_N_PHYSICAL-1)
            H_Y_POINT(1,K)= H_Y_PHYSICAL(H_N_PHYSICAL)
            H_Y_POINT(2,K)= H_Y_PHYSICAL(1)
    ELSE
       H_X_POINT(0,K)= H_X_PHYSICAL(K-1)
       H_X_POINT(1,K)= H_X_PHYSICAL(K)
       H_X_POINT(2,K)= H_X_PHYSICAL(K+1)
       H_Y_POINT(0,K)= H_Y_PHYSICAL(K-1)
       H_Y_POINT(1,K)= H_Y_PHYSICAL(K)
       H_Y_POINT(2,K)= H_Y_PHYSICAL(K+1)
    END IF
!---------- FINDING THE COEFFICIENT OF THE INTERFACIAL MARKERS IN ARCLENGTH COORDINATES -----------
    CALL COEFFICIENT (K)
!---------- FINDING THE 3*3 PACKAGES OF THE NODES AROUND EVERY PHYSICAL BOUNDARY POINTS -----------
    X_BOX= H_X_PHYSICAL(K)
    Y_BOX= H_Y_PHYSICAL(K)
!
    CALL BOX(X_BOX,Y_BOX)
!-- WE HAVE TO FINED THE NUMERICAL BOUNDARY POINTS FOR X AND Y COMPONENTS OF VELOCITY SEPARATELY --
!---------------------------------- X-COMPONENT OF VELOCITY ---------------------------------------
!--------- FINDING THE INTERSECTION BETWEEN THE NORMAL LINE TO THE PHYSICAL BOUNDARY --------------
!--------------- FROM THE EVERY POINT IN OUT 3*3 PACKAGE AND PHYSICAL BOUNDARY  -------------------
    CALL H_INTERSECTION (NX,NY,X12_GRID,Y_GRID,HU_S_INTER)
!------------ DELTA IS THE INNER PRODUCT OF THE UNIT NORMAL VECTOR OF THE PHYSICAL ----------------
!----------------- BOUNDARY AND THE UNIT VECTOR STARTED FROM THE POINTS IN OUR --------------------
!------------------------- 3*3 PACHAGE NORMAL TO THE PHYSICAL BOUNDARY ----------------------------
    CALL H_INNER_PRODUCT(NX,NY,X12_GRID,Y_GRID,XU_INTER,YU_INTER,HU_S_INTER,HU_LANDA_X,            &
                         HU_LANDA_Y,HU_SLOPE,HU_DISTANT,HU_DELTA,HU_B)
!---------- THE NUMERICAL BOUNDARY POINTS WHO MEET THE REQUIRED CONDITIONS ARE FOUND HERE ---------
    CALL H_FNB(NX,NY,X12_GRID,Y_GRID,XU_INTER,YU_INTER,HU_S_INTER,HU_LANDA_X,HU_LANDA_Y,           &
               HU_SLOPE,HU_DISTANT,HU_DELTA,HU_B,UFLAG,                                            &
               HU_N_NB,HU_I_NB,HU_J_NB,HU_X_NB,HU_Y_NB,HU_LAN_X,HU_LAN_Y,HU_B_NB,                  &
               HU_SLOPE_NB,UDIS_1,XU_INT,YU_INT)
!----------------------------------- Y-COMPONENT OF VELOCITY -------------------------------------
    CALL H_INTERSECTION (NX,NY,X_GRID,Y12_GRID,HV_S_INTER)
    CALL H_INNER_PRODUCT(NX,NY,X_GRID,Y12_GRID,XV_INTER,YV_INTER,HV_S_INTER,HV_LANDA_X,            &
                         HV_LANDA_Y,HV_SLOPE,HV_DISTANT,HV_DELTA,HV_B)
    CALL H_FNB(NX,NY,X_GRID,Y12_GRID,XV_INTER,YV_INTER,HV_S_INTER,HV_LANDA_X,HV_LANDA_Y,           &
               HV_SLOPE,HV_DISTANT,HV_DELTA,HV_B,VFLAG,                                            &
               HV_N_NB,HV_I_NB,HV_J_NB,HV_X_NB,HV_Y_NB,HV_LAN_X,HV_LAN_Y,HV_B_NB,                  &
               HV_SLOPE_NB,VDIS_1,XV_INT,YV_INT)
!------------------------- REPEATED PROCESS FOR PRESSURE OR TEMPRATURE ----------------------------
    CALL H_INTERSECTION (NX,NY,X_GRID,Y_GRID,HT_S_INTER)
    CALL H_INNER_PRODUCT(NX,NY,X_GRID,Y_GRID,X_INTER,Y_INTER,HT_S_INTER,HT_LANDA_X,                &
                         HT_LANDA_Y,HT_SLOPE,HT_DISTANT,HT_DELTA,HT_B)
    CALL H_FNB(NX,NY,X_GRID,Y_GRID,X_INTER,Y_INTER,HT_S_INTER,HT_LANDA_X,HT_LANDA_Y,               &
               HT_SLOPE,HT_DISTANT,HT_DELTA,HT_B,TFLAG,                                            &
               HT_N_NB,HT_I_NB,HT_J_NB,HT_X_NB,HT_Y_NB,HT_LAN_X,HT_LAN_Y,HT_B_NB,                  &
               HT_SLOPE_NB,TDIS_1,XT_INT,YT_INT)
 END DO
!------------ FINDING THE VIRTUAL POINTS IN THE FLUID REGION TO BE USED FOR INTERPOLATION ---------
CALL H_INTERPOLATION_DATA(NX,HU_N_NB,HU_X_NB,HU_Y_NB,HU_SLOPE_NB,HU_B_NB,UDIS_1,HU_LAN_X,          &
                     HU_LAN_Y,HU_DIS_INPOL,IU_IP,JU_IP,HU_I_NB,HU_J_NB,HU_X_INT,HU_Y_INT)
CALL H_INTERPOLATION_DATA(NX,HV_N_NB,HV_X_NB,HV_Y_NB,HV_SLOPE_NB,HV_B_NB,VDIS_1,HV_LAN_X,          &
                     HV_LAN_Y,HV_DIS_INPOL,IV_IP,JV_IP,HV_I_NB,HV_J_NB,HV_X_INT,HV_Y_INT)
CALL H_INTERPOLATION_DATA(NX,HT_N_NB,HT_X_NB,HT_Y_NB,HT_SLOPE_NB,HT_B_NB,TDIS_1,HT_LAN_X,          &
                     HT_LAN_Y,HT_DIS_INPOL,IT_IP,JT_IP,HT_I_NB,HT_J_NB,HT_X_INT,HT_Y_INT)
!
END SUBROUTINE
