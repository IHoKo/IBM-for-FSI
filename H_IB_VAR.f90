!The variables that are used in implementation of the boundary onditions are defined in this madule.
MODULE IB_VARIABLES
USE SOLVER_VARIABLES, ONLY: NX,NY
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------
INTEGER,PARAMETER               :: N_ORDER= 2
!----------------------------- REQUIRED PARAMETERS TO DEFINE GEOMETRY -----------------------------
REAL(8),DIMENSION(1:NX)         :: H_X_PHYSICAL,H_Y_PHYSICAL
INTEGER                         :: H_N_PHYSICAL
!--------------------------------------------- GRID -----------------------------------------------
REAL(8),DIMENSION(0:NY+2)       :: Y12_GRID,Y_GRID
REAL(8),DIMENSION(0:NX+2)       :: X12_GRID,X_GRID
REAL(8)                         :: H_DELTA_X,H_DEL
REAL(8),DIMENSION(0:NY+2)       :: H_DELTA_Y12,H_DELTA_Y
INTEGER,DIMENSION(0:NX+2,0:NY+2):: I_GRID,I12_GRID
INTEGER,DIMENSION(0:NX+2,0:NY+2):: J_GRID,J12_GRID
!------------------------------------------- GEOMETRY ---------------------------------------------
REAL(8)                         :: H_RADIUS
REAL(8)                         :: H_X_CENTER,H_Y_CENTER
!--------------------------------------- INITIAL CONDITION ----------------------------------------
REAL(8),DIMENSION(1:NX)         :: H_U_BODY
REAL(8),DIMENSION(1:NX)         :: H_V_BODY
REAL(8),DIMENSION(1:NX)         :: H_T_BODY
!------------------------------------------ COEFICIENTS -------------------------------------------
REAL(8),DIMENSION(0:2)          :: H_A_X,H_A_Y
REAL(8),DIMENSION(0:2,1:NX)     :: H_X_POINT,H_Y_POINT
!--------------------------------------------- BOX ------------------------------------------------
INTEGER                         :: H_M_START,H_N_START,H_M_END,H_N_END
!-------------------------------------------- DELTA -----------------------------------------------
REAL(8),DIMENSION(0:NX+2,0:NY+2):: HV_LANDA_X,HV_LANDA_Y,HV_DISTANT,                               &
                                   HV_DELTA,HV_SLOPE,HV_B,HV_S_INTER
REAL(8),DIMENSION(0:NX+2,0:NY+2):: HU_LANDA_X,HU_LANDA_Y,HU_DISTANT,                               &
                                   HU_DELTA,HU_SLOPE,HU_B,HU_S_INTER
REAL(8),DIMENSION(0:NX+2,0:NY+2):: XU_INTER,YU_INTER,XV_INTER,YV_INTER
REAL(8),DIMENSION(0:NX+2,0:NY+2):: X_INTER,Y_INTER
REAL(8),DIMENSION(0:NX+2,0:NY+2):: HT_LANDA_X,HT_LANDA_Y,HT_DISTANT,                               &
                                   HT_DELTA,HT_SLOPE,HT_B,HT_S_INTER
!----------------------------------- FIND NUMERICAL BOUNDARY --------------------------------------
INTEGER,DIMENSION(0:NX+2,0:NY+2):: UFLAG,VFLAG,TFLAG
REAL(8),DIMENSION(1:NX)         :: HU_LAN_X,HU_LAN_Y,HU_B_NB,HU_SLOPE_NB
REAL(8),DIMENSION(1:NX)         :: HV_LAN_X,HV_LAN_Y,HV_B_NB,HV_SLOPE_NB
REAL(8),DIMENSION(1:NX)         :: HT_LAN_X,HT_LAN_Y,HT_B_NB,HT_SLOPE_NB
REAL(8),DIMENSION(1:NX)         :: XU_INT,YU_INT,XV_INT,YV_INT
REAL(8),DIMENSION(1:NX)         :: XT_INT,YT_INT
REAL(8),DIMENSION(1:NX)         :: HU_X_NB,HU_Y_NB
REAL(8),DIMENSION(1:NX)         :: HV_X_NB,HV_Y_NB
REAL(8),DIMENSION(1:NX)         :: HT_X_NB,HT_Y_NB
REAL(8),DIMENSION(1:NX)         :: UDIS_1,VDIS_1,TDIS_1
INTEGER,DIMENSION(1:NX)         :: HV_I_NB,HV_J_NB
INTEGER,DIMENSION(1:NX)         :: HU_I_NB,HU_J_NB
INTEGER,DIMENSION(1:NX)         :: HT_I_NB,HT_J_NB
INTEGER                         :: HU_N_NB,HV_N_NB,HT_N_NB
!---------------------------------------- INTERPOLATION -------------------------------------------
REAL(8),DIMENSION(1:2,1:NX)     :: HU_X_INT,HU_Y_INT
REAL(8),DIMENSION(1:2,1:NX)     :: HV_X_INT,HV_Y_INT
REAL(8),DIMENSION(1:2,1:NX)     :: HT_X_INT,HT_Y_INT
REAL(8),DIMENSION(1:3,1:NX)     :: HU_DIS_INPOL,HV_DIS_INPOL,HT_DIS_INPOL
INTEGER,DIMENSION(1:4,1:2,1:NX) :: IU_IP,JU_IP
INTEGER,DIMENSION(1:4,1:2,1:NX) :: IV_IP,JV_IP
INTEGER,DIMENSION(1:4,1:2,1:NX) :: IT_IP,JT_IP
!--------------------------------------------------------------------------------------------------
REAL(8),DIMENSION(0:NX+2,0:NY+2):: H_U,H_V,H_T
!-------------------------- FIND NUMERICAL BOUNDARY TO CALCULATE FORCE ----------------------------
REAL(8)                         :: CD,CL,CX,CY
INTEGER,DIMENSION(0:NX+2,0:NY+2):: FLAG_RBF
END MODULE
!--------------------------------------------------------------------------------------------------
