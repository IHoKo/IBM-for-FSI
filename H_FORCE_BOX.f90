!With this subroutine, we can calculate the lift and drag coefficients.
!Note thta the length characteristics is assumed to be '1' here and
!characteristics velocity, UM, density, RO, and dynamic viscosity,MU have to be determined
!For further information about the method of calculating forces, onemay refer to
!'A moving control volume approach to computing hydrodynamic forces and torques on immersed bodies',
!Nangia et al.
!Journal of Computational Physics, Volume 347, 15 October 2017, Pages 437-462
SUBROUTINE H_FORCE_BOX
USE IB_VARIABLES, ONLY: H_DELTA_X,H_RADIUS,CX,CY
USE SOLVER_VARIABLES, ONLY: NX,NY,X12,Y12,P,UX,UY,UC,UC_OLD,UM,DT,UX_OLD,UY_OLD
IMPLICIT NONE
INTEGER:: I,J,K
INTEGER:: H_MS_BOX,H_ME_BOX,H_NS_BOX,H_NE_BOX
REAL(8):: H_L_BOX,H_W_BOX,H_LC_BOX,H_WC_BOX
REAL(8):: RO,MU,A,PI
REAL(8):: H_SUM_BY_P,H_SUM_TY_P,H_SUM_LX_P,H_SUM_RX_P
REAL(8):: H_SUM_BY_UROU,H_SUM_BX_UROU,H_SUM_TX_UROU,H_SUM_TY_UROU
REAL(8):: H_SUM_LY_UROU,H_SUM_LX_UROU,H_SUM_RX_UROU,H_SUM_RY_UROU
REAL(8):: H_SUM_BY_MU,H_SUM_BX_MU,H_SUM_TX_MU,H_SUM_TY_MU
REAL(8):: H_SUM_LY_MU,H_SUM_LX_MU,H_SUM_RX_MU,H_SUM_RY_MU
REAL(8):: H_SUM_X_ROU,H_SUM_Y_ROU,H_SUM_X_ROU_OLD,H_SUM_Y_ROU_OLD
REAL(8):: BODY_MOMENTUM
REAL(8):: FX,FY
!----------------------------------------- Initialization -----------------------------------------
PI             = DACOS(-1.0D0)
H_L_BOX        = 8.0D0
H_W_BOX        = 4.0D0
H_LC_BOX       = X12(NX/2)
H_WC_BOX       = Y12(NY/2)
H_MS_BOX       = INT((H_LC_BOX-0.5D0*H_L_BOX)/H_DELTA_X)
H_ME_BOX       = INT(H_L_BOX/H_DELTA_X)+H_MS_BOX
H_NS_BOX       = NY/2-INT(DABS(H_WC_BOX-0.5D0*H_W_BOX)/H_DELTA_X)
H_NE_BOX       = NY/2+INT(DABS(H_WC_BOX-0.5D0*H_W_BOX)/H_DELTA_X)
!------------------------------------------------
K              = 1
H_SUM_BY_P     = 0.0D0
H_SUM_TY_P     = 0.0D0
H_SUM_BY_UROU  = 0.0D0
H_SUM_TY_UROU  = 0.0D0
H_SUM_LY_UROU  = 0.0D0
H_SUM_RY_UROU  = 0.0D0
H_SUM_BY_MU    = 0.0D0
H_SUM_TY_MU    = 0.0D0
H_SUM_LY_MU    = 0.0D0
H_SUM_RY_MU    = 0.0D0
H_SUM_Y_ROU    = 0.0D0
H_SUM_LX_P     = 0.0D0
H_SUM_RX_P     = 0.0D0
H_SUM_BX_UROU  = 0.0D0
H_SUM_TX_UROU  = 0.0D0
H_SUM_LX_UROU  = 0.0D0
H_SUM_RX_UROU  = 0.0D0
H_SUM_BX_MU    = 0.0D0
H_SUM_TX_MU    = 0.0D0
H_SUM_LX_MU    = 0.0D0
H_SUM_RX_MU    = 0.0D0
H_SUM_X_ROU    = 0.0D0
H_SUM_X_ROU_OLD= 0.0D0
H_SUM_Y_ROU_OLD= 0.0D0
RO             = 1.0D0
MU             = 0.01D0
!------------------------------------ NUMERICAL DISCRITIZATION ------------------------------------
!--------------------------- DISCRETE APPRPXIMATION TO SURFACE INTEGRALS --------------------------
DO I= H_MS_BOX,H_ME_BOX
   H_SUM_BY_P= 0.5D0*(P(I,H_NS_BOX,K)+P(I,H_NS_BOX+1,K))*H_DELTA_X+H_SUM_BY_P
   H_SUM_TY_P=-0.5D0*(P(I,H_NE_BOX,K)+P(I,H_NE_BOX+1,K))*H_DELTA_X+H_SUM_TY_P
END DO
!
DO I= H_NS_BOX,H_NE_BOX
   H_SUM_LX_P= 0.5D0*(P(H_MS_BOX,I,K)+P(H_MS_BOX+1,I,K))*H_DELTA_X+H_SUM_LX_P
   H_SUM_RX_P=-0.5D0*(P(H_ME_BOX,I,K)+P(H_ME_BOX+1,I,K))*H_DELTA_X+H_SUM_RX_P
END DO
!------------------------------------------------
DO I= H_MS_BOX,H_ME_BOX
   H_SUM_BY_UROU= RO*(UY(I,H_NS_BOX,K)*UY(I,H_NS_BOX,K))*H_DELTA_X+H_SUM_BY_UROU
   H_SUM_BX_UROU= RO*(UY(I,H_NS_BOX,K)*(0.25D0*(UX(I-1,H_NS_BOX,K)+UX(I,H_NS_BOX,K)+               &
                      UX(I-1,H_NS_BOX+1,K)+UX(I,H_NS_BOX+1,K))))*H_DELTA_X+H_SUM_BX_UROU
!
   H_SUM_TY_UROU=-RO*(UY(I,H_NE_BOX,K)*UY(I,H_NE_BOX,K))*H_DELTA_X+H_SUM_TY_UROU
   H_SUM_TX_UROU=-RO*(UY(I,H_NE_BOX,K)*(0.25D0*(UX(I-1,H_NE_BOX,K)+UX(I,H_NE_BOX,K)+               &
                      UX(I-1,H_NE_BOX+1,K)+UX(I,H_NE_BOX+1,K))))*H_DELTA_X+H_SUM_TX_UROU
END DO
!
DO I= H_NS_BOX,H_NE_BOX
   H_SUM_LX_UROU= RO*(UX(H_MS_BOX,I,K)*UX(H_MS_BOX,I,K))*H_DELTA_X+H_SUM_LX_UROU
   H_SUM_LY_UROU= RO*(UX(H_MS_BOX,I,K)*(0.25D0*(UY(H_MS_BOX,I-1,K)+UY(H_MS_BOX,I,K)+               &
                      UY(H_MS_BOX+1,I-1,K)+UY(H_MS_BOX+1,I,K))))*H_DELTA_X+H_SUM_LY_UROU
!
   H_SUM_RX_UROU=-RO*(UX(H_ME_BOX,I,K)*UX(H_ME_BOX,I,K))*H_DELTA_X+H_SUM_RX_UROU
   H_SUM_RY_UROU=-RO*(UX(H_ME_BOX,I,K)*(0.25D0*(UY(H_ME_BOX,I-1,K)+UY(H_ME_BOX,I,K)+               &
                      UY(H_ME_BOX+1,I-1,K)+UY(H_ME_BOX+1,I,K))))*H_DELTA_X+H_SUM_RY_UROU
END DO
!------------------------------------------------
DO I= H_MS_BOX,H_ME_BOX
   H_SUM_BX_MU=-MU*((UX(I,H_NS_BOX+1,K)-UX(I,H_NS_BOX,K)+UX(I-1,H_NS_BOX+1,K)-                     &
                     UX(I-1,H_NS_BOX,K))/(2.0D0*H_DELTA_X)+(UY(I+1,H_NS_BOX,K)-                    &
                     UY(I-1,H_NS_BOX,K))/(2.0D0*H_DELTA_X))*H_DELTA_X+H_SUM_BX_MU
   H_SUM_BY_MU=-2.0D0*MU*((UY(I,H_NS_BOX+1,K)-UY(I,H_NS_BOX-1,K))                                  &
              /(2.0D0*H_DELTA_X))*H_DELTA_X+H_SUM_BY_MU
!
   H_SUM_TX_MU= MU*((UX(I,H_NE_BOX+1,K)-UX(I,H_NE_BOX,K)+UX(I-1,H_NE_BOX+1,K)-                     &
                     UX(I-1,H_NE_BOX,K))/(2.0D0*H_DELTA_X)+(UY(I+1,H_NE_BOX,K)-                    &
                     UY(I-1,H_NE_BOX,K))/(2.0D0*H_DELTA_X))*H_DELTA_X+H_SUM_TX_MU
   H_SUM_TY_MU= 2.0D0*MU*((UY(I,H_NE_BOX+1,K)-UY(I,H_NE_BOX-1,K))                                  &
              /(2.0D0*H_DELTA_X))*H_DELTA_X+H_SUM_TY_MU
END DO
!
DO I= H_NS_BOX,H_NE_BOX
   H_SUM_LY_MU=-MU*((UY(H_MS_BOX+1,I,K)-UY(H_MS_BOX,I,K)+UY(H_MS_BOX+1,I-1,K)-                     &
                     UY(H_MS_BOX,I-1,K))/(2.0D0*H_DELTA_X)+(UX(H_MS_BOX,I+1,K)-                    &
                     UX(H_MS_BOX,I-1,K))/(2.0D0*H_DELTA_X))*H_DELTA_X+H_SUM_LY_MU
   H_SUM_LX_MU=-2.0D0*MU*((UX(H_MS_BOX+1,I,K)-UX(H_MS_BOX-1,I,K))                                  &
              /(2.0D0*H_DELTA_X))*H_DELTA_X+H_SUM_LX_MU
!
   H_SUM_RY_MU= MU*((UY(H_ME_BOX+1,I,K)-UY(H_ME_BOX,I,K)+UY(H_ME_BOX+1,I-1,K)-                     &
                     UY(H_ME_BOX,I-1,K))/(2.0D0*H_DELTA_X)+(UX(H_ME_BOX,I+1,K)-                    &
                     UX(H_ME_BOX,I-1,K))/(2.0D0*H_DELTA_X))*H_DELTA_X+H_SUM_RY_MU
   H_SUM_RX_MU= 2.0D0*MU*((UX(H_ME_BOX+1,I,K)-UX(H_ME_BOX-1,I,K))/                                 &
               (2.0D0*H_DELTA_X))*H_DELTA_X+H_SUM_RX_MU
END DO
!-------------------------- CHANGE IN CONTROL VOLUME MOMENTUM -------------------------------------
DO J= H_NS_BOX,H_NE_BOX
    DO I= H_MS_BOX,H_ME_BOX
        IF (J==H_NS_BOX.OR.J==H_NE_BOX.OR.I==H_MS_BOX.OR.I==H_ME_BOX) THEN
            A= 0.5D0
        ELSE
            A= 1.0D0
        END IF
!
        H_SUM_X_ROU    =-RO*UX(I,J,K)*A*H_DELTA_X*H_DELTA_X+H_SUM_X_ROU
        H_SUM_X_ROU_OLD= RO*UX_OLD(I,J,K)*A*H_DELTA_X*H_DELTA_X+H_SUM_X_ROU_OLD
!
        H_SUM_Y_ROU    =-RO*UY(I,J,K)*A*H_DELTA_X*H_DELTA_X+H_SUM_Y_ROU
        H_SUM_Y_ROU_OLD= RO*UY_OLD(I,J,K)*A*H_DELTA_X*H_DELTA_X+H_SUM_Y_ROU_OLD
!
    END DO
END DO
H_SUM_X_ROU= H_SUM_X_ROU/DT+H_SUM_X_ROU_OLD/DT
H_SUM_Y_ROU= H_SUM_Y_ROU/DT+H_SUM_Y_ROU_OLD/DT
!---------------------------------- CHANGE IN BODY MOMENTUM ---------------------------------------
BODY_MOMENTUM= RO*PI*H_RADIUS*H_RADIUS*(UC-UC_OLD)/DT
!--------------------------------- CALCULATING THE FORCES ---
FX= H_SUM_LX_P+H_SUM_RX_P+H_SUM_BX_UROU+H_SUM_TX_UROU+H_SUM_LX_UROU+H_SUM_RX_UROU+                 &
    H_SUM_BX_MU+H_SUM_TX_MU+H_SUM_LX_MU+H_SUM_RX_MU+H_SUM_X_ROU+BODY_MOMENTUM
!
FY= H_SUM_BY_P+H_SUM_TY_P+H_SUM_BY_UROU+H_SUM_TY_UROU+H_SUM_LY_UROU+H_SUM_RY_UROU+                 &
    H_SUM_BY_MU+H_SUM_TY_MU+H_SUM_LY_MU+H_SUM_RY_MU+H_SUM_Y_ROU
!
CX= FX/(0.5D0*RO*UM*UM)
CY= FY/(0.5D0*RO*UM*UM)
!
END SUBROUTINE
