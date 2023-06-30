!------ The auxillary points in the interpolation stencil are identified in this subroutine -------
!--------------------------------------------------------------------------------------------------
SUBROUTINE H_INTERPOLATION_DATA(H_NX,H_N_NB,H_X_NB,H_Y_NB,H_SLOPE_NB,H_B_NB,DIS_1,H_LAN_X,H_LAN_Y,      &
                           H_DIS_INPOL,I_IP,J_IP,H_I_NB,H_J_NB,H_X_INT,H_Y_INT)
USE IB_VARIABLES, ONLY: H_DELTA_X,H_DEL
!
IMPLICIT NONE
INTEGER                          :: I
INTEGER                          :: H_NX
REAL(8)                          :: PI
REAL(8)                          :: DELTA_1,DELTA_2,DELTA_3,DELTA_4
REAL(8)                          :: MIN_1,MIN_2,MIN_3,MIN_4
INTEGER                          :: H_N_NB
REAL(8),DIMENSION(1:4)           :: X_IP,Y_IP
REAL(8),DIMENSION(1:H_NX)        :: DIS_1,H_LAN_X,H_LAN_Y,H_X_NB,H_Y_NB,H_B_NB,                    &
                                    H_SLOPE_NB
INTEGER,DIMENSION(1:H_NX)        :: H_I_NB,H_J_NB
REAL(8),DIMENSION(1:3,1:H_NX)    :: H_DIS_INPOL
REAL(8),DIMENSION(1:2,1:H_NX)    :: H_X_INT,H_Y_INT
INTEGER,DIMENSION(1:4,1:2,1:H_NX):: I_IP,J_IP
!
PI= ACOS(-1.0D0)
DO I=1,H_N_NB
!
   H_DIS_INPOL(1,I)= DIS_1(I)
!--------------------------------------------------------------------------------------------------
   X_IP(1)= H_X_NB(I)+DBLE(NINT(DSIGN(1.0D0,H_LAN_X(I))))*H_DELTA_X
   Y_IP(1)= H_SLOPE_NB(I)*X_IP(1)+H_B_NB(I)
   DELTA_1= DSQRT((X_IP(1)-H_X_NB(I))*(X_IP(1)-H_X_NB(I))+                                         &
                  (Y_IP(1)-H_Y_NB(I))*(Y_IP(1)-H_Y_NB(I)))
!
   Y_IP(2)= H_Y_NB(I)+DBLE(NINT(DSIGN(1.0D0,H_LAN_Y(I))))*H_DEL
   X_IP(2)=(Y_IP(2)-H_B_NB(I))/H_SLOPE_NB(I)
   DELTA_2= DSQRT((X_IP(2)-H_X_NB(I))*(X_IP(2)-H_X_NB(I))+                                         &
                  (Y_IP(2)-H_Y_NB(I))*(Y_IP(2)-H_Y_NB(I)))
!
   X_IP(3)= H_X_NB(I)+DBLE(NINT(DSIGN(2.0D0,H_LAN_X(I))))*H_DELTA_X
   Y_IP(3)= H_SLOPE_NB(I)*X_IP(3)+H_B_NB(I)
   DELTA_3= DSQRT((X_IP(3)-H_X_NB(I))*(X_IP(3)-H_X_NB(I))+                                         &
                  (Y_IP(3)-H_Y_NB(I))*(Y_IP(3)-H_Y_NB(I)))
!
   Y_IP(4)= H_Y_NB(I)+DBLE(NINT(DSIGN(2.0D0,H_LAN_Y(I))))*H_DEL
   X_IP(4)=(Y_IP(4)-H_B_NB(I))/H_SLOPE_NB(I)
   DELTA_4= DSQRT((X_IP(4)-H_X_NB(I))*(X_IP(4)-H_X_NB(I))+                                         &
                  (Y_IP(4)-H_Y_NB(I))*(Y_IP(4)-H_Y_NB(I)))
!
   MIN_1  = MIN(DELTA_1,DELTA_2)
   MIN_2  = MIN(DELTA_2,DELTA_3)
   MIN_3  = MIN(DELTA_1,DELTA_4)
   MIN_4  = MIN(DELTA_3,DELTA_4)
!
   IF (DABS(MIN_1-DELTA_1)<=10.E-10.AND.DABS(MIN_2-DELTA_3)<=10.E-10) THEN
       I_IP(1,1,I) = H_I_NB(I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(1,1,I) = H_J_NB(I)
       I_IP(3,1,I) = H_I_NB(I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(3,1,I) = H_J_NB(I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
!
       I_IP(2,1,I) = I_IP(1,1,I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(2,1,I) = J_IP(1,1,I)
       I_IP(4,1,I) = I_IP(3,1,I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(4,1,I) = J_IP(3,1,I)
       H_X_INT(1,I)= X_IP(1)
       H_Y_INT(1,I)= Y_IP(1)
!
       I_IP(1,2,I) = H_I_NB(I)+NINT(DSIGN(2.0D0,H_LAN_X(I)))
       J_IP(1,2,I) = H_J_NB(I)
       I_IP(3,2,I) = H_I_NB(I)+NINT(DSIGN(2.0D0,H_LAN_X(I)))
       J_IP(3,2,I) = H_J_NB(I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
!
       I_IP(2,2,I) = I_IP(1,2,I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(2,2,I) = J_IP(1,2,I)
       I_IP(4,2,I) = I_IP(3,2,I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(4,2,I) = J_IP(3,2,I)
       H_X_INT(2,I)= X_IP(3)
       H_Y_INT(2,I)= Y_IP(3)
!
   ELSEIF (DABS(MIN_1-DELTA_1)<=10.E-10.AND.DABS(MIN_2-DELTA_2)<=10.E-10)THEN
       I_IP(1,1,I) = H_I_NB(I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(1,1,I) = H_J_NB(I)
       I_IP(3,1,I) = H_I_NB(I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(3,1,I) = H_J_NB(I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
!
       I_IP(2,1,I) = I_IP(1,1,I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(2,1,I) = J_IP(1,1,I)
       I_IP(4,1,I) = I_IP(3,1,I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(4,1,I) = J_IP(3,1,I)
       H_X_INT(1,I)= X_IP(1)
       H_Y_INT(1,I)= Y_IP(1)
!
       I_IP(1,2,I) = H_I_NB(I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(1,2,I) = H_J_NB(I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       I_IP(2,2,I) = H_I_NB(I)+NINT(DSIGN(2.0D0,H_LAN_X(I)))
       J_IP(2,2,I) = H_J_NB(I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
!
       I_IP(3,2,I) = I_IP(1,2,I)
       J_IP(3,2,I) = J_IP(1,2,I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       I_IP(4,2,I) = I_IP(2,2,I)
       J_IP(4,2,I) = J_IP(2,2,I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       H_X_INT(2,I)= X_IP(2)
       H_Y_INT(2,I)= Y_IP(2)
!
   ELSEIF (DABS(MIN_1-DELTA_2)<=10.E-10.AND.DABS(MIN_3-DELTA_1)<=10.E-10) THEN
       I_IP(1,1,I) = H_I_NB(I)
       J_IP(1,1,I) = H_J_NB(I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       I_IP(2,1,I) = H_I_NB(I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(2,1,I) = H_J_NB(I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
!
       I_IP(3,1,I) = I_IP(1,1,I)
       J_IP(3,1,I) = J_IP(1,1,I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       I_IP(4,1,I) = I_IP(2,1,I)
       J_IP(4,1,I) = J_IP(2,1,I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       H_X_INT(1,I)= X_IP(2)
       H_Y_INT(1,I)= Y_IP(2)
!
       I_IP(1,2,I) = H_I_NB(I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(1,2,I) = H_J_NB(I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       I_IP(3,2,I) = H_I_NB(I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(3,2,I) = H_J_NB(I)+NINT(DSIGN(2.0D0,H_LAN_Y(I)))
!
       I_IP(2,2,I) = I_IP(1,2,I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(2,2,I) = J_IP(1,2,I)
       I_IP(4,2,I) = I_IP(3,2,I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(4,2,I) = J_IP(3,2,I)
       H_X_INT(2,I)= X_IP(1)
       H_Y_INT(2,I)= Y_IP(1)
!
   ELSEIF (H_SLOPE_NB(I)==TAN(PI/2.)) THEN
       I_IP(1,1,I) = H_I_NB(I)
       J_IP(1,1,I) = H_J_NB(I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       I_IP(2,1,I) = H_I_NB(I)
       J_IP(2,1,I) = H_J_NB(I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
!
       I_IP(3,1,I) = I_IP(1,1,I)
       J_IP(3,1,I) = J_IP(1,1,I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       I_IP(4,1,I) = I_IP(2,1,I)
       J_IP(4,1,I) = J_IP(2,1,I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       H_X_INT(1,I)= H_X_NB(I)
       H_Y_INT(1,I)= H_Y_NB(I)+DBLE(NINT(DSIGN(1.0D0,H_LAN_Y(I))))*H_DELTA_X
        !
       I_IP(1,2,I) = H_I_NB(I)
       J_IP(1,2,I) = H_J_NB(I)+NINT(DSIGN(2.0D0,H_LAN_Y(I)))
       I_IP(2,2,I) = H_I_NB(I)
       J_IP(2,2,I) = H_J_NB(I)+NINT(DSIGN(2.0D0,H_LAN_Y(I)))
!
       I_IP(3,2,I) = I_IP(1,2,I)
       J_IP(3,2,I) = J_IP(1,2,I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       I_IP(4,2,I) = I_IP(2,2,I)
       J_IP(4,2,I) = J_IP(2,2,I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       H_X_INT(2,I)= H_X_NB(I)
       H_Y_INT(2,I)= H_Y_NB(I)+DBLE(NINT(DSIGN(2.0D0,H_LAN_Y(I))))*H_DELTA_X
!
   ELSEIF (DABS(MIN_1-DELTA_2)<=10.E-10.AND. DABS(MIN_3-DELTA_4)<=10.E-10) THEN
       I_IP(1,1,I) = H_I_NB(I)
       J_IP(1,1,I) = H_J_NB(I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       I_IP(2,1,I) = H_I_NB(I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(2,1,I) = H_J_NB(I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))

       I_IP(3,1,I) = I_IP(1,1,I)
       J_IP(3,1,I) = J_IP(1,1,I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       I_IP(4,1,I) = I_IP(2,1,I)
       J_IP(4,1,I) = J_IP(2,1,I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       H_X_INT(1,I)= X_IP(2)
       H_Y_INT(1,I)= Y_IP(2)
!
       I_IP(1,2,I) = H_I_NB(I)
       J_IP(1,2,I) = H_J_NB(I)+NINT(DSIGN(2.0D0,H_LAN_Y(I)))
       I_IP(2,2,I) = H_I_NB(I)+NINT(DSIGN(1.0D0,H_LAN_X(I)))
       J_IP(2,2,I) = H_J_NB(I)+NINT(DSIGN(2.0D0,H_LAN_Y(I)))
!
       I_IP(3,2,I) = I_IP(1,2,I)
       J_IP(3,2,I) = J_IP(1,2,I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       I_IP(4,2,I) = I_IP(2,2,I)
       J_IP(4,2,I) = J_IP(2,2,I)+NINT(DSIGN(1.0D0,H_LAN_Y(I)))
       H_X_INT(2,I)= X_IP(4)
       H_Y_INT(2,I)= Y_IP(4)
   END IF
!
   H_DIS_INPOL(2,I)= DSQRT((H_X_INT(1,I)-H_X_NB(I))*(H_X_INT(1,I)-H_X_NB(I))+                      &
                           (H_Y_INT(1,I)-H_Y_NB(I))*(H_Y_INT(1,I)-H_Y_NB(I)))
   H_DIS_INPOL(3,I)= DSQRT((H_X_INT(2,I)-H_X_NB(I))*(H_X_INT(2,I)-H_X_NB(I))+                      &
                           (H_Y_INT(2,I)-H_Y_NB(I))*(H_Y_INT(2,I)-H_Y_NB(I)))
!--------------------------------------------------------------------------------------------------
END DO
END SUBROUTINE
