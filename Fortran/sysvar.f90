MODULE sysvar

IMPLICIT NONE

REAL, DIMENSION(2,2)                      :: Mrb         ! Rigid-body inertia matrix
REAL, DIMENSION(2,2)                      :: A0          ! Zero-frequency added mass matrix
REAL, DIMENSION(2,2)                      :: Dl          ! Linear viscous damping matrix
REAL, DIMENSION(2,2)                      :: Chs         ! Hydrostatic stiffness matrix
REAL, DIMENSION(2,2)                      :: Cmr         ! Mooring stiffness matrix

CONTAINS 

SUBROUTINE readsysprop

OPEN (UNIT=110, FILE='inp_sysprop.txt', STATUS='old', ACTION='read')

! Read Mrb
READ(110,*) Mrb(1,:)
READ(110,*) Mrb(2,:)

CLOSE (110)

END SUBROUTINE readsysprop


END MODULE sysvar