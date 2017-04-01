PROGRAM FWTsim

IMPLICIT                        NONE

REAL :: dt ! Time step [s]
REAL :: Time(16001)  ! Time array [s]
REAL :: Mrb(2,2) ! Rigid-body inertia matrix
REAL :: Ainf(2,2) ! Infinite added mass matrix
REAL :: Tret ! Retardation function time step
REAL :: K11(1,200) K15(1,200) K51(1,200) K55(1,200)  ! Retardation functions
REAL :: Dl(2,2) ! Linear damping matrix
REAL :: Dq(2,2) ! Quadratic damping matrix
REAL :: Chs(2,2) ! Hydrostatic stiffness matrix
REAL :: Cmr(2,2) ! Mooring stiffness matrix

END PROGRAM FWTsim