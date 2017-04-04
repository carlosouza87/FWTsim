FUNCTION sysdyn(x,t,dt,iter)

IMPLICIT NONE

INTEGER                                   :: iter        ! Current RK iteration

REAL                                      :: t           ! Time instant
REAL                                      :: dt          ! Time step
REAL, DIMENSION(:), ALLOCATABLE           :: x           ! Current state
REAL, DIMENSION(:), ALLOCATABLE           :: sysdyn      ! Time derivative of x
REAL, DIMENSION(:), ALLOCATABLE           :: eta         ! Position state
REAL, DIMENSION(:), ALLOCATABLE           :: nu          ! Velocity state
REAL, DIMENSION(:), ALLOCATABLE           :: eta_p        ! Time derivative of position state
REAL, DIMENSION(:), ALLOCATABLE           :: nu_p         ! Time derivative of velocity state
REAL, DIMENSION(1,1)                      :: Mrb         ! Rigid-body inertia matrix
REAL, DIMENSION(1,1)                      :: A0          ! Zero-frequency added mass matrix
REAL, DIMENSION(1,1)                      :: Dl          ! Linear viscous damping matrix
REAL, DIMENSION(1,1)                      :: Cmr         ! Mooring stiffness matrix


Mrb = 7794048
A0 = 7.99e6
Dl = 1e5
Cmr = 41180

eta = x(1)
nu = x(2)

eta_p = nu;
nu_p = (-Dl*nu-Cmr*eta)/(Mrb+A0)

sysdyn = [eta_p;nu_p]