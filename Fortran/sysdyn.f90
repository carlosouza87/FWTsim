FUNCTION sysdyn(x,t,dt,iter)

IMPLICIT NONE

INTEGER                                   :: iter        ! Current RK iteration
INTEGER                                   :: N           ! Order of system

REAL                                      :: t           ! Time instant
REAL                                      :: dt          ! Time step
REAL, DIMENSION(2)                        :: x           ! Current state
REAL, DIMENSION(2)                        :: sysdyn      ! Time derivative of x
REAL                                      :: eta         ! Position state
REAL                                      :: nu          ! Velocity state
!REAL, DIMENSION(:), ALLOCATABLE           :: eta         ! Position state
!REAL, DIMENSION(:), ALLOCATABLE           :: nu          ! Velocity state
REAL                                      :: eta_p       ! Time derivative of position state
REAL                                      :: nu_p        ! Time derivative of velocity state
!REAL, DIMENSION(:), ALLOCATABLE           :: eta_p       ! Time derivative of position state
!REAL, DIMENSION(:), ALLOCATABLE           :: nu_p        ! Time derivative of velocity state
REAL                                      :: Mrb         ! Rigid-body inertia matrix
REAL                                      :: A0          ! Zero-frequency added mass matrix
REAL                                      :: Dl          ! Linear viscous damping matrix
REAL                                      :: Cmr         ! Mooring stiffness matrix
!REAL, DIMENSION(1,1)                      :: Mrb         ! Rigid-body inertia matrix
!REAL, DIMENSION(1,1)                      :: A0          ! Zero-frequency added mass matrix
!REAL, DIMENSION(1,1)                      :: Dl          ! Linear viscous damping matrix
!REAL, DIMENSION(1,1)                      :: Cmr         ! Mooring stiffness matrix

N = 2

Mrb = 7794048
A0 = 7.99e6
Dl = 1e5
Cmr = 41180

eta = x(1)
nu = x(2)

eta_p = nu
nu_p = (-Dl*nu-Cmr*eta)/(Mrb+A0)

sysdyn(1) = eta_p
sysdyn(N/2+1:N) = nu_p

END FUNCTION sysdyn