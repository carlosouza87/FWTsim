PROGRAM FWTsim

IMPLICIT NONE

INTEGER                                   :: Nsys        ! Order of system
INTEGER                                   :: Nsteps      ! Number of simulation time steps
INTEGER                                   :: k1          ! Count variable for DO construct
INTEGER                                   :: i           ! Counter

REAL                                      :: dt          ! Time step [s]
REAL                                      :: ti, tf      ! Initial and final time [s]
REAL, DIMENSION(:),   ALLOCATABLE         :: Time        ! Time array [s]
REAL                                      :: t           ! Current time during simulation [s]
REAL, DIMENSION(2,1)                      :: eta0        ! Initial positions [m;rad]
REAL, DIMENSION(2,1)                      :: nu0         ! Initial velocities [m/s;rad/s]
REAL, DIMENSION(:,:), ALLOCATABLE         :: X           ! System state [multi-dimensional]
!REAL, DIMENSION(2,2) 					  :: Mrb         ! Rigid-body inertia matrix [kg kg.m;kg.m kg.m^2]
!REAL, DIMENSION(2,2)                      :: Ainf        ! Infinite added mass matrix [kg kg.m;kg.m kg.m^2]
!REAL                                      :: Tret        ! Retardation function time step [s]
!REAL, DIMENSION(1,:), ALLOCATABLE         :: K11 K15 K51 K55  ! Retardation functions [N/m N N N.m]
!REAL, DIMENSION(2,2)                      :: Dl          ! Linear damping matrix [N.s/m -;- N.s]
!REAL, DIMENSION(2,2)                      :: Dq          ! Quadratic damping matrix [N.s^2/m^2 -;- N.s^2]
!REAL, DIMENSION(2,2)                      :: Chs         ! Hydrostatic stiffness matrix [N/m N;N N.m]
!REAL, DIMENSION(2,2)                      :: Cmr         ! Mooring stiffness matrix [N/m N;N N.m]

!EXTERNAL SYSDYN

ti = 0
tf = 2000
dt = 0.1
Nsteps = (ti-tf+dt)/dt    ! Number of time steps in simulation
Time = [(ti*(i-1)*dt, i=1,Nsteps)]

eta0 = -20.0
nu0 = 0

Nsys = 2                  ! Only eta and nu make the system

X(1,1) = eta0
X(2,1) = nu0

simloop: DO k1=2,Nsteps
    t = Time(k1)
	call RK4th(X(1:2,k1-1),t,dt,Nsys,sysdyn,X(1:2,k1))
END DO simloop

CONTAINS

SUBROUTINE RK4th(x0,t,dt,N,f,x)

    IMPLICIT NONE

    INTEGER                                   :: iter        ! Current RK iteration
    INTEGER                                   :: N           ! Order of system

    REAL                                      :: t           ! Time instant
    REAL                                      :: dt          ! Time step
    REAL, DIMENSION(:), ALLOCATABLE           :: x0          ! Current state
    REAL, DIMENSION(:), ALLOCATABLE           :: xp          ! Integrated state
    REAL, DIMENSION(:), ALLOCATABLE           :: x           ! Integrated state
    REAL, DIMENSION(:), ALLOCATABLE           :: k1          ! State of first iteration
    REAL, DIMENSION(:), ALLOCATABLE           :: k2          ! State of second iteration
    REAL, DIMENSION(:), ALLOCATABLE           :: k3          ! State of third iteration
    REAL, DIMENSION(:), ALLOCATABLE           :: k4          ! State of fourth iteration

    EXTERNAL                                  :: f           ! Function to be integrated             

    iter = 1
    k1(1:N) = f(x0(1:N),t,dt,iter)

    iter = 2
    k2(1:N) = f(x0(1:N)+k1(1:N)*dt/2,t+dt/2,dt,iter)

    iter = 3
    k3(1:N) = f(x0(1:N)+k2(1:N)*dt/2,t+dt/2,dt,iter)

    iter = 4
    k4(1:N) = f(x0(1:N)+k3(1:N)*dt,t+dt,dt,iter)

    x(1:N) = x0(1:N) + (dt/6)*(k1(1:N)+2*k2(1:N)+2*k3(1:N)+k4(1:N))

END SUBROUTINE RK4th
    	
END PROGRAM FWTsim

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

sysdyn(1:N/2) = eta_p
sysdyn(N/2+1:N) = nu_p

END FUNCTION sysdyn