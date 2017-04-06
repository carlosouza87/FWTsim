PROGRAM FWTsim

USE sysvar

IMPLICIT NONE

INTEGER                                   :: Nsys        ! Order of system
INTEGER                                   :: Nsteps      ! Number of simulation time steps
INTEGER                                   :: k1          ! Count variable for DO construct
INTEGER                                   :: i           ! Counter

REAL                                      :: dt          ! Time step [s]
REAL                                      :: ti, tf      ! Initial and final time [s]
REAL, DIMENSION(20001)                    :: Time        ! Time array [s]
REAL                                      :: t           ! Current time during simulation [s]
REAL                                      :: eta0        ! Initial positions [m;rad]
REAL                                      :: nu0         ! Initial velocities [m/s;rad/s]
!REAL, DIMENSION(2,1)                      :: eta0        ! Initial positions [m;rad]
!REAL, DIMENSION(2,1)                      :: nu0         ! Initial velocities [m/s;rad/s]
REAL, DIMENSION(2,20001)                  :: X           ! System state [multi-dimensional]
! REAL, EXTERNAL                            :: sysdyn      ! Output of sysdyn function
!REAL, DIMENSION(2,2) 					  :: Mrb         ! Rigid-body inertia matrix [kg kg.m;kg.m kg.m^2]'
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
Nsteps = (tf-ti+dt)/dt    ! Number of time steps in simulation
!Nsteps = 20001
Time = [(ti+(i-1)*dt, i=1,Nsteps)]

CALL readsysprop
WRITE(*,*) Mrb


eta0 = -20.0
nu0 = 0

Nsys = 2                  ! Only eta and nu make the system

X(1,1) = eta0
X(2,1) = nu0

simloop: DO k1=2,Nsteps
    t = Time(k1)
!	CALL RK4th(X(1:2,k1-1),t,dt,Nsys,sysdyn,X(1:2,k1))
    CALL RK4th(X(1:2,k1-1),t,dt,Nsys,X(1:2,k1))
END DO simloop



OPEN (UNIT=100,FILE="results.txt",ACTION="write",STATUS="replace")
writeloop: DO k1 = 1, Nsteps
    WRITE (100,*) Time(k1), X(1,k1), X(2,k1)
END DO writeloop

CLOSE (100)

CONTAINS

!SUBROUTINE RK4th(x0,t,dt,N,f,x)
SUBROUTINE RK4th(x0,t,dt,N,x)

    IMPLICIT NONE

    INTEGER                                   :: iter        ! Current RK iteration
    INTEGER                                   :: N           ! Order of system

    REAL                                      :: t           ! Time instant
    REAL                                      :: dt          ! Time step
    REAL, DIMENSION(2)                        :: x0          ! Current state
    REAL, DIMENSION(2)                        :: x_p         ! State to integrate
    REAL, DIMENSION(2)                        :: x           ! Integrated state
    REAL, DIMENSION(2)                        :: k1          ! State of first iteration
    REAL, DIMENSION(2)                        :: k2          ! State of second iteration
    REAL, DIMENSION(2)                        :: k3          ! State of third iteration
    REAL, DIMENSION(2)                        :: k4          ! State of fourth iteration
   ! REAL, DIMENSION(:), ALLOCATABLE           :: f           ! Input from function
   ! EXTERNAL                                  :: f           ! Function to be integrated    



    iter = 1
	
!	write(*,*) x0, N, t, dt
!   k1(1:N) = f(x0(1:N),t,dt,iter)
!   k1(1:N) = sysdyn(x0(1:N),t,dt,iter,x_p(1:N))
    CALL sysdyn(x0(1:N),t,dt,iter,k1(1:N))

    iter = 2
!    k2(1:N) = f(x0(1:N)+k1(1:N)*dt/2,t+dt/2,dt,iter)
!    k2(1:N) = sysdyn(x0(1:N)+k1(1:N)*dt/2,t+dt/2,dt,iter,x_p(1:N))
    CALL sysdyn(x0(1:N)+k1(1:N)*dt/2,t+dt/2,dt,iter,k2(1:N))
    iter = 3
!    k3(1:N) = f(x0(1:N)+k2(1:N)*dt/2,t+dt/2,dt,iter)
!    k3(1:N) = sysdyn(x0(1:N)+k2(1:N)*dt/2,t+dt/2,dt,iter,x_p(1:N))
    CALL sysdyn(x0(1:N)+k2(1:N)*dt/2,t+dt/2,dt,iter,k3(1:N))
    iter = 4
!    k4(1:N) = f(x0(1:N)+k3(1:N)*dt,t+dt,dt,iter)
!    k4(1:N) = sysdyn(x0(1:N)+k3(1:N)*dt,t+dt,dt,iter,x_p(1:N))
    CALL sysdyn(x0(1:N)+k3(1:N)*dt,t+dt,dt,iter,k4(1:N))

    x(1:N) = x0(1:N) + (dt/6)*(k1(1:N)+2*k2(1:N)+2*k3(1:N)+k4(1:N))

END SUBROUTINE RK4th

SUBROUTINE sysdyn(x,t,dt,iter,x_p)

IMPLICIT NONE

	INTEGER                                   :: iter        ! Current RK iteration
	INTEGER                                   :: N           ! Order of system

	REAL                                      :: t           ! Time instant
	REAL                                      :: dt          ! Time step
	REAL, DIMENSION(2)                        :: x           ! Current state
	REAL, DIMENSION(2)                        :: x_p         ! State to integrate
	REAL                                      :: eta         ! Position state
	REAL                                      :: nu          ! Velocity state
	REAL                                      :: eta_p       ! Time derivative of position state
	REAL                                      :: nu_p        ! Time derivative of velocity state
	REAL                                      :: Mrb         ! Rigid-body inertia matrix
	REAL                                      :: A0          ! Zero-frequency added mass matrix
	REAL                                      :: Dl          ! Linear viscous damping matrix
	REAL                                      :: Cmr         ! Mooring stiffness matrix

	N = 2

	Mrb = 7794048
	A0 = 7.99e6
	Dl = 1e5
	Cmr = 41180

	eta = x(1)
	nu = x(2)

	eta_p = nu
	nu_p = (-Dl*nu-Cmr*eta)/(Mrb+A0)

	x_p(1) = eta_p
	x_p(N/2+1:N) = nu_p

END SUBROUTINE sysdyn
    	
END PROGRAM FWTsim
