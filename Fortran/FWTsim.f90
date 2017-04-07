PROGRAM FWTsim

USE sysvar, syssub

IMPLICIT NONE

! INTEGER                                   :: Nsys        ! Order of system
! INTEGER                                   :: Nsteps      ! Number of simulation time steps
! INTEGER                                   :: k1          ! Count variable for DO construct
! INTEGER                                   :: i           ! Counter

! REAL                                      :: dt          ! Time step [s]
! REAL                                      :: ti, tf      ! Initial and final time [s]
! REAL, DIMENSION(20001)                    :: Time        ! Time array [s]
! REAL                                      :: t           ! Current time during simulation [s]
! REAL                                      :: eta0        ! Initial positions [m;rad]
! REAL                                      :: nu0         ! Initial velocities [m/s;rad/s]
REAL, DIMENSION(2,:)                  :: X           ! System state [multi-dimensional]
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
WRITE(*,*) 'Mrb'
WRITE(*,*) Mrb(1,:)
WRITE(*,*) Mrb(2,:)

WRITE(*,*) 'Cmr'
WRITE(*,*) Cmr(1,:)
WRITE(*,*) Cmr(2,:)


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

    	
END PROGRAM FWTsim
