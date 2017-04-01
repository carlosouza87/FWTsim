PROGRAM FWTsim

IMPLICIT NONE

INTEGER                                   :: Nsys        ! Number of system states
INTEGER                                   :: Nsteps      ! Number of simulation time steps
INTEGER                                   :: k1          ! Count variable for DO construct

REAL                                      :: dt          ! Time step [s]
REAL                                      :: ti tf       ! Initial and final time [s]
REAL, DIMENSION(1,:), ALLOCATABLE         :: Time        ! Time array [s]
REAL                                      :: t           ! Current time during simulation [s]
REAL, DIMENSION(2,1)                      :: eta0        ! Initial positions [m;rad]
REAL, DIMENSION(2,1)                      :: nu0         ! Initial velocities [m/s;rad/s]
REAL, DIMENSION(:,:), ALLOCATABLE         :: X           ! System state [multi-dimensional]
REAL, DIMENSION(2,2) 					  :: Mrb         ! Rigid-body inertia matrix [kg kg.m;kg.m kg.m^2]
REAL, DIMENSION(2,2)                      :: Ainf        ! Infinite added mass matrix [kg kg.m;kg.m kg.m^2]
REAL                                      :: Tret        ! Retardation function time step [s]
REAL, DIMENSION(1,:), ALLOCATABLE         :: K11 K15 K51 K55  ! Retardation functions [N/m N N N.m]
REAL, DIMENSION(2,2)                      :: Dl          ! Linear damping matrix [N.s/m -;- N.s]
REAL, DIMENSION(2,2)                      :: Dq          ! Quadratic damping matrix [N.s^2/m^2 -;- N.s^2]
REAL, DIMENSION(2,2)                      :: Chs         ! Hydrostatic stiffness matrix [N/m N;N N.m]
REAL, DIMENSION(2,2)                      :: Cmr         ! Mooring stiffness matrix [N/m N;N N.m]


Nsteps = (ti-tf+dt)/dt    ! Number of time steps in simulation
Nsys = 4                  ! Only eta and nu make the system

simloop: DO k1=1,Nsteps
    t = Time(k1)
END DO simloop
    	
END PROGRAM FWTsim