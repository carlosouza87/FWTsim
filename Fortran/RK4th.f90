FUNCTION RK4th(xk,t,dt,sysdyn)

REAL, DIMENSION(:,1), ALLOCATABLE         :: RK4th          ! Integrated state
REAL, DIMENSION(:,1), ALLOCATABLE         :: xk             ! State to be integrated
REAL                                      :: dt             ! Time step [s]

INTEGER                                   :: iter           ! Iteration number              

EXTERNAL                                  :: sysdyn         ! Function that provides the derivative of the state based on the system dynamics
