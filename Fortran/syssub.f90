MODULE syssub

SUBROUTINE read_inp
! Reads input files (system properties and simulation parameters)

OPEN (UNIT=110, FILE='inp_sysprop.txt', STATUS='old', ACTION='read')

    ! Read number of DOF
    READ(110,*) Ndof

    ! Read Mrb
    READ(110,*) Mrb(1,:)
    READ(110,*) Mrb(2,:)
	
	! Read Add
	READ(110,*) Add(1,:)
    READ(110,*) Add(2,:)
	
    ! Read Dl
	READ(110,*) Dl(1,:)
    READ(110,*) Dl(2,:)
	
    ! Read Dq
	READ(110,*) Dq(1,:)
    READ(110,*) Dq(2,:)
	
	! Read Chs
	READ(110,*) Chs(1,:)
    READ(110,*) Chs(2,:)
	
	! Read Cmr
	READ(110,*) Cmr(1,:)
    READ(110,*) Cmr(2,:)

CLOSE (110)

OPEN (UNIT=120, FILE='inp_simpar.txt', STATUS='old', ACTION='read')

	! ***********************************************************
	! FALTA LER ARQUIVO COM PARÂMETROS DA SIMULACAO!!!
	! ***********************************************************

CLOSE (120)

END SUBROUTINE read_inp

SUBROUTINE RK4th(x0,t,dt,N,x)

    IMPLICIT NONE

    INTEGER                                   :: iter        ! Current RK iteration
    INTEGER                                   :: N           ! Order of system

    REAL                                      :: t           ! Time instant
    REAL                                      :: dt          ! Time step
    REAL, DIMENSION(:)                        :: x0          ! Current state 
    REAL, DIMENSION(:)                        :: x           ! Integrated state
    REAL, DIMENSION(:), ALLOCATABLE           :: k1          ! State of first iteration
    REAL, DIMENSION(:), ALLOCATABLE           :: k2          ! State of second iteration
    REAL, DIMENSION(:), ALLOCATABLE           :: k3          ! State of third iteration
    REAL, DIMENSION(:), ALLOCATABLE           :: k4          ! State of fourth iteration

	ALLOCATE (k1(N), k2(n), k3(N), k4(N))
	
    iter = 1	
    CALL sysdyn(x0(1:N),t,dt,iter,k1(1:N))
	
    iter = 2
    CALL sysdyn(x0(1:N)+k1(1:N)*dt/2,t+dt/2,dt,iter,k2(1:N))
	
    iter = 3 
    CALL sysdyn(x0(1:N)+k2(1:N)*dt/2,t+dt/2,dt,iter,k3(1:N))
	
    iter = 4
    CALL sysdyn(x0(1:N)+k3(1:N)*dt,t+dt,dt,iter,k4(1:N))

    x(1:N) = x0(1:N) + (dt/6)*(k1(1:N)+2*k2(1:N)+2*k3(1:N)+k4(1:N))
	
	DEALLOCATE (k1, k2, k3, k4)

END SUBROUTINE RK4th

SUBROUTINE sysdyn(x,t,dt,iter,x_p)

IMPLICIT NONE

	INTEGER                                   :: iter        ! Current RK iteration
	INTEGER                                   :: N           ! Order of system

	REAL, DIMENSION(:)                        :: x           ! Current state
	REAL, DIMENSION(:)                        :: x_p         ! State to integrate
	REAL, DIMENSION(:), ALLOCATABLE           :: eta         ! Position state
	REAL, DIMENSION(:), ALLOCATABLE           :: nu          ! Velocity state
	REAL, DIMENSION(:), ALLOCATABLE           :: eta_p       ! Time derivative of position state
	REAL, DIMENSION(:), ALLOCATABLE           :: nu_p        ! Time derivative of velocity state

	ALLOCATE (eta(N/2), nu(N/2), eta_p(N/2), nu_p(N/2))
	
	eta = x(1:N/2)
	nu = x(N/2+1:N)

	eta_p = nu
	! ***********************************************************
	! FALTA LIDAR COM A MULTIPLICACÃO!!!
	! ***********************************************************
	nu_p = (-Dl*nu-Dq*ABS(nu)*nu-Cmr*eta-Chs*eta)/(Mrb+A0)

	x_p(1) = eta_p
	x_p(N/2+1:N) = nu_p
	
	DEALLOCATE (eta, nu, eta_p, nu_p)

END SUBROUTINE sysdyn

END MODULE syssub