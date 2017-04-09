MODULE syssub
USE sysvar

IMPLICIT NONE

CONTAINS

SUBROUTINE read_inp
! Reads input files (system properties and simulation parameters)

OPEN (UNIT=200, FILE='inp_sysprop.txt', STATUS='old', ACTION='read')

    ! Read number of DOF
    READ(200,*) Ndof
	
	testNdof: IF (Ndof > 2) THEN
                WRITE(*,*) 'Sorry, max. 2 DOF so far. '
        STOP
    END IF testNdof
	
	! Read Mrb
	readMrb: DO k1 = 1,Ndof
	    READ(200,*) Mrb(k1,:)
	END DO readMrb		
	
	! Read Add
	readAdd: DO k1 = 1,Ndof
	    READ(200,*) Add(k1,:)
	END DO readAdd

	! Read Dl
	readDl: DO k1 = 1,Ndof
	    READ(200,*) Dl(k1,:)
	END DO readDl
	
   ! Read Dq
	readDq: DO k1 = 1,Ndof
	    READ(200,*) Dq(k1,:)
	END DO readDq
	
	! Read Chs
	readChs: DO k1 = 1,Ndof
	    READ(200,*) Chs(k1,:)
	END DO readChs
	
	! Read Cmr
	readCmr: DO k1 = 1,Ndof
	    READ(200,*) Cmr(k1,:)
	END DO readCmr	

    ! ! Read Mrb
    ! READ(200,*) Mrb(1,:)
    ! READ(200,*) Mrb(2,:)
	
	! ! Read Add
	! READ(200,*) Add(1,:)
    ! READ(200,*) Add(2,:)
	
    ! ! Read Dl
	! READ(200,*) Dl(1,:)
    ! READ(200,*) Dl(2,:)
	
    ! ! Read Dq
	! READ(200,*) Dq(1,:)
    ! READ(200,*) Dq(2,:)
	
	! ! Read Chs
	! READ(200,*) Chs(1,:)
    ! READ(200,*) Chs(2,:)
	
	! ! Read Cmr
	! READ(200,*) Cmr(1,:)
    ! READ(200,*) Cmr(2,:)

CLOSE (200)

OPEN (UNIT=120, FILE='inp_simpar.txt', STATUS='old', ACTION='read')

    ! Read initial time
	READ(120,*) ti
	
	! Read final time
	READ(120,*) tf
		
	! Read time step
	READ(120,*) dt	
	
	! Read eta0
	readeta0: DO k1 = 1,Ndof
	    READ(120,*) eta0(k1)
	END DO readeta0
	
	! Read nu0
	readnu0: DO k1 = 1,Ndof
	    READ(120,*) nu0(k1)
	END DO readnu0

CLOSE (120)

END SUBROUTINE read_inp

SUBROUTINE RK4th(x0,t,dt,x)

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

	N = SIZE(x0,1)
	
	ALLOCATE (k1(N), k2(n), k3(N), k4(N))
	
    iter = 1	
    CALL sysdyn(x0(1:N),t,dt,iter,k1(1:N))
	
    iter = 2
    CALL sysdyn(x0(1:N)+k1(1:N)*dt/2,t+dt/2,dt,iter,k2(1:N))
	
    iter = 3 
    CALL sysdyn(x0(1:N)+k2(1:N)*dt/2,t+dt/2,dt,iter,k3(1:N))
	
    iter = 4
    CALL sysdyn(x0(1:N)+k3(1:N)*dt,t+dt,dt,iter,k4(1:N))
	
	! write (*,*) 'Fora Temer!'

    x(1:N) = x0(1:N) + (dt/6)*(k1(1:N)+2*k2(1:N)+2*k3(1:N)+k4(1:N))
	
	DEALLOCATE (k1, k2, k3, k4)

END SUBROUTINE RK4th

SUBROUTINE sysdyn(x,t,dt,iter,x_p)

IMPLICIT NONE

	INTEGER                                   :: iter        ! Current RK iteration
	INTEGER                                   :: N           ! Order of system

	REAL                                      :: t           ! Time instant
    REAL                                      :: dt          ! Time step
	REAL, DIMENSION(4)                        :: x           ! Current state
	REAL, DIMENSION(4)                        :: x_p         ! State to integrate
	REAL, DIMENSION(2,1)                      :: eta         ! Position state
	REAL, DIMENSION(2,1)                      :: nu          ! Velocity state
	REAL, DIMENSION(2,1)                      :: eta_p       ! Time derivative of position state
	REAL, DIMENSION(2,1)                      :: nu_p        ! Time derivative of velocity state
	REAL, DIMENSION(2,1)                      :: Brhs        ! Right-hand-side matrix of eq. of motions
	REAL, DIMENSION(2,2)                      :: Minv        ! Inverse of Mrb+Add
	! REAL, DIMENSION(2,2)                      :: J           ! Transformation matrix (body-fixed to inertial system)
	! REAL, DIMENSION(:)                        :: x           ! Current state
	! REAL, DIMENSION(:)                        :: x_p         ! State to integrate
	! REAL, DIMENSION(:), ALLOCATABLE           :: eta         ! Position state
	! REAL, DIMENSION(:), ALLOCATABLE           :: nu          ! Velocity state
	! REAL, DIMENSION(:), ALLOCATABLE           :: eta_p       ! Time derivative of position state
	! REAL, DIMENSION(:), ALLOCATABLE           :: nu_p        ! Time derivative of velocity state
	! REAL, DIMENSION(:,:), ALLOCATABLE         :: Brhs        ! Right-hand-side matrix of eq. of motions
	! REAL, DIMENSION(:,:), ALLOCATABLE         :: Minv        ! Inverse of Mrb+Add

	! ALLOCATE (eta(N/2), nu(N/2), eta_p(N/2), nu_p(N/2),Brhs(N/2,N/2),Minv(N/2,N/2))
	
	! Initialize position and velocity vectors from input state
	
	N = SIZE(x,1)
	
	eta(:,1) = x(1:N/2)
	nu(:,1) = x(N/2+1:N)
	
	! Calculate derivative of eta, eta_p = nu
	eta_p = nu
	
	! Inversion of Mrb + Add
	Minv = Mrb + Add
	CALL matinv2(Minv)
	
	! Calculate right-hand-side of equations of motions
	Brhs = (-MATMUL(Dl,nu)-MATMUL(Dq,ABS(nu)*nu)-MATMUL(Cmr,eta)-MATMUL(Chs,eta))

	! Calculate derivative of nu, nu_p = Minv*Brhs
	nu_p = MATMUL(Minv,Brhs)

	x_p(1:N/2) = eta_p(:,1)
	x_p(N/2+1:N) = nu_p(:,1)
	
	! DEALLOCATE (eta, nu, eta_p, nu_p)

END SUBROUTINE sysdyn

SUBROUTINE matinv2(A) 
    ! Performs a direct calculation of the inverse of a 2×2 matrix. 
	! Source: http://fortranwiki.org/fortran/show/Matrix+inversion
	
    REAL                    :: A(2,2)   !! Matrix
    REAL                    :: B(2,2)   !! Temp matrix
    REAL                    :: detinv   !! Determinant

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
	
	A(:,:) = B(:,:)
END SUBROUTINE

END MODULE syssub