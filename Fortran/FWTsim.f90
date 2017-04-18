PROGRAM FWTsim

USE sysvar
USE syssub

IMPLICIT NONE

REAL, DIMENSION(:), ALLOCATABLE            :: Time        ! Time array [s]
REAL, DIMENSION(:,:), ALLOCATABLE          :: X           ! System state [multi-dimensional]


CALL read_inp

test: DO k1 = 1,Nelem
    WRITE(*,*) Blade_dim(k1,:)
END DO test 


Nsys = Ndof*2 ! Order of system
Nsteps = (tf-ti+dt)/dt    ! Number of time steps in simulation

ALLOCATE (X(Nsys,Nsteps), Time(Nsteps))

Time = [(ti+(i-1)*dt, i=1,Nsteps)] ! Time vector


X(1:Ndof,1) = eta0(:)
X(Ndof+1:2*Ndof,1) = nu0(:)

simloop: DO k1=2,Nsteps
    t = Time(k1)
    CALL RK4th(X(1:Nsys,k1-1),t,dt,X(1:Nsys,k1))
END DO simloop



OPEN (UNIT=100,FILE="results.txt",ACTION="write",STATUS="replace")
writeloop: DO k1 = 1, Nsteps
    WRITE (100,*) Time(k1), X(1,k1), X(2,k1), X(3,k1), X(4,k1)
END DO writeloop

CLOSE (100)

DEALLOCATE (X, Time, Blade_dim, Foil_prop)

    	
END PROGRAM FWTsim
