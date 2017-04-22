program fwtsim

use sysvar
use syssub

implicit none

real(8), dimension(:), allocatable         :: Time        ! Time array [s]
real(8), dimension(:,:), allocatable       :: X           ! System state [multi-dimensional]


call read_inp

! test: do k1 = 1,Nelem
    ! write(*,*) Blade_dim(k1,:)
! end do test 


Nsys = Ndof*2 ! Order of system
Nsteps = (tf-ti+dt)/dt    ! Number of time steps in simulation

allocate (X(Nsys,Nsteps), Time(Nsteps))

Time = [(ti+(i-1)*dt, i=1,Nsteps)] ! Time vector


X(1:Ndof,1) = eta0(:)
X(Ndof+1:2*Ndof,1) = nu0(:)

simloop: do k1=2,Nsteps
    t = Time(k1)
    call RK4th(X(1:Nsys,k1-1),t,dt,X(1:Nsys,k1))
end do simloop

open (unit=100,file="results.txt",action="write",status="replace")
writeloop: do k1 = 1, Nsteps
    write (100,*) time(k1), x(1,k1), x(2,k1), x(3,k1), x(4,k1)
end do writeloop

close (100)

deallocate (x, time, blade_dim, foil_prop)

    	
end program FWTsim
