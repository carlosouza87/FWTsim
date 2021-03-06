program fwtsim

use sysvar
use syssub

implicit none

real(8), dimension(:), allocatable         :: Time        ! Time array [s]
real(8), dimension(:,:), allocatable       :: X           ! System state [multi-dimensional]
real(8)                                    :: simtime     ! Simulation time

integer,dimension(8) :: time_values

write(*,*) 'Running...'
call date_and_time(VALUES=time_values)
simtime = real(time_values(5))*3600 + real(time_values(6))*60 + real(time_values(7)) &
    + real(time_values(8))/1000
! write(*,*) time_values(5), time_values(6), time_values(7), time_values(8)

call read_inp

Nsys = Ndof*2 + 1 ! Order of system
Nsteps = (tf-ti+dt)/dt    ! Number of time steps in simulation

allocate (X(Nsys,Nsteps), Time(Nsteps), Th_hist(Nsteps), &
    Qrt_hist(Nsteps), Qgen_hist(Nsteps), beta_hist(Nsteps))

Time = [(ti+(i-1)*dt, i=1,Nsteps)] ! Time vector


X(1:Ndof,1) = eta0(:)
X(Ndof+1:2*Ndof,1) = nu0(:)
X(2*Ndof+1,1) = Omg_rt

simloop: do k_time=2,Nsteps
    t = Time(k_time)
    call RK4th(X(1:Nsys,k_time-1),t,dt,X(1:Nsys,k_time))
    ! check_progress: if (mod(t,100.) < 1e-6) then
        ! write(*,*) 't = ', t
    ! end if check_progress
end do simloop

open (unit=100,file="results.txt",action="write",status="replace")
!write(100,*) 'Time ', 'x1 ', 'x2 ', 'x1p ', 'x2p ', 'Omg_rt ', 'Th ', 'Qaer '
writeloop: do k_time = 1, Nsteps
    write (100,*) Time(k_time), x(1,k_time), x(2,k_time), x(3,k_time), x(4,k_time), &
        x(5,k_time), beta_hist(k_time), Th_hist(k_time), Qrt_hist(k_time), Qgen_hist(k_time)
end do writeloop

close (100)

deallocate (x, time, blade_dim, foil_prop, Th_hist, Qrt_hist, Qgen_hist, beta_hist)

call date_and_time(VALUES=time_values)
! write(*,*) time_values(5), time_values(6), time_values(7), time_values(8)
simtime = real(time_values(5))*3600 + real(time_values(6))*60 + real(time_values(7)) &
    + real(time_values(8))/1000 - simtime
    
    write(*,*) 'Execution time: ', simtime, 's'

    	
end program FWTsim
