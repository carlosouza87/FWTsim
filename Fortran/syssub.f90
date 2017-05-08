module syssub
use sysvar
use lib_array

implicit none

contains

subroutine read_inp
! Reads input files for the simulation
    
    ! System properties 
    open (unit=200, file='inp_sysprop.dat', status='old', action='read')

        ! read number of doF
        read(200,*) Ndof
        
        testNdof: if (Ndof > 2) then
                    write(*,*) 'Sorry, max. 2 rigid-body DOF so far. '
            stop
        end if testNdof
        
        ! read Mrb
        readMrb: do k1 = 1,Ndof
            read(200,*) Mrb(k1,:)
        end do readMrb		
        
        ! read Add
        readAdd: do k1 = 1,Ndof
            read(200,*) Add(k1,:)
        end do readAdd

        ! read Dl
        readDl: do k1 = 1,Ndof
            read(200,*) Dl(k1,:)
        end do readDl
        
       ! read Dq
        readDq: do k1 = 1,Ndof
            read(200,*) Dq(k1,:)
        end do readDq
        
        ! read Chs
        readChs: do k1 = 1,Ndof
            read(200,*) Chs(k1,:)
        end do readChs
        
        ! read Cmr
        readCmr: do k1 = 1,Ndof
            read(200,*) Cmr(k1,:)
        end do readCmr	

        ! ! read Mrb
        ! read(200,*) Mrb(1,:)
        ! read(200,*) Mrb(2,:)
        
        ! ! read Add
        ! read(200,*) Add(1,:)
        ! read(200,*) Add(2,:)
        
        ! ! read Dl
        ! read(200,*) Dl(1,:)
        ! read(200,*) Dl(2,:)
        
        ! ! read Dq
        ! read(200,*) Dq(1,:)
        ! read(200,*) Dq(2,:)
        
        ! ! read Chs
        ! read(200,*) Chs(1,:)
        ! read(200,*) Chs(2,:)
        
        ! ! read Cmr
        ! read(200,*) Cmr(1,:)
        ! read(200,*) Cmr(2,:)

    close (200)

    ! Simulation parameters
    open (unit=120, file='inp_simpar.dat', status='old', action='read')

        ! read ti
        read(120,*) ti
        
        ! read tf
        read(120,*) tf
            
        ! read dt
        read(120,*) dt	
        
        ! read t_clutch
        read(120,*) t_clutch	
        
        ! read eta0
        readeta0: do k1 = 1,Ndof
            read(120,*) eta0(k1)
        end do readeta0
        
        ! read nu0
        readnu0: do k1 = 1,Ndof
            read(120,*) nu0(k1)
        end do readnu0
        
        ! read Uinf
        read(120,*) Uinf
        
        ! Read initial Omg_rt
        read(120,*) Omg_rt

    close (120)

    ! Rotor properties
    open (unit=130, file='inp_rotorprop.dat', status='old', action='read')

        ! Read number of foils
        read(130,*) Nfoils
        
        ! Read all foil files, and store the length of each of them	
        foil_length: do k1 = 1,Nfoils
        
            read(130,*) foilfilename(k1)	
        
            open (unit=140, file=foilfilename(k1), status='old', action='read')
                read(140,*) Ninc(k1)
            close (140)
        end do foil_length
        
        ! Determine the maximum	length of all foil descriptions
        Max_inc = maxval(Ninc)
        
        ! Read Bl
        read(130,*) Bl
        
        ! Read initial value for beta
        read(130,*) beta
        
        ! Read Rtip
        read(130,*) Rtip
        
        ! Read Rhub
        read(130,*) Rhub
        
        ! Read Zhub
        read(130,*) Zhub
        
        ! Read Nelem
        read(130,*) Nelem
        
        ! Read matrix with blade elements properties
        allocate (Blade_dim(Nelem,5))  
        
        read_blddim: do k1 = 1,Nelem
            read(130,*) Blade_dim(k1,:)	  	
        end do read_blddim		
    close (130)

    ! 3D matrix with foil data
    allocate( Foil_prop(Max_inc,4,Nfoils))

    read_foildata: do k1 = 1,Nfoils
        open (unit=150, file=foilfilename(k1), status='old', action='read')
            read(150,*) 
            write_currfoil: do k2 = 1,Ninc(k1)
                read(150,*) Foil_prop(k2,:,k1)
            end do write_currfoil
        close (150)
    end do read_foildata
    
    ! Generator properties
    open (unit=160, file='inp_generprop.dat', status='old', action='read')
    
        ! Read Igen
        read(160,*) Igen
        
        ! Read Ngr
        read(160,*) Ngr
        
        ! Read PC_KI
        read(160,*) PC_KI
        
        ! Read PC_KK
        read(160,*) PC_KK
        
        ! Read PC_KP
        read(160,*) PC_KP
        
        ! Read PC_MaxPit
        read(160,*) PC_MaxPit
        
        ! Read PC_MaxRat
        read(160,*) PC_MaxRat
        
        ! Read PC_MinPit
        read(160,*) PC_MinPit
        
        ! Read PC_RefSpd
        read(160,*) PC_RefSpd
        
        ! Read VS_CtInSp
        read(160,*) VS_CtInSp
        
        ! Read VS_DT
        read(160,*) VS_DT
        
        ! Read VS_MaxRat
        read(160,*) VS_MaxRat    

        ! Read VS_MaxTq
        read(160,*) VS_MaxTq
        
        ! Read VS_Rgn2K
        read(160,*) VS_Rgn2K
        
        ! Read VS_Rgn2Sp
        read(160,*) VS_Rgn2Sp
        
        ! Read VS_Rgn3MP
        read(160,*) VS_Rgn3MP
        
        ! Read VS_RtGnSp
        read(160,*) VS_RtGnSp
        
        ! Read VS_RtPwr
        read(160,*) VS_RtPwr
        
        ! Read VS_SlPc
        read(160,*) VS_SlPc       
        
    close(160)
    
end subroutine read_inp

subroutine RK4th(x0,t,dt,x)
! Performs numerical integration with a 4-th order Runge-Kutta algorithm

    implicit none

    integer                                   :: iter        ! Current RK iteration
    integer                                   :: N           ! Order of system

    real(8)                                   :: t           ! Time instant
    real(8)                                   :: dt          ! Time step
    real(8), dimension(:)                     :: x0          ! Current state 
    real(8), dimension(:)                     :: x           ! Integrated state
    real(8), dimension(:), allocatable        :: r1          ! State of first iteration
    real(8), dimension(:), allocatable        :: r2          ! State of second iteration
    real(8), dimension(:), allocatable        :: r3          ! State of third iteration
    real(8), dimension(:), allocatable        :: r4          ! State of fourth iteration

	N = size(x0,1)
	
	allocate (r1(N), r2(n), r3(N), r4(N))
	
    iter = 1	
    call sysdyn(x0(1:N),t,dt,iter,r1(1:N))
	
    iter = 2
    call sysdyn(x0(1:N)+r1(1:N)*dt/2,t+dt/2,dt,iter,r2(1:N))
	
    iter = 3 
    call sysdyn(x0(1:N)+r2(1:N)*dt/2,t+dt/2,dt,iter,r3(1:N))
	
    iter = 4
    call sysdyn(x0(1:N)+r3(1:N)*dt,t+dt,dt,iter,r4(1:N))
	
	x(1:N) = x0(1:N) + (dt/6)*(r1(1:N)+2*r2(1:N)+2*r3(1:N)+r4(1:N))
	
	deallocate (r1, r2, r3, r4)

end subroutine RK4th

subroutine sysdyn(x,t,dt,iter,x_p)
! Calculates the states derivatives, based on the system dynamics

    implicit none

	integer                                   :: iter        ! Current RK iteration
    ! integer                                   :: idx         ! Current index in state vector
	
    real(8), dimension(2,1)                   :: Brhs        ! Right-hand-side vector of eq. of motions
    real(8)                                   :: dt          ! Time step
    real(8), dimension(2,1)                   :: eta         ! Position state
    real(8), dimension(2,1)                   :: eta_p       ! Time derivative of position state
    real(8), dimension(2,1)                   :: Fwind       ! Vector with wind loads
    real(8)                                   :: Idt         ! Drivetrain inertia about rotor axis [kg.m^2]
    real(8), dimension(2,2)                   :: Minv        ! Inverse of Mrb+Add
	real(8), dimension(2,1)                   :: nu          ! Velocity state    
	real(8), dimension(2,1)                   :: nu_p        ! Time derivative of velocity state
    real(8)                                   :: Omg_gen     ! Generator rotational speed [rad/s]
    real(8)                                   :: Omg_rt      ! Rotor rotational speed [rad/s]
    real(8)                                   :: Omg_rt_d    ! Time derivative of rotor rotational speed [rad/s^2]
    real(8)                                   :: Qaero       ! Aerodynamic torque [N.m]
    real(8)                                   :: Qgen        ! Generator torque [N.m]
	real(8)                                   :: t           ! Time instant        
    ! real(8)                                   :: Qm          ! Rotor moment
    ! real(8)                                   :: Th          ! Rotor thrust
    real(8)                                   :: Uhub        ! Hub velocity
	real(8), dimension(:)                     :: x           ! Current state
	real(8), dimension(:)                     :: x_p         ! State to integrate

	
	
	! real(8), dimension(2,2)                        :: J           ! Transformation matrix (body-fixed to inertial system)
	! real(8), dimension(:)                          :: x           ! Current state
	! real(8), dimension(:)                          :: x_p         ! State to integrate
	! real(8), dimension(:), allocatable             :: eta         ! Position state
	! real(8), dimension(:), allocatable             :: nu          ! Velocity state
	! real(8), dimension(:), allocatable             :: eta_p       ! Time derivative of position state
	! real(8), dimension(:), allocatable             :: nu_p        ! Time derivative of velocity state
	! real(8), dimension(:,:), allocatable           :: Brhs        ! Right-hand-side matrix of eq. of motions
	! real(8), dimension(:,:), allocatable   

	! allocate (eta(N/2), nu(N/2), eta_p(N/2), nu_p(N/2),Brhs(N/2,N/2),Minv(N/2,N/2))
	
    ! ! Initialize idx
	! idx = 0
    
    ! Initialize position and velocity vectors from input state
	eta(:,1) = x(1:2)
    nu(:,1) = x(3:4)
    ! Omg_rt = x(5)
    ! **************************************
    ! GAMBIARRA !!!
    ! **************************************
    Omg_rt = 1.3
    
    check_iter: if (iter == 1) then           
        ! Call controller for calculating blade pitch and generator torque        
        Omg_gen = Omg_rt*Ngr
        Qgen = Qgen_hist(k_time-1)         
        call control(beta,Omg_gen,Qgen)
        ! **************************************
        ! GAMBIARRA !!!
        ! **************************************
        beta = 0
        Qgen_hist(k_time) = Qgen
        beta_hist(k_time) = beta
        
        ! Calculate rotor loads
        Uhub = nu(1,1) + nu(2,1)*Zhub          
        call BEM_ning(Uhub,Th_hist(k_time),Qrt_hist(k_time))      
        Qaero = Qrt_hist(k_time)      
    else
        Qgen = Qgen_hist(k_time)        
        Qaero = Qrt_hist(k_time)
    end if check_iter

    
    Fwind(1,1) = -Th_hist(k_time)
    Fwind(2,1) = -Th_hist(k_time)*Zhub  
    
    ! Inversion of Mrb + Add
	Minv = Mrb + Add
	call matinv2(Minv)
	
	! Calculate right-hand-side of equations of motions
	Brhs = -matmul(Dl,nu)-matmul(Dq,abs(nu)*nu)-matmul(Cmr,eta)-matmul(Chs,eta)+Fwind
    
    ! Calculate the derivative of eta, eta_p = nu
    eta_p = nu
    ! Calculate the derivative of nu, nu_p = Minv*Brhs
	nu_p = matmul(Minv,Brhs)
    
    check_clutch: if (t <= t_clutch) then
        eta_p(:,1) = eta_p(:,1)*0
        nu_p(:,1) = nu_p(:,1)*0
    end if check_clutch    
 
    ! Evaluate drivetrain dynamics
    Idt = Irt + (Ngr**2)*Igen ! Total drivetrain inertia
    Omg_rt_d = (Qaero-Ngr*Qgen)/Idt ! Derivative of rotor speed    

    ! Output time-derivative of states
    x_p(1:2) = eta_p(:,1)
    x_p(3:4) = nu_p(:,1)
    ! x_p(5) = Omg_rt_d   
    ! **************************************
    ! GAMBIARRA !!!
    ! **************************************    
    x_p(5) = 0    
	
	! deallocate (eta, nu, eta_p, nu_p)

end subroutine sysdyn

subroutine BEM_ning(Uhub,Th,Qm)
! Blade Element Momentum theory following Ning (2013)

    use sysvar
    use lib_array

    implicit none
    
    integer                                   :: k_elem         ! Counter for elements along the blade []
    
    real(8), parameter                        :: tlr = 1e-6     ! Tolerance []
    
    real(8)                                   :: a, a_p         ! Axial and tangential induction factors []
    real(8)                                   :: alpha          ! Angle of attack [deg]
    real(8)                                   :: ch             ! Chord length [m]
    real(8)                                   :: cd, cl         ! Drag and lift coefficients []
    real(8)                                   :: cn, ct         ! Normal and tangential force coefficients []
    real(8)                                   :: dQm, dTh       ! Incremental moment and thrust [N]
    real(8)                                   :: dr             ! Section length [m]
    real(8)                                   :: F              ! Tip/hub loss correction factor []
    real(8)                                   :: Fhub, Ftip     ! Hub and tip loss correction factors []
    real(8)                                   :: gamma1, gamma2, gamma3 ! Auxiliary factors
    real(8)                                   :: kappa, kappa_p ! Convenience parameters for axial and tangential inductions, as defined in Ning (2013)    
    ! real(8)                                   :: lbd_tip        ! Tip speed ratio []    
    real(8)                                   :: phi            ! Inflow angle [rad]
    real(8)                                   :: Qm             ! Moment [N]     
    real(8)                                   :: Th             ! Thrust [N]
    real(8)                                   :: Uhub           ! Hub velocity [m/s]
    real(8)                                   :: Urel           ! Relative flow velocity [m/s]
    
    ! Initialize Th and Qm 
    Th = 0
    Qm = 0    
    
    ! Relative wind velocity [m/s]
    Urel = Uinf + Uhub
    ! Urel = Uinf
    
    do_along_elem: do k_elem = 1,Nelem
        r = Blade_dim(k_elem,1)
        twist = Blade_dim(k_elem,2)*PI/180
        dr = Blade_dim(k_elem,3)
        ch = Blade_dim(k_elem,4)
        foil_id = int(Blade_dim(k_elem,5))
        
        lbd_r = Omg_rt*r/Urel ! Local speed ratio for current element
        sigma_p = Bl*ch/(2*PI*r) ! Solidity for current element
        
        phi = zero(tlr,PI-tlr,macheps,tlr/1000,f_ning) ! Inflow angle
        
        ! Calculate tip and hub losses    
        Ftip = 2/PI * acos(exp(-(Bl/2*(Rtip-r)/(r*sin(phi))))) ! Tip loss correction (Prandt)
        Fhub = 2/PI * acos(exp(-(Bl/2*(r-Rhub)/(r*sin(phi))))) ! Hub loss correction
        F =  Ftip*Fhub;	! Total loss
                
        ! Lift and drag coefficients for calculated phi
        alpha = (phi - (twist+beta))*180/pi
        cl = interp1d(Foil_prop(1:Ninc(foil_id),1,foil_id),Foil_prop(1:Ninc(foil_id),2,foil_id),alpha)
	    cd = interp1d(Foil_prop(1:Ninc(foil_id),1,foil_id),Foil_prop(1:Ninc(foil_id),3,foil_id),alpha)
    
        ! Normal and tangential coefficients for calculated phi
        cn = cl*cos(phi) + cd*sin(phi)
	    ct = cl*sin(phi) - cd*cos(phi)
        
        ! Convenience factors (see Ning, (2013))
        kappa  = (sigma_p*cn)/(4*F*sin(phi)**2)
        kappa_p = (sigma_p*ct)/(4*F*sin(phi)*cos(phi))
      
        ! Check if kappa exceeds limit of 2/3 (corresponding to a = 0.4)
        test_kappa: if (kappa < 2./3) then
            a = kappa/(1+kappa) ! Axial induction factor, ordinary BEM theory
        else	
            avoid_sing: if (kappa == 25./(18*F)) then
                kappa = kappa + 1e-5
            end if avoid_sing
            gamma1 = 2*F*kappa - (10./9-F)
            gamma2 = 2*F*kappa - F*(4./3-F) 
            gamma3 = 2*F*kappa - (25./9-2*F)
            a = (gamma1-sqrt(gamma2))/gamma3 ! Axial induction factor, with modifications proposed in Ning (2013)
        end if test_kappa	
        
        a_p = kappa_p/(1-kappa_p) ! Tangential induction factor
        ! a_p = (-1+4*F*sin(phi)*cos(phi)/(sigma_p*ct))**(-1) ! Tangential induction factor
        
        ! Calculate incremental thrust and moment
        dTh = 4*PI*r*rho_a*Urel**2*(1-a)*a*dr        
        dQm = 4*PI*r**3*rho_a*Urel*Omg_rt*(1-a)*a_p*dr
        
        ! Add increments to total thrust and moment
        Th = Th + dTh
        Qm = Qm + dQm       
        	
    end do do_along_elem

    
! Reference:
! Ning, S. A. - "A simple solution method for the blade element momentum equations with guaranteed convergence."
! Wind Energy, 2013 , 17 , 1327-1345    
end subroutine BEM_ning

function f_ning (phi)
! Calculates function whose zero corresponds to the actual phi, according to Ning (2013).

    implicit none
    
    real(8)                                   :: a              ! Axial induction factor []
    real(8)                                   :: alpha          ! Angle of attack [deg]
    real(8)                                   :: cd, cl         ! Drag and lift coefficients []
    real(8)                                   :: cn, ct         ! Normal and tangential force coefficients []
    real(8)                                   :: F              ! Tip/hub loss correction factor []
    real(8)                                   :: Fhub, Ftip     ! Hub and tip loss correction factors []
    real(8)                                   :: f_ning         ! Function output
    real(8)                                   :: gamma1, gamma2, gamma3 ! Auxiliary factors
    real(8)                                   :: kappa, kappa_p ! Convenience parameters for axial and tangential inductions, as defined in Ning (2013)
	real(8)                                   :: phi            ! Inflow angle    
    
    alpha = (phi - (twist+beta))*180/pi
			
	cl = interp1d(Foil_prop(1:Ninc(foil_id),1,foil_id),Foil_prop(1:Ninc(foil_id),2,foil_id),alpha)
	cd = interp1d(Foil_prop(1:Ninc(foil_id),1,foil_id),Foil_prop(1:Ninc(foil_id),3,foil_id),alpha)
    
    cn = cl*cos(phi) + cd*sin(phi)
	ct = cl*sin(phi) - cd*cos(phi)
    
   	Ftip = 2/PI * acos(exp(-(Bl/2*(Rtip-r)/(r*sin(phi))))) ! Tip loss correction (Prandt)
	Fhub = 2/PI * acos(exp(-(Bl/2*(r-Rhub)/(r*sin(phi))))) ! Hub loss correction

    F =  Ftip*Fhub;	
    
    kappa  = (sigma_p*cn)/(4*F*sin(phi)**2)
    kappa_p = (sigma_p*ct)/(4*F*sin(phi)*cos(phi))  	
        
	test_kappa: if (kappa < 2./3) then    
	    a = kappa/(1+kappa)
    else
	    avoid_sing: if (kappa == 25./(18*F)) then
		    kappa = kappa + 1e-5
	    end if avoid_sing	
	    gamma1 = 2*F*kappa - (10./9-F)
		gamma2 = 2*F*kappa - F*(4./3-F) 
		gamma3 = 2*F*kappa - (25./9-2*F)
	    a = (gamma1-sqrt(gamma2))/gamma3
    end if test_kappa	
    
    test_phi: if (phi > 0.) then
        f_ning = sin(phi)/(1-a) - cos(phi)*(1-kappa_p)/lbd_r 
    else
        f_ning = sin(phi)*(1-kappa) - cos(phi)*(1-kappa_p)/lbd_r
    end if test_phi

! Reference:
! Ning, S. A. - "A simple solution method for the blade element momentum equations with guaranteed convergence."
! Wind Energy, 2013 , 17 , 1327-1345      
end function f_ning

subroutine control(BlPitch,GenSpeed,GenTrq)
! Control system, adapted from DISCON (FAST)

    implicit none

    real(8)                                       :: BlPitch        ! Blade pitch angle [rad]
    real(8)                                       :: GenSpeed       ! Current  HSS (generator) speed [rad/s]
    real(8)                                       :: GenSpeedF      ! Filtered HSS (generator) speed [rad/s]
    real(8)                                       :: GenTrq         ! Generator torque [N.m]
    real(8)                                       :: GK             ! Gain correction factor for scheduled PI controller []
    real(8), save                                 :: IntSpdErr      ! Integral of generator speed error [rad]
    real(8)                                       :: LastGenTrq     ! Commanded electrical generator torque the last time the controller was called [N.m]
    real(8)                                       :: PitCom         ! Commanded blade pitch angle [rad]
    real(8)                                       :: PitComI        ! Integral term of command pitch [rad]
    real(8)                                       :: PitComP        ! Proportional term of command pitch [rad]
    real(8)                                       :: PitComT        ! Total command pitch based on the sum of the proportional and integral terms [rad]
    real(8)                                       :: PitRate        ! Pitch rates of each blade based on the current pitch angles and current pitch command [rad/s]
    real(8)                                       :: SpdErr         ! Current speed error [rad/s]
    real(8)                                       :: TrqRate        ! Torque rate based on the current and last torque commands [N.m/s]
    real(8), save                                 :: VS_Slope15     ! Torque/speed slope of region 1 1/2 cut-in torque ramp [N.m.s/rad]
    real(8), save                                 :: VS_Slope25     ! Torque/speed slope of region 2 1/2 induction generator [N.m.s/rad]
    real(8), save                                 :: VS_SySp        ! Synchronous speed of region 2 1/2 induction generator [rad/s]
    real(8), save                                 :: VS_TrGnSp      ! Transitional generator speed (HSS side) between regions 2 and 2 1/2 [rad/s]

    VS_SySp    = VS_RtGnSp/( 1.0 +  0.01*VS_SlPc )
    VS_Slope15 = ( VS_Rgn2K*VS_Rgn2Sp*VS_Rgn2Sp )/( VS_Rgn2Sp - VS_CtInSp )
    VS_Slope25 = ( VS_RtPwr/PC_RefSpd           )/( VS_RtGnSp - VS_SySp   )
    if ( VS_Rgn2K == 0.0 )  then  ! .TRUE. if the Region 2 torque is flat, and thus, the denominator in the ELSE condition is zero
       VS_TrGnSp = VS_SySp
    else                          ! .TRUE. if the Region 2 torque is quadratic with speed
        VS_TrGnSp = ( VS_Slope25 - sqrt( VS_Slope25*( VS_Slope25 - 4.0*VS_Rgn2K*VS_SySp ) ) )/( 2.0*VS_Rgn2K )
    endif

    PitCom     = BlPitch                      ! This will ensure that the variable speed controller picks the correct control region and the pitch controller pickes the correct gain on the first call
    GK         = 1.0/( 1.0 + PitCom/PC_KK )   ! This will ensure that the pitch angle is unchanged if the initial SpdErr is zero
    IntSpdErr  = PitCom/( GK*PC_KI )          ! This will ensure that the pitch angle is unchanged if the initial SpdErr is zero

    LastGenTrq = GenTrq  ! Initialize the value of LastGenTrq 
    GenSpeedF = GenSpeed ! No filtering for now   

    if ( (   GenSpeedF >= VS_RtGnSp ) .OR. (  PitCom >= VS_Rgn3MP ) )  then ! We are in region 3 - power is constant
        GenTrq = VS_RtPwr/PC_RefSpd   !JASON:MAKE REGION 3 CONSTANT TORQUE INSTEAD OF CONSTANT POWER FOR HYWIND:    
        ! GenTrq = VS_RtPwr/GenSpeedF
    elseif ( GenSpeedF <= VS_CtInSp )  then                                    ! We are in region 1 - torque is zero
        GenTrq = 0.0
    elseif ( GenSpeedF <  VS_Rgn2Sp )  then                                    ! We are in region 1 1/2 - linear ramp in torque from zero to optimal
        GenTrq = VS_Slope15*( GenSpeedF - VS_CtInSp )
    elseif ( GenSpeedF <  VS_TrGnSp )  then                                    ! We are in region 2 - optimal torque is proportional to the square of the generator speed
        GenTrq = VS_Rgn2K*GenSpeedF*GenSpeedF
    else                                                                       ! We are in region 2 1/2 - simple induction generator transition region
        GenTrq = VS_Slope25*( GenSpeedF - VS_SySp   )
    endif
    
    ! Saturate the commanded torque using the maximum torque limit:
    
    GenTrq  = min( GenTrq, VS_MaxTq)   ! Saturate the command using the maximum torque limit
    
    ! Saturate the commanded torque using the torque rate limit:
    TrqRate = ( GenTrq*(1+eps) - LastGenTrq )/dt               ! Torque rate (unsaturated)
    TrqRate = min( max( TrqRate, -VS_maxRat ), VS_maxRat )   ! Saturate the torque rate using its maximum absolute value
    GenTrq  = LastGenTrq + TrqRate*dt                  ! Saturate the command using the torque rate limit    

    
    ! ****************************************************************************************************
    ! Blade pitch controller
    ! ****************************************************************************************************
    ! Compute the gain scheduling correction factor based on the previously
    !   commanded pitch angle for blade 1:

    GK = 1.0/( 1.0 + PitCom/PC_KK )


    ! Compute the current speed error and its integral w.r.t. time; saturate the
    !   integral term using the pitch angle limits:

    SpdErr    = GenSpeedF - PC_RefSpd                                 ! Current speed error
    IntSpdErr = IntSpdErr + SpdErr*dt                           ! Current integral of speed error w.r.t. time
    IntSpdErr = min( max( IntSpdErr, PC_minPit/( GK*PC_KI ) ), &
                                   PC_maxPit/( GK*PC_KI )      )    ! Saturate the integral term using the pitch angle limits, converted to integral speed error limits

    ! Compute the pitch commands associated with the proportional and integral
    !   gains:

    PitComP   = GK*PC_KP*   SpdErr                                    ! Proportional term
    PitComI   = GK*PC_KI*IntSpdErr                                    ! Integral term (saturated)

    ! Superimpose the individual commands to get the total pitch command;
    !   saturate the overall command using the pitch angle limits:

    PitComT   = PitComP + PitComI                                     ! Overall command (unsaturated)
    PitComT   = min( max( PitComT, PC_minPit ), PC_maxPit )           ! Saturate the overall command using the pitch angle limits

    ! Saturate the overall commanded pitch using the pitch rate limit:
    PitRate = (PitComT - BlPitch)/dt
    PitRate = min(max(PitRate,-PC_maxRat),PC_maxRat)
    PitCom = BlPitch + PitRate*dt

    BlPitch = PitCom


end subroutine control

subroutine matinv2(A) 
    ! Performs a direct calculation of the inverse of a 2×2 matrix. 
	! Source: http://fortranwiki.org/fortran/show/Matrix+inversion
	
    real(8)                 :: A(2,2)   !! Matrix
    real(8)                 :: B(2,2)   !! Temp matrix
    real(8)                 :: detinv   !! Determinant

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
	
	A(:,:) = B(:,:)
end subroutine

function zero ( a, b, machep, t, f )

!*****************************************************************************80
!
!! ZERO seeks the root of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the zero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2013
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the change of 
!    sign interval.
!
!    Input, real ( kind = 8 ) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real ( kind = 8 ) T, a positive error tolerance.
!
!    Input, external real ( kind = 8 ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose zero is being sought.
!
!    Output, real ( kind = 8 ) ZERO, the estimated value of a zero of
!    the function F.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fc
  real ( kind = 8 ) m
  real ( kind = 8 ) machep
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) zero
!
!  Make local copies of A and B.
!
  sa = a
  sb = b
  fa = f ( sa )
  fb = f ( sb )

  c = sa
  fc = fa
  e = sb - sa
  d = e

  do

    if ( abs ( fc ) < abs ( fb ) ) then

      sa = sb
      sb = c
      c = sa
      fa = fb
      fb = fc
      fc = fa

    end if

    tol = 2.0D+00 * machep * abs ( sb ) + t
    m = 0.5D+00 * ( c - sb )

    if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
      exit
    end if

    if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

      e = m
      d = e

    else

      s = fb / fa

      if ( sa == c ) then

        p = 2.0D+00 * m * s
        q = 1.0D+00 - s

      else

        q = fa / fc
        r = fb / fc
        p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
        q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

      end if

      if ( 0.0D+00 < p ) then
        q = - q
      else
        p = - p
      end if

      s = e
      e = d

      if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
        p < abs ( 0.5D+00 * s * q ) ) then
        d = p / q
      else
        e = m
        d = e
      end if

    end if

    sa = sb
    fa = fb

    if ( tol < abs ( d ) ) then
      sb = sb + d
    else if ( 0.0D+00 < m ) then
      sb = sb + tol
    else
      sb = sb - tol
    end if

    fb = f ( sb )

    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if

  end do

  zero = sb

  return
end function zero

end module syssub