module sysvar

implicit none

character(256), dimension(20)                :: foilfilename  ! Names of foil data files

integer                                      :: Bl            ! Number of rotor blades 
integer                                      :: foil_id       ! Foil identifier
integer                                      :: i             ! Counter for general purposes
integer                                      :: k_time        ! Counter for time 
integer                                      :: k1, k2, k3    ! Count variables for DO constructs
integer                                      :: Max_inc       ! Maximum number of incidence directions over all foils considered
integer                                      :: Ndof          ! Number of rigid-body degrees of freedom
integer                                      :: Nelem         ! Number of elements on blades
integer                                      :: Nfoils        ! Number of foils for describing blade profile
integer, dimension(20)                       :: Ninc          ! Number of incidence directions for foil aerodynamics properties
integer                                      :: Nsteps        ! Number of simulation time steps
integer                                      :: Nsys          ! Order of system

real(8)                                      :: Ngear = 97.0   ! Rotor-generator transmission ratio [] - TO BE DEFINED IN INPUT FILE
real(8)                                      :: mu_a = 1.77e-5 ! Dynamic viscosity of air [kg.s/m]
real(8)                                      :: PI = 4.D0*datan(1.D0) ! Pi []
real(8)                                      :: macheps = epsilon(PI) ! Machine tolerance
real(8)                                      :: rho_a = 1.225  ! Air density (kg/m^3)
real(8)                                      :: rho_w = 1025.0 ! Water density (kg/m^3)
real(8)                                      :: Uinf = 0.1 ! Wind velocity [m/s]

real(8), dimension(2,2)                      :: Add           ! Constant-frequency added mass matrix
real(8)                                      :: beta          ! Blade pitch angle [rad]
real(8), dimension(:), allocatable           :: beta_hist     ! Blade pitch angle (variable for storage) [rad]
real(8), dimension(:,:), allocatable         :: Blade_dim     ! Matrix with balde dimensions
real(8), dimension(2,2)                      :: Chs           ! Hydrostatic stiffness matrix
real(8), dimension(2,2)                      :: Cmr           ! Mooring stiffness matrix
real(8), parameter                           :: eps = epsilon(eps) ! Precision []
real(8), dimension(2,2)                      :: Dl            ! Linear viscous damping matrix
real(8), dimension(2,2)                      :: Dq            ! Quadratic viscous damping matrix
real(8)                                      :: dt            ! Time step [s]
real(8), dimension(2)                        :: eta0          ! Initial positions [m;rad]
real(8), dimension(:,:,:), allocatable       :: Foil_prop     ! 3D matrix for foil aerodynamic properties
! real(8), dimension(2,1)                      :: Fwind       ! Vector with wind loads
real(8)                                      :: Igen          ! Generator inertia about HSS [kg.m^2]
real(8)                                      :: Irt           ! Blades + hub inertia about rotor axis [kg.m^2]
real(8)                                      :: lbd_r         ! Local speed ratio []
real(8), dimension(2,2)                      :: Mrb           ! Rigid-body inertia matrix
real(8)                                      :: Ngr           ! Gear ratio []
real(8), dimension(2)                        :: nu0           ! Initial velocities [m/s;rad/s]
real(8)                                      :: Omg_rt        ! Rotor angular velocity [rad/s] - TO BE CALCULATED
! real(8), dimension(:), allocatable           :: Omg_hist      ! Rotor angular velocity (variable for storage) [rad/s]
real(8)                                      :: PC_KI         ! Integral gain for pitch controller [s]
real(8)                                      :: PC_KK         ! Pitch angle were the the derivative of the aerodynamic power w.r.t. pitch has increased by a factor of two relative to the derivative at rated pitch (zero) [rad]
real(8)                                      :: PC_KP         ! Proportional gain for pitch controller at rated pitch (zero) [s]
real(8)                                      :: PC_MaxPit     ! Maximum pitch setting in pitch controller [rad]
real(8)                                      :: PC_MaxRat     ! Maximum pitch  rate (in absolute value) in pitch  controller [rad/s]
real(8)                                      :: PC_MinPit     ! Minimum pitch setting in pitch controller [rad]
real(8)                                      :: PC_RefSpd     ! Desired (reference) HSS speed for pitch controller [rad/s]
real(8), dimension(:), allocatable           :: Qgen_hist     ! Generator torque (variable for storage) [N]
real(8), dimension(:), allocatable           :: Qrt_hist      ! Rotor torque (variable for storage) [N]
real(8)                                      :: r             ! Local blade radius [m]
real(8)                                      :: Rhub          ! Hub radius   [m]
real(8)                                      :: Rtip          ! Rotor radius [m]
real(8)                                      :: sigma_p       ! Blade element solidity []
real(8)                                      :: t             ! Current time during simulation [s]
real(8)                                      :: t_clutch      ! Time for clutching system dynamics [s]
real(8)                                      :: ti, tf        ! Initial and final time [s]
real(8), dimension(:), allocatable           :: Th_hist       ! Rotor thrust (variable for storage) [N]
real(8)                                      :: twist         ! Blade element twist angle [rad]
real(8)                                      :: VS_CtInSp     ! Transitional generator speed (HSS side) between regions 1 and 1 1/2 [rad/s]
real(8)                                      :: VS_DT         ! Communication interval for torque controller [s]
real(8)                                      :: VS_MaxRat     ! Maximum torque rate (in absolute value) in torque controller [N.m/s]
real(8)                                      :: VS_MaxTq      ! Maximum generator torque in Region 3 (HSS side) [N.m]
real(8)                                      :: VS_Rgn2K      ! Generator torque constant in Region 2 (HSS side) [N.m.s^2/rad^2]
real(8)                                      :: VS_Rgn2Sp     ! Transitional generator speed (HSS side) between regions 1 1/2 and 2 [rad/s]
real(8)                                      :: VS_Rgn3MP     ! Minimum pitch angle at which the torque is computed as if we are in region 3 regardless of the generator speed [rad]
real(8)                                      :: VS_RtGnSp     ! Rated generator speed (HSS side) [rad/s]
real(8)                                      :: VS_RtPwr      ! Rated generator generator power in Region 3 [W]
real(8)                                      :: VS_SlPc       ! Rated generator slip percentage in Region 2 1/2 [%]
real(8)                                      :: Zhub          ! Hub heigth [m]

! real(8), dimension(:), allocatable           :: eta0          ! Initial positions [m;rad]
! real(8), dimension(:), allocatable           :: nu0           ! Initial velocities [m/s;rad/s]
! real(8), dimension(:,:), allocatable         :: Mrb           ! Rigid-body inertia matrix
! real(8), dimension(:,:), allocatable         :: Add           ! Constant-frequency added mass matrix
! real(8), dimension(:,:), allocatable         :: Dl            ! Linear viscous damping matrix
! real(8), dimension(:,:), allocatable         :: Dq            ! Quadratic viscous damping matrix
! real(8), dimension(:,:), allocatable         :: Chs           ! Hydrostatic stiffness matrix
! real(8), dimension(:,:), allocatable         :: Cmr           ! Mooring stiffness matrix

end module sysvar

