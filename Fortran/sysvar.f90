module sysvar

implicit none

character(256), dimension(20)                :: foilfilename  ! Names of foil data files

integer                                      :: Bl            ! Number of rotor blades 
integer                                      :: foil_id       ! Foil identifier
integer                                      :: i             ! Counter for general purposes
integer                                      :: k1, k2, k3    ! Count variables for DO constructs
integer                                      :: Max_inc       ! Maximum number of incidence directions over all foils considered
integer                                      :: Ndof          ! Number of degrees of freedom
integer                                      :: Nelem         ! Number of elements on blades
integer                                      :: Nfoils        ! Number of foils for describing blade profile
integer, dimension(20)                       :: Ninc          ! Number of incidence directions for foil aerodynamics properties
integer                                      :: Nsteps        ! Number of simulation time steps
integer                                      :: Nsys          ! Order of system

real(8)                                      :: mu_a = 1.77e-5 ! Dynamic viscosity of air [kg.s/m]
real(8)                                      :: PI = 4.D0*datan(1.D0) ! Pi []
real(8)                                      :: macheps = epsilon(PI) ! Machine tolerance
real(8)                                      :: Omg = 1.26 ! Rotor angular velocity - TO BE CALCULATED
real(8)                                      :: rho_a = 1.18  ! Air density (kg/m^3)
real(8)                                      :: rho_w = 1025.0 ! Water density (kg/m^3)
real(8)                                      :: Uinf = 11.4 ! Wind velocity - TO BE DECLARED IN INPUT FILE

real(8), dimension(2,2)                      :: Add           ! Constant-frequency added mass matrix
real(8)                                      :: beta          ! Blade pitch angle [rad]
real(8), dimension(:,:), allocatable         :: Blade_dim     ! Matrix with balde dimensions
real(8), dimension(2,2)                      :: Chs           ! Hydrostatic stiffness matrix
real(8), dimension(2,2)                      :: Cmr           ! Mooring stiffness matrix
real(8), dimension(2,2)                      :: Dl            ! Linear viscous damping matrix
real(8), dimension(2,2)                      :: Dq            ! Quadratic viscous damping matrix
real(8)                                      :: dt            ! Time step [s]
real(8), dimension(2)                        :: eta0          ! Initial positions [m;rad]
real(8), dimension(:,:,:), allocatable       :: Foil_prop     ! 3D matrix for foil aerodynamic properties
real(8), dimension(2,1)                      :: Fwind       ! Vector with wind loads
real(8)                                      :: lbd_r         ! Local speed ratio []
real(8), dimension(2,2)                      :: Mrb           ! Rigid-body inertia matrix
real(8), dimension(2)                        :: nu0           ! Initial velocities [m/s;rad/s]
!real(8)                                      :: Omg           ! Rotor angular velocity [rad/s]
real(8)                                      :: r             ! Local blade radius [m]
real(8)                                      :: Rhub          ! Hub radius   [m]
real(8)                                      :: Rtip          ! Rotor radius [m]
real(8)                                      :: sigma_p       ! Blade element solidity []
real(8)                                      :: t             ! Current time during simulation [s]
real(8)                                      :: ti, tf        ! Initial and final time [s]
real(8)                                      :: twist         ! Blade element twist angle [rad]
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

