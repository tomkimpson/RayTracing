module parameters
implicit none


!Define float precision
integer, parameter :: dp = selected_real_kind(15, 307)
!integer, parameter :: dp = selected_real_kind(33, 4931)

!Define useful stuff
real(kind=dp), parameter :: PI = 4.D0*ATAN(1.D0) 


!BH intrinsic parameters
real(kind=dp), parameter :: MBH = 4.31d6!BH mass in solar masses
real(kind=dp), parameter :: a= 0.90_dp !BH spin parameter Now set later


!Observer location
real(kind=dp), parameter :: r_obs = 1.0e3
real(kind=dp), parameter :: theta_obs = PI/2.0_dp

!Plasma density profile normalisation
!Set = 0 for vacuum
real(kind=dp), parameter :: N0 = 0.0e7 !1.0e8 

!Integration parameters
character(len=20), parameter :: IntegrationType = 'Backwards' !Forwards/Backwards

!Type for integrating multiple trajectories
character(len=20), parameter :: mode = 'equator' ! 'frequency', 'image', 'equator', 'single' , 'shoot'



!Parameters used if forward ray tracing
real(kind=dp), parameter :: r_init = 200.0_dp , theta_init = PI/2.0_dp, phi_init = 3.10_dp
real(kind=dp), parameter :: dir_theta = PI/2.0_dp , dir_phi = 0.010_dp !0.0_dp


!Parameters used if shooting
integer(kind=dp), parameter :: load = 1 !Do you want to load the target points?
character(len=200) :: targets_file = "/Users/tomkimpson/Data/ThesisData/MPD/targets.txt"


integer(kind=dp), parameter :: secondary_rays = 0 !Turn on/off search for seconary rays
character(len=20), parameter :: optimizer = 'PS' ! Set the optimisation method used 
!'CGD' - conjugate gradient descent. !Use with high precision
!'PS' - pattern serach. USe with low precision

!IO parameters
integer(kind=dp) :: Nrows = 6, Ncols = 1e6
character(len=200) :: IOpath = "/Users/tomkimpson/Data/ThesisData/RT/"





end module parameters
