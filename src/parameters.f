module parameters
implicit none


!Define float precision
!integer, parameter :: dp = selected_real_kind(15, 307)
integer, parameter :: dp = selected_real_kind(33, 4931)

!Define useful stuff
real(kind=dp), parameter :: PI = 4.D0*ATAN(1.D0) 


!BH intrinsic parameters
real(kind=dp), parameter :: MBH = 4.31d6!BH mass in solar masses
real(kind=dp), parameter :: a= +0.60_dp !BH spin parameter Now set later


!Observer location
real(kind=dp), parameter :: r_obs = 1000.0_dp
real(kind=dp), parameter :: theta_obs = PI/2.0_dp

!Plasma density profile normalisation
!Set = 0 for vacuum
real(kind=dp), parameter :: N0 =  0.0_dp !3.50e7 !1.0e8 

!Integration parameters
character(len=20), parameter :: IntegrationType = 'Backwards' !Forwards/Backwards

!Type for integrating multiple trajectories
character(len=20), parameter :: mode = 'shoot' ! 'frequency', 'image', 'equator', 'single' , 'shoot'



!Parameters used if forward ray tracing
real(kind=dp), parameter :: r_init = 800.0_dp , theta_init = PI/2.0_dp, phi_init = 0.0_dp
real(kind=dp), parameter :: dir_theta = PI/4.0_dp , dir_phi = 0.0_dp !0.0_dp



!IO parameters
integer(kind=dp) :: Nrows = 6, Ncols = 1e6
character(len=200) :: IOpath = "/Users/tomkimpson/Data/ThesisData/RT/"





end module parameters
