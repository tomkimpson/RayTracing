module constants

use parameters

implicit none




!Universal constants

real(kind=dp), parameter :: Newton_g = 6.67408d-11 
real(kind=dp), parameter :: Msolar = 1.989d30 
real(kind=dp), parameter :: mu = Newton_g*MBH*Msolar
real(kind=dp), parameter :: light_c = 3.0d8
real(kind=dp), parameter :: convert_m = light_c**2/(Newton_g*MBH*Msolar) !Multiply by this to go TO Natural units
real(kind=dp), parameter :: convert_s = light_c**3/(Newton_g*MBH*Msolar) !Multiply by this to go TO Natural units
real(kind=dp), parameter :: convert_spin= light_c/(Newton_g*(MBH*Msolar)**2.0_dp) !Multiply by this to go TO Natural units
real(kind=dp), PARAMETER :: electron_charge = 4.803204250e-10 !CGS
real(kind=dp), PARAMETER :: electron_mass = 9.109383560e-28 !CGS


!Useful constants
real(kind=dp), parameter :: Rhor = 1.0_dp + sqrt(1.0_dp - a**2) + 1.0d-2 !Horizon + eps
real(kind=dp), parameter :: a2 = a**2


!Target points for intersection search
real(kind=dp) :: xTarget, yTarget, zTarget
real(kind=dp) :: rTarget,thetaTarget, phiTarget

!Set the 'x precision' for intersection search
real(kind=dp),parameter :: dx_eps = epsilon(Rhor)

!Set the target intersection precision
real(kind=dp) :: ds_eps 

!Set the gradient bit
real(kind=dp) :: dg 

!Used in pattern search
real(kind=dp) :: decay_factor



!Temporay
real(kind=dp) :: global_t = 10.0_dp



!Integration constants
real(kind=dp) :: escal !This will be made higher / lower for double vs quad precision
real(kind=dp), parameter :: S = 0.90_dp
real(kind=dp), parameter :: Pgrow = -0.20_dp
real(kind=dp), parameter :: Pshrink = -0.250_dp
real(kind=dp), parameter :: errcon = (5.0_dp/S)**(1.0_dp/Pgrow)

!Related to plasma



!Cash-Karp parameters
real(kind = dp) :: B21=1.0_dp/5.0_dp
real(kind = dp) :: B31 = 3.0_dp/40.0_dp , B32 = 9.0_dp/40.0_dp
real(kind = dp) :: B41 = 3.0_dp/10.0_dp, B42 = -9.0_dp/10.0_dp, B43 = 6.0_dp/5.0_dp 
real(kind = dp) :: B51 = -11.0_dp/54.0_dp, B52 = 5.0_dp/2.0_dp, B53 = -70.0_dp/27.0_dp
real(kind = dp) :: B54 = 35.0_dp/27.0_dp
real(kind = dp) :: B61 = 1631.0_dp/55296.0_dp, B62 = 175.0_dp/512.0_dp, B63 = 575.0_dp/13824.0_dp
real(kind = dp) :: B64 = 44275.0_dp/110592.0_dp, B65 = 253.0_dp/4096.0_dp
real(kind = dp) :: c1 = 37.0_dp/378.0_dp, c3 = 250.0_dp/621.0_dp, c4 = 125.0_dp/594.0_dp
real(kind = dp) :: c6=512.0_dp/1771.0_dp
real(kind = dp) :: cbar1 = 2825.0_dp/27648.0_dp, cbar3 = 18575.0_dp/48384.0_dp
real(kind = dp) :: cbar4=13525.0_dp/55296.0_dp, cbar5 = 277.0_dp/14336.0_dp, cbar6 = 1.0_dp/4.0_dp





end module constants
