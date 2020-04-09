module odes

use parameters
use constants

implicit none

private 


public derivs,density_profile

contains


subroutine derivs(y,c, dy)

!Arguments
real(kind=dp), intent(IN), dimension(6) :: y
real(kind=dp), intent(OUT), dimension(6) :: dy !deriatives
real(kind=dp), intent(in), dimension(4) :: c

!Other
real(kind=dp) :: r, theta,pr,ptheta
real(kind=dp) :: L, kappa
real(kind=dp) :: delta, sigma,r2,ctheta, stheta, SD, B2
real(kind=dp) :: f,g,f1,g1
real(kind=dp) :: plasma_pr, plasma_ptheta


!Load the variables and constants
r = y(1) ; theta = y(2) ; pr = y(5) ; ptheta=y(6)
L = c(1) ; kappa = c(2) ; B2 = c(4)

!Define some useful stuff
r2 = r**2
ctheta = cos(theta)
stheta = sin(theta)


delta = r2 - 2.0_dp*r + a2
sigma = r2 -a2*ctheta
SD = sigma*delta

!Set of ODEs to solve
dy(1) = pr*delta/sigma
dy(2) = ptheta / sigma
dy(3) = (2.0_dp*a*r + (sigma-2.0_dp*r)*L/stheta**2) / (SD)
dy(4) = 1.0_dp + (2.0_dp*r*(r2 + a2)-2.0_dp*a*r*L) /SD


call density_profile(r,theta,f,g,f1,g1)

plasma_pr = -B2*(f1*delta/2.0_dp + (r-1.0_dp)*f) / SD
plasma_ptheta = -B2*g1/(2.0_dp*sigma)

dy(5) = (-kappa*(r-1.0_dp) +2.0_dp*r*(r2 + a2) - 2.0_dp*a*L)/SD - 2.0_dp*pr**2*(r-1.0_dp)/sigma + plasma_pr
dy(6) = ctheta*stheta*(L**2/stheta**4 - a2) / sigma + plasma_ptheta



end subroutine derivs


subroutine density_profile(r,theta,f,g,f1,g1)
!Arguments
real(kind=dp), intent(in) :: r,theta
real(kind=dp), intent(out) :: f,g,f1,g1

f = r**(0.90_dp)
g = 0.0_dp

f1 = 0.90_dp*r**(-0.10_dp)
g1 = 0.0_dp


end subroutine density_profile

end module odes
