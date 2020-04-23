
module IC


use parameters 
use constants
use odes

implicit none

private backwards

public initial_conditions

contains



subroutine initial_conditions(alpha,beta,nu_obs,x,c)
!Arguments
!image plane coordinates
real(kind=dp),intent(in) :: alpha,beta,nu_obs 
!coordinate vector 
real(kind=dp), intent(out), dimension(6) :: x 
!constants vector 
real(kind=dp), intent(out), dimension(4) :: c 

!Other

!Position/Velocity vector
real(kind=dp), dimension(6) :: PV_vector
real(kind=dp) :: r,theta,phi,rdot0, thetadot0,phidot0
!Time coordinate
real(kind=dp) :: t
!momenta
real(kind=dp) :: sig,del, pr, ptheta 
!energy
real(kind=dp) :: s1,E2,En, Eprime, Eobs
!Angular momentum
real(kind=dp) :: Lz 
!Constant kappa
real(kind=dp) :: kappa
!Plasma frequency, number density terms (f,g) and derivatives (f1,g1)
real(kind=dp) :: omega2, f,g,f1,g1
real(kind=dp) :: B2 


if (IntegrationType .EQ. 'Backwards') then

    call backwards(alpha, beta, PV_vector)

else if (IntegrationType .EQ. 'Forwards') then

    call forwards(PV_vector)

endif

!Etract info from vector
r = PV_vector(1) ; theta = PV_vector(2) ; phi = PV_vector(3)
rdot0 = PV_vector(4) ; thetadot0 = PV_vector(5) ; phidot0 = PV_vector(6)


!Define some useful quantitites
sig = r**2.0_dp +(a*cos(theta))**2.0_dp
del = r**2.0_dp - 2.0_dp*r +a**2



!Time coord
t=0.0_dp



!Plasma
B2= 4.0_dp*PI*electron_charge**2 *N0 / electron_mass
call density_profile(r,theta,f,g,f1,g1)
omega2 = B2 *(f+g) / sig


!Calculate ray velocity magnitude
s1 = sig-2.0*r
Eobs = 2.0_dp*PI*nu_obs*1e9
Eobs = 1.0_dp
Eprime = sqrt( (Eobs**2 - omega2 + 2.0_dp*r*omega2/sig) / (s1*(rdot0**2 /del + thetadot0**2) + del*sin(theta)**2*phidot0**2) )

!... and correct velocity components
rdot0 = rdot0 * Eprime ; thetadot0 = thetadot0 * Eprime ; phidot0 = phidot0 * Eprime

!...and double check the derived energy matches the one you want
E2 = s1*(rdot0**2.0/del +thetadot0**2.0 +omega2/sig) + del*sin(theta)**2.0*phidot0**2.0
En  = sqrt(E2)
!print *, En
!print *, Eobs


!Momenta
pr = rdot0 *sig/del
ptheta = sig*thetadot0



!Angular momentum
Lz = ((sig*del*phidot0 - 2.0*a*r*En) * sin(theta)**2.0/s1)

!Normalize to E=1
pr = pr/Eobs
ptheta = ptheta/Eobs
Lz = Lz/Eobs
B2 = B2 / Eobs**2

!One more constant to define
kappa = ptheta**2.0 + a**2*sin(theta)**2.0 + Lz**2.0/sin(theta)**2.0

!Save to vectors for output
x(1) = r
x(2) = theta
x(3) = phi
x(4) = t
x(5) = pr
x(6) = ptheta

c(1) = Lz
c(2) = kappa
c(3) = 1.0e-6 !initial stepsize
c(4) = B2



!Write stuff to global
Eobs_global = Eobs



end subroutine initial_conditions


subroutine backwards(alpha,beta,PV)
!Arguments
real(kind=dp), intent(in) :: alpha,beta
real(kind=dp), intent(out),dimension(6) :: PV
!Other

!Cartesian coords
real(kind=dp) :: xprime, yprime, zprime 

!BL coords
real(kind=dp) :: w,r,theta,phi


!derivative
real(kind=dp) :: sig,u,zdot,v,rdot0,thetadot0,phidot0,del

!Covert to primed Cartesian frame
xprime = sqrt(r_obs**2.0_dp +a**2) * sin(theta_obs) - beta*cos(theta_obs)
yprime = alpha
zprime = r_obs*cos(theta_obs) + beta*sin(theta_obs)


! Boyer-Lindquist
w = xprime**2.0_dp +yprime**2.0_dp +zprime **2.0_dp - a**2
r = sqrt((w+sqrt(w**2.0+4.0*a**2*zprime**2.0_dp))/2.0_dp)
theta = acos(zprime/r)
phi = atan2(yprime,xprime)


!Derivatives
sig = r**2.0_dp +(a*cos(theta))**2.0_dp
del = r**2.0_dp - 2.0_dp*r +a**2
u = sqrt(r**2.0_dp+a**2)
zdot = -1.0_dp
v= -sin(theta_obs)*cos(phi)

rdot0 = -zdot*(-u**2.0*cos(theta_obs)*cos(theta)+r*u*v*sin(theta))/sig
thetadot0 = -zdot*(cos(theta_obs)*r*sin(theta)+u*v*cos(theta))/sig
phidot0 = -zdot*sin(theta_obs)*sin(phi)/(u*sin(theta))



!write for output
PV(1) = r ;PV(2) = theta ;PV(3) = phi
PV(4) = rdot0 ;PV(5) = thetadot0 ;PV(6) = phidot0


end subroutine backwards


subroutine forwards(PV)
!Arguments
real(kind=dp), intent(out),dimension(6) :: PV

!Other
real(kind=dp) :: r,theta,phi
real(kind=dp) :: xdot, ydot, zdot
real(kind=dp) :: rdot0, thetadot0, phidot0
real(kind=dp) :: mm, sigma

!Load initial position from parameters file
r= r_init ; theta = theta_init ; phi = phi_init

!Cartesian direction of ray from angles in parameters.f
xdot = sin(dir_theta)*cos(dir_phi)
ydot = sin(dir_theta)*sin(dir_phi)
zdot = cos(dir_theta)


!Convert to BL
sigma = r**2.0_dp +(a*cos(theta))**2.0_dp
mm = sqrt(r**2 + a**2)
rdot0 = mm*r*sin(theta)*cos(phi)*xdot/sigma &
       +mm*r*sin(theta)*sin(phi)*ydot/sigma &
       +mm**2*cos(theta)*zdot/sigma



thetadot0 = (mm*cos(theta)*cos(phi) * xdot &
           +mm*cos(theta)*sin(phi) * ydot &
           -r*sin(theta)* zdot&
           )/sigma


phidot0 = (-sin(phi)*xdot + cos(phi)*ydot)/(mm*sin(theta))


print *, xdot,ydot,zdot
print *, r,theta,phi
print *, rdot0,thetadot0,phidot0



!write for output
PV(1) = r ;PV(2) = theta ;PV(3) = phi
PV(4) = rdot0 ;PV(5) = thetadot0 ;PV(6) = phidot0

end subroutine forwards


end module IC




