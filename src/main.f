program main

use parameters
use constants
use IO

implicit none
real(kind=dp) :: aL, aU, bL ,bU !image plane limits
real(kind=dp) :: nu0, nu1 !Frequency limits
real(kind=dp) :: alpha, beta, nu_obs
integer(kind=dp),parameter :: N = 1000
integer(kind=dp) :: i,j
real(kind=dp) :: ds
real(kind=dp),dimension(4) :: globals

!Array to store optimisation info
real(kind=dp), dimension(N,2) :: optim

!Some housekeeeping
call setup()




if (mode .EQ. 'frequency') then
    alpha = -7.0_dp; beta = 0.0_dp !Only used if integration mode = backwards
    nu0 = 0.180_dp ; nu1 = 6.0_dp !GHz
    do i = 0,N
    print *, N-i
    nu_obs = nu0 + (nu1-nu0)/real(N,kind=dp) * real(i,kind=dp)
    call run(alpha,beta,nu_obs,ds,1)
    enddo


else if (mode .EQ. 'shoot') then


    if (IntegrationType .EQ. 'Forwards') then
    print *, "Error: Cannot do forward integration in shooting mode"
    stop
    endif


    if (N0 .ne. 0.0_dp) then
    print *, 'You are no longer ray tracing in vacuum. Check the frequency'
    stop
    endif


    !Set the target point
    !Note - will likely just read this in from file in future
    
    rTarget = 160.0_dp ; thetaTarget = PI/2.0_dp ; phiTarget = 2.760_dp
    
    xTarget = sqrt(rTarget**2 + a2) * sin(thetaTarget) * cos(phiTarget)
    yTarget = sqrt(rTarget**2 + a2) * sin(thetaTarget) * sin(phiTarget)
    zTarget = rTarget * cos(thetaTarget)
    

    !Initial guess
    alpha = -30.0_dp !yTarget !0.01_dp
    beta = zTarget !0.01_dp
   
    !Set the ray frequency
    nu_obs = 100.0_dp !Doesnt matter for vacuum
    
    !Setup for optimisation
    ds = 100.0_dp
    globals = 0.0_dp
    
    !Do the first run
    call run(alpha,beta,nu_obs,ds,1)
    print *, 'Origin:', alpha,beta,ds

   !Then optimise to find the minimum
    !do while (ds .GT. ds_eps)
    do i=1,20
    call optimise_alpha_beta(alpha,beta,nu_obs,ds,globals)
    enddo


    call run(alpha,beta,nu_obs,ds,1)
    


    print *, 'completed'
    print *, ds
!    stop

else if (mode .EQ. 'image') then

    if (IntegrationType .EQ. 'Forwards') then
    print *, "Error: Cannot do forward integration in image mode"
    stop
    endif

    aL = -6.0_dp ; aU = 6.0_dp ; bL = -6.0_dp ;bU = 6.0_dp
    nu_obs = 100.0_dp !GHz
    
    do i = 0,N
    alpha = aL + (aU-aL)/real(N,kind=dp) * real(i,kind=dp)
    print *, i
    do j = 0,N
    beta = bL + (bU-bL)/real(N,kind=dp) * real(j,kind=dp)
    call run(alpha,beta,nu_obs,ds,1)

    enddo
    enddo



else if (mode .EQ. 'equator') then


    if (IntegrationType .EQ. 'Forwards') then
    print *, "Error: Cannot do forward integration in image mode"
    stop
    endif

    aL = -6.0_dp ; aU = 6.0_dp ;
    nu_obs = 100.0_dp !GHz
    
    do i = 0,N
    alpha = aL + (aU-aL)/real(N,kind=dp) * real(i,kind=dp)
    beta = 0.0_dp
    call run(alpha,beta,nu_obs,ds,1)
    enddo




else if (mode .EQ. 'single') then

    alpha = -7.0_dp ; beta = 0.0_dp ; nu_obs = 0.210_dp
    call run(alpha,beta,nu_obs,ds,1)


endif


end program main



subroutine run(alpha,beta,nu,ds,plot)

use parameters
use IC
use numerical_methods
use IO
implicit none


!Arguments



real(kind=dp),intent(in) :: alpha,beta,nu
real(kind=dp),intent(out) :: ds
integer(kind=4), intent(in) :: plot !plot overrride. if = 1, always plot
!Other
real(kind=dp), dimension(6) :: x 
real(kind=dp), dimension(4) :: c 
integer(kind=dp) :: counts
real(kind=dp), dimension(Nrows, Ncols) :: output !Array to save to. Column major
integer(kind=dp) :: code

!Set the intial conditions
call initial_conditions(alpha,beta,nu,x,c)

code = 0
counts = 1
output(:,counts) = x

!Integrate numerically
do while (x(1) .GT. Rhor)

    call RKF(x,c)


    if (counts .GT. Ncols) then
    print *, 'Error: counts > Ncols', counts, Ncols
    stop
    endif


    counts = counts + 1
    output(:,counts) = x
 
    

   
    !Exit condition. Intersection search
    if (c(3) .LT. 0.0_dp) then
 
    call calculate_ds(x,ds)


    code = 2

    exit
    endif

    !Exit condition. !General ray tracing
    if (x(1) .GT. r_obs .and. x(5) .GT. 0.0_dp) then
    code = 1
    Exit
    endif


enddo


!I/O
if (plot .eq. 1) then
call ToFile(output,counts,alpha,beta,nu,c(1))
endif





end subroutine run



subroutine optimise_alpha_beta(alpha,beta,nu_obs,ds,globals)

use parameters
use constants

implicit none

!Arguments
real(kind=dp), intent(inout) :: alpha,beta,ds
real(kind=dp),intent(in) :: nu_obs
real(kind=dp), intent(inout),dimension(4) :: globals

!Other
real(kind=dp) :: ds_alpha, ds_beta
real(kind=dp) :: gA, gB
real(kind=dp) :: zeta,hA,hB
real(kind=dp) :: eta, t, ds_trial
real(kind=dp) :: ds1, alpha1, beta1
real(kind=dp) :: trial_alpha, trial_beta


!Get the gradients in the alpha/beta directions

!Gradient alpha
call run(alpha+dg,beta,nu_obs,ds_alpha,0)
gA = - ((ds_alpha - ds)/dg)


!Gradient beta
call run(alpha,beta+dg,nu_obs,ds_beta,0)
gB = - ((ds_beta - ds)/dg)


if (globals(1) .EQ. 0.0_dp) then
!First go
zeta = 0.0_dp
else
zeta = (gA*gA +gB*gB)/(globals(1)**2 + globals(2)**2)
endif

hA = gA +zeta*globals(3)
hB = gB +zeta*globals(4)


!Normalise the direction. Just tells you higher or lower info now
hB = gB / abs(gA)
hA = gA / abs(gA)




!Got the direction, now get the stepsizze by performing a line search
eta = 2.0_dp !0.50_dp
t = global_t !5e-9
ds1 = 1e20


11 do 

    !Trial values
    trial_alpha = alpha + t*hA
    trial_beta = beta + t*hB

    !Try a step in the direction
    call run(trial_alpha,trial_beta,nu_obs,ds_trial,0)

    print *, 'Tr:', trial_alpha,ds_trial,t

    if (ds_trial .LT. ds) then

    !Update the best values
    alpha1 = alpha+t*hA
    beta1 = beta+t*hB
    ds1 = ds_trial
    global_t = t 



    exit


    else
    !Has stopped improving.
    !Use this stepsize going forwards
    
    !exit
    if (t .LT. 1e-7) then
    print *, 'no iszeable change', ds1,ds
    exit    
    endif


    t = t/eta
    goto 11
    endif






enddo




!Update before exiting subroutine 

if (ds1 .LT. ds) then
globals(1) = gA ; globals(2) = gB ; globals(3) = hA ; globals(4) = hB
alpha = alpha1 ; beta = beta1 ; ds = ds1

else
!Reset, gives it a kick
print *, 'Failure- needs reset'
print *, 'Gradients', gA, gB
print *, 'attempt = ', alpha,t*hA
globals = 0.0_dp
!stop
!global_t = 1e-3 !1.0_dp !global_t/10.0_dp
endif



!print *, 'out = ', alpha,beta,ds
print *, '---------------------'
!stop

end subroutine optimise_alpha_beta



subroutine calculate_ds(v,ds)

use parameters
use constants

implicit none

real(kind=dp), dimension(6) :: v !input vector
real(kind=dp) :: ds,x,y,z,r,theta,phi,m

!Load data
r = v(1) ; theta = v(2) ; phi = v(3) + 2.0_dp * PI
m = sqrt(r**2 + a2)

!Convert to cartesian
x = m*sin(theta)*cos(phi)
y = m*sin(theta)*sin(phi)
z = r*cos(theta)



!Calcuale the square of the difference


ds = (x - xTarget)**2 + (y - yTarget)**2 + (z-zTarget)**2


print *, 'cartesian norm:', (x - xTarget)**2 + (y - yTarget)**2 + (z-zTarget)**2


!ds = (r - rTarget)**2 + (theta - thetaTarget)**2 + (phi-phiTarget)**2



end subroutine calculate_ds



subroutine setup()
use parameters
use constants
!Some useful setup and print statements

if (dp .EQ. 8) then
escal = 1.0e15
dg = 1e-5
ds_eps = 1e-1
else if (dp .EQ. 16) then
escal = 1.0e19
dg = 1.0e-14 !Numerical gradient
ds_eps = (1.0e-6)**2
endif

end subroutine setup



