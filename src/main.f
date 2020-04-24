program main

use parameters
use constants
!use IO
!use ray
!use optimization

implicit none
integer(kind=dp) :: i


!Some housekeeeping
call setup()

!Run the program
if (mode .EQ. 'frequency') then
    call subroutine_frequency()
    
else if (mode .EQ. 'shoot') then
    call subroutine_shoot()

else if (mode .EQ. 'image') then
    call subroutine_image()
    

else if (mode .EQ. 'equator') then
    call subroutine_equator()
    
else if (mode .EQ. 'single') then

    call subroutine_single()

endif



print *, 'Program Completed'
end program main


subroutine subroutine_equator()


use parameters
use constants
use ray

implicit none

real(kind=dp) :: alpha, beta,nu_obs
real(kind=dp) :: aL, aU, bL ,bU !image plane limits
integer(kind=dp) :: i,j
integer(kind=dp),parameter :: N = 100
real(kind=dp) :: ds


if (IntegrationType .EQ. 'Forwards') then
print *, "Error: Cannot do forward integration in equator mode"
stop
endif



aL = -6.0_dp ; aU = 6.0_dp ;
nu_obs = 100.0_dp !GHz
beta = 0.0_dp


do i = 0,N
    alpha = aL + (aU-aL)/real(N,kind=dp) * real(i,kind=dp)
    beta = 0.0_dp
    call run(alpha,beta,nu_obs,ds,1)
enddo



end subroutine subroutine_equator





subroutine subroutine_image()


use parameters
use constants
use ray

implicit none

real(kind=dp) :: alpha, beta,nu_obs
real(kind=dp) :: aL, aU, bL ,bU !image plane limits
integer(kind=dp) :: i,j
integer(kind=dp),parameter :: N = 100
real(kind=dp) :: ds


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



end subroutine subroutine_image




subroutine subroutine_single()

use parameters
use constants
use ray

implicit none
real(kind=dp) :: alpha, beta, nu_obs
real(kind=dp) :: ds

alpha = -7.0_dp ; beta = 0.0_dp ; nu_obs = 0.210_dp
call run(alpha,beta,nu_obs,ds,1)

end subroutine subroutine_single




subroutine subroutine_frequency()


use parameters
use constants
use ray

implicit none


real(kind=dp) :: alpha, beta, nu_obs
integer(kind=dp),parameter :: N = 100
real(kind=dp) :: nu0, nu1 !Frequency limits
real(kind=dp) :: ds
integer(kind=dp) :: i


alpha = -7.0_dp; beta = 0.0_dp !Only used if integration mode = backwards
nu0 = 0.180_dp ; nu1 = 6.0_dp !GHz



do i = 0,N
    print *, N-i
    nu_obs = nu0 + (nu1-nu0)/real(N,kind=dp) * real(i,kind=dp)
    call run(alpha,beta,nu_obs,ds,1)
enddo



end subroutine subroutine_frequency










subroutine subroutine_shoot()

use parameters
use constants
implicit none

!For iterating
integer(kind=dp) :: i

!For reading
integer(kind=dp) :: IOstatus
real(kind=dp), dimension(11) :: read_row
real(kind=dp) :: nu_obs

    if (IntegrationType .EQ. 'Forwards') then
    print *, "Error: Cannot do forward integration in shooting mode"
    stop
    endif



    if (N0 .ne. 0.0_dp) then
    print *, 'You are no longer ray tracing in vacuum. Be sure to check the frequency'
    endif



    if (load .eq. 1) then
        !Load the target points from file
        IOstatus = 0
        open(unit=15, file = targets_file,status='old',access='sequential',form='formatted',action='read')


        do
            read(15,*,IOSTAT=IOstatus) read_row
        
            if (IOstatus .EQ. 0)then
                rTarget = read_row(9); thetaTarget = read_row(10); phiTarget = read_row(11)
                uvector(1:4) = read_row(5:8)
                print *, 'U in=', uvector
                print *, ' polar Target = ', rTarget,thetaTarget,phiTarget            
                
                nu_obs = 1.0_dp
                call find_intersecting_rays(nu_obs,1)
            else
                print *, 'reached end of file'
                close(15)
                exit
            endif


        
        enddo
        
        close(15)


    else
        !Just specify by hand    
        rTarget = 160.0_dp ; thetaTarget = PI/2.0_dp ; phiTarget = 2.80_dp 
        nu_obs = 0.10_dp
        
        do i = 1,10
            nu_obs = nu_obs + 1.0_dp
            call find_intersecting_rays(nu_obs,1)
        enddo
   endif



end subroutine subroutine_shoot





subroutine find_intersecting_rays(nu_obs,ray_order)

use parameters
use constants
use ray
use optimization

implicit none

!priamry/secondary ray
integer(kind=4), intent(in) :: ray_order

!Imsge plane coords
real(kind=dp) :: alpha,beta
!f(alpha,beta) = ds
real(kind=dp) :: ds
!Frequency
real(kind=dp),intent(in) :: nu_obs

!Integers
integer(kind=dp) :: i,stat


!Globals - used for high precision gradient descent
real(kind=dp),dimension(4) :: globals


!Convert it to Cartesian
xTarget = sqrt(rTarget**2 + a2) * sin(thetaTarget) * cos(phiTarget)
yTarget = sqrt(rTarget**2 + a2) * sin(thetaTarget) * sin(phiTarget)
zTarget = rTarget * cos(thetaTarget)
print *, 'Cartesian Targets = ', xTarget, yTarget, zTarget



!Initial guess - primary ray

!This should really be a subroutine

if (N0 .ne. 0.0_dp) then
!If youre not working in vacuum assumes that you are just looking at a single target point
!Useful as can then use previous info for initial alpha,beta guess


if (alpha_previous .ne. 0.0_dp) then

alpha = alpha_previous ; beta = beta_previous

else

if (ray_order .eq. 1) then
    alpha = yTarget
elseif (ray_order .eq. 2) then
    alpha = -yTarget
endif
beta = zTarget


endif





endif




print *, alpha_previous, beta_previous
print *, alpha,beta




stat = 0
!Do the first run
call run(alpha,beta,nu_obs,ds,0)



!Setup for optimisation
if (optimizer .eq. 'CGD') then
globals = 0.0_dp
global_t = 5.0d-9 

print *, 'Begig optimization for ray order =', ray_order, ' and optimizer= ', optimizer


    do while (ds .GT. ds_eps)
    call optimise_alpha_beta(alpha,beta,nu_obs,ds,stat,globals)
 

    if (stat .ne. 0) then
    exit
    endif


    enddo


elseif (optimizer .eq. 'PS') then
dg = 2.0_dp
decay_factor = 2.0_dp

    do while (ds .GT. ds_eps)
    call pattern_search(alpha,beta,nu_obs,ds,stat)
    


    
    if (stat .ne. 0) then
    exit
    endif



    enddo

else
print *, 'Optimiser type:', optimizer, ' not recognized'
stop

endif




!Do a run with the found parameters
if (stat .eq. 0) then
print *, 'Optimisation converged successfully for ray ',ray_order , ' with alpha/beta/ds = ',alpha,beta, ds
call run(alpha,beta,nu_obs,ds,1)

alpha_previous = alpha
beta_previous = beta

endif
    

end subroutine find_intersecting_rays












subroutine setup()

use parameters
use constants
!Some useful setup and print statements

!Precision specific 

if (dp .EQ. 8) then
escal = 1.0e15
dg = 2.0_dp
ds_eps = 1e-6
dx_eps = 1.0e-12 
else if (dp .EQ. 16) then
escal = 1.0e19
bit = 1.0e-16 !Numerical gradient
dx_eps = 1.0d-19 
ds_eps = (1.0e-6)**2
!ds_eps = 1e-2
endif




end subroutine setup



