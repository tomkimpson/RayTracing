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


    if (IntegrationType .EQ. 'Forwards') then
    print *, "Error: Cannot do forward integration in shooting mode"
    stop
    endif



    if (N0 .ne. 0.0_dp) then
    print *, 'You are no longer ray tracing in vacuum. Check the frequency'
    stop
    endif



    if (load .eq. 1) then
        !Load the target points from file
        IOstatus = 0
        open(unit=15, file = targets_file,status='old',access='sequential',form='formatted',action='read')


        do
            read(15,*,IOSTAT=IOstatus) read_row
        
            if (IOstatus .EQ. 0)then
                rTarget = read_row(9); thetaTarget = read_row(10); phiTarget = read_row(11)
 !               rTarget = 219.96460677885887    ;thetaTarget= 1.5707963267880871 ; phiTarget=  3.1329373667919755 
            !    rTarget = 197.60074342080654; thetaTarget=1.5707963269538816 ; phiTarget= 4.8012548074405341
                uvector(1:4) = read_row(5:8)
                print *, 'U in=', uvector
                print *, ' polar Target = ', rTarget,thetaTarget,phiTarget            
                call find_intersecting_rays()
            else
                print *, 'reached end of file'
                close(15)
                exit
            endif


        
        enddo
        
        close(15)


    else
        !Just specify by hand    
        rTarget = 160.0_dp ; thetaTarget = PI/2.0_dp ; phiTarget = 1.90_dp 
        call find_intersecting_rays()
    endif



end subroutine subroutine_shoot





subroutine find_intersecting_rays()

use parameters
use constants
use ray
use optimization

implicit none

!Imsge plane coords
real(kind=dp) :: alpha,beta
!f(alpha,beta) = ds
real(kind=dp) :: ds
!Frequency - doesnt matter here but needs to be set
real(kind=dp) :: nu_obs
!Integer for iterating if needed
integer(kind=dp) :: i,stat


!temporary for debugginh
real(kind=dp) :: dl

!Convert it to Cartesian
xTarget = sqrt(rTarget**2 + a2) * sin(thetaTarget) * cos(phiTarget)
yTarget = sqrt(rTarget**2 + a2) * sin(thetaTarget) * sin(phiTarget)
zTarget = rTarget * cos(thetaTarget)

print *, 'Cartesian Targets = ', xTarget, yTarget, zTarget



!Initial guess - primary ray
alpha = yTarget
beta = zTarget


!Set the ray frequency
nu_obs = 100.0_dp !Doesnt matter for vacuum

!Setup for optimisation
ds = 100.0_dp
dg = 2.0_dp
decay_factor = 2.0_dp


stat = 0
!Do the first run
call run(alpha,beta,nu_obs,ds,0)



!Then optimise to find the minimum
do while (ds .GT. ds_eps)
    call pattern_search(alpha,beta,nu_obs,ds,stat)
    
    if (stat .ne. 0) then
    exit
    endif

enddo


!Do a run with the found parameters
    if (stat .eq. 0) then
        print *, 'Optimisation converged successfully for primary ray with alpha/beta/ds = ',alpha,beta, ds
        call run(alpha,beta,nu_obs,ds,1)
    

    elseif (stat .eq. 5) then
    print *, '------------------------Fell into BH-----------------------------'
    print *, '------------------------Fell into BH-----------------------------'
    print *, '------------------------Fell into BH-----------------------------'
    print *, 'targets =', rTarget, thetaTarget, phiTarget
    call run(alpha,beta,nu_obs,ds,3)

    else

        print *, 'Optimisation reached a precision limit for the primary ray with alpha/beta/ds= ',alpha,beta,ds
        call run(alpha,beta,nu_obs,ds,10)

    endif




if (secondary_rays .EQ. 1 .and. xTarget .LT. 0.0_dp) then
    !Reset and search for a secondary ray    
    print *, 'Searching for secondary ray'

    alpha = -yTarget
    beta = zTarget
 
    dg = 2.0_dp
    ds = 100.0_dp    
    decay_factor = 2.0_dp

    stat = 0
    call run(alpha,beta,nu_obs,ds,0)


    ds = 1e20
    do while (ds .GT. ds_eps)
    call pattern_search(alpha,beta,nu_obs,ds,stat)

    if (stat .ne. 0) then
    !AHs reached precision limit
    exit
    endif


    enddo

    if (stat .eq. 0) then
        print *, 'Optimisation converged successfully for secondary ray with alpha/beta/ds = ',alpha,beta, ds
        call run(alpha,beta,nu_obs,ds,2)
 


    elseif (stat .eq. 5) then

    print *, 'Fell into BH'
    
    print *, '------------------------Fell into BH-----------------------------'
    print *, '------------------------Fell into BH-----------------------------'
    print *, '------------------------Fell into BH-----------------------------'
    print *, 'targets =', rTarget, thetaTarget, phiTarget


    call run(alpha,beta,nu_obs,ds,3)

   else
        print *, 'Optimisation reached a precision limit for the secondary ray with alpha/beta/ds= ',alpha,beta,ds
        call run(alpha,beta,nu_obs,ds,20)

    endif

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
!ds_eps = 1e-13
ds_eps = 1e-6
else if (dp .EQ. 16) then
escal = 1.0e19
dg = 1.0e-14 !Numerical gradient
ds_eps = (1.0e-6)**2
endif




end subroutine setup



