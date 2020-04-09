program main

use parameters

implicit none
real(kind=dp) :: aL, aU, bL ,bU !image plane limits
real(kind=dp) :: nu0, nu1 !Frequency limits
real(kind=dp) :: alpha, beta, nu_obs
integer(kind=dp),parameter :: N = 100
integer(kind=dp) :: i,j




!Some housekeeeping
call setup()




if (mode .EQ. 'frequency') then
    alpha = -7.0_dp; beta = 0.0_dp !Only used if integration mode = backwards
    nu0 = 0.180_dp ; nu1 = 6.0_dp !GHz
    do i = 0,N
    print *, N-i
    nu_obs = nu0 + (nu1-nu0)/real(N,kind=dp) * real(i,kind=dp)
    call run(alpha,beta,nu_obs)
    enddo

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
    call run(alpha,beta,nu_obs)

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
    call run(alpha,beta,nu_obs)
    enddo




else if (mode .EQ. 'single') then

    alpha = -7.0_dp ; beta = 0.0_dp ; nu_obs = 0.210_dp
    call run(alpha,beta,nu_obs)


endif


end program main



subroutine run(alpha,beta,nu)

use parameters
use IC
use numerical_methods
use IO
implicit none


!Arguments



real(kind=dp),intent(in) :: alpha,beta,nu
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
    counts = counts + 1
    output(:,counts) = x
    
    !Exit condition
    if (x(1) .GT. r_obs .and. x(5) .GT. 0.0_dp) then
    code = 1
    Exit
    endif


enddo


!I/O
call ToFile(output,counts,alpha,beta,nu,c(1))






end subroutine run








subroutine setup()
use parameters
use constants
!Some useful setup and print statements

if (dp .EQ. 8) then
escal = 1.0e15
else if (dp .EQ. 16) then
escal = 1.0e19
endif

end subroutine setup



