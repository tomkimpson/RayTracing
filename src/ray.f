module ray


use parameters 
use constants
use IC
use numerical_methods
use IO


implicit none

private calculate_ds

public run

contains



subroutine run(alpha,beta,nu,ds,plot)
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


!call calculate_ds_cartesian(output(:,counts),ds)

endif





end subroutine run






subroutine calculate_ds(v,ds)

real(kind=dp), dimension(6) :: v !input vector
real(kind=dp) :: ds,x,y,z,r,theta,phi,m

!Load data
r = v(1) ; theta = v(2) ; phi = v(3) !+ 2.0_dp * PI
m = sqrt(r**2 + a2)

!Convert to cartesian
x = m*sin(theta)*cos(phi)
y = m*sin(theta)*sin(phi)
z = r*cos(theta)


!Calcuale the square of the difference
ds = (x - xTarget)**2 + (y - yTarget)**2 + (z-zTarget)**2

!ds = (r - rTarget)**2 + (cos(theta) - cos(thetaTarget))**2 + (sin(phi)-sin(phiTarget))**2

end subroutine calculate_ds


subroutine calculate_ds_cartesian(v,ds)


real(kind=dp), dimension(6) :: v !input vector
real(kind=dp) :: ds,x,y,z,r,theta,phi,m

!Load data
r = v(1) ; theta = v(2) ; phi = v(3) !+ 2.0_dp * PI
m = sqrt(r**2 + a2)

!Convert to cartesian
x = m*sin(theta)*cos(phi)
y = m*sin(theta)*sin(phi)
z = r*cos(theta)


!Calcuale the square of the difference
ds = (x - xTarget)**2 + (y - yTarget)**2 + (z-zTarget)**2


end subroutine calculate_ds_cartesian




end module ray




