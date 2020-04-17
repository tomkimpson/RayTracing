module numerical_methods

use parameters
use constants
use odes

implicit none

private GrowStepsize, ShrinkStepsize


public RKF

contains


subroutine RKF(y,c)
!Integrate using the Runge-Kutta-Fehlberg algorithm

!Arguments
real(kind=dp), dimension(6), intent(inout) :: y
real(kind=dp), dimension(4) ,intent(inout) :: c

!Other
real(kind=dp), dimension(size(y)) :: y1,y2,y3,y4,y5,y6
real(kind=dp), dimension(size(y)) :: k1,k2,k3,k4,k5,k6
real(kind=dp), dimension(size(y)) :: dy1, dy2, dy3, dy4, dy5,dy6
real(kind=dp), dimension(size(y)) :: ynew, yerr
real(kind=dp), dimension(size(y)) :: deltaErr, yscal, ratio
real(kind=dp) :: errmax, h

!Some reals used for intersection search

!The x cooordinate at this timestep
real(kind=dp) :: x
! the difference between the x coodinate and the xtarget
real(kind=dp) :: dx
!shows we are now iterating close to intersection points
integer(kind=dp) :: intersecting


intersecting = 0


11 continue


h = c(3)


! Y1
y1 = y
call derivs(y1,c,dy1)
k1 = h * dy1



!Y2
y2 = y1 + B21*k1
call derivs(y2,c,dy2)
k2 = h * dy2



!Y3
y3 = y1 + B31*k1 + B32*k2
call derivs(y3,c,dy3)
k3 = h * dy3


!Y4
y4 = y1 + B41*k1 + B42*k2 + B43*k3
call derivs(y4,c,dy4)
k4 = h * dy4


!Y5
y5 = y1 + B51*k1 + B52*k2 + B53*k3 + B54*k4 
call derivs(y5,c,dy5)
k5 = h * dy5


!Y6
y6 = y1 + B61*k1 + B62*k2 + B63*k3 + B64*k4 + B65*k5
call derivs(y6,c,dy6)
k6 = h * dy6


!Update
ynew = y1 + c1*k1  + c3*k3 + c4*k4  +c6*k6 
yerr = y1 + cbar1*k1 + cbar3*k3 + cbar4*k4 + cbar5*k5 + cbar6*k6

deltaErr = abs(ynew - yerr)
yscal = abs(y1) + abs(k1) + 1.0d-3
ratio = deltaErr/yscal
errmax = escal * maxval(ratio)




!Some conditions that are only used in shooting mode




if (mode .EQ. 'shoot') then



    !Search for overstep c.f. target point
    !This is currently specialised to the observer in the equatoiral plane.
    x = sqrt(ynew(1)**2 + a2) * sin(ynew(2))*cos(ynew(3))
    dx = x - xTarget

    
    if (y(1) .eq. ynew(1)) then
    !stepsize is so small that variables are no longer updating    
    !Take what youve got and exit
    y = ynew
    c(3) = -1.0_dp !Signifies to outer routine to quit
    return
    endif



    if (abs(dx) .LT. 1e-15) then
    y = ynew
    c(3) = -1.0_dp !Signifies to outer routine to quit
    return
    endif


    if (abs(dx) .LT. dx_eps) then
    !Update and exit
    y = ynew
    c(3) = -1.0_dp !Signifies to outer routine to quit
    return
    endif

    if (dx .LT. 0.0_dp) then
    !Overstepped
    !Step again with smaller h
    !set intersecting = 1 to signify to the rest of the program that we are near the intersection point
    !This switches us to a fixed stepsize which is only updated below
    intersecting = 1
    c(3) = c(3) / 2.0_dp
    goto 11
    endif



endif






if (intersecting .EQ. 0) then

!Just update adaptively as usual

    if (errmax .GT. 1.0_dp) then
    !This is not good. Do not update yOUT and reduce the stepsize
    call ShrinkStepsize(errmax,h)
    c(3) = h
    goto 11
    else
    !This is good. Update yOUT and try to increase the stepsize a little bit
    call GrowStepsize(errmax,h)
    c(3) = h
    y = ynew
    endif


else 
!switch to a fixed stepsize
    
    y = ynew
    c(3) = h
    goto 11
endif




end subroutine RKF



subroutine GrowStepsize(errmax,h)
real(kind=dp) :: errmax,h
if (errmax .GT. errcon) then
h = S*h*errmax**Pgrow
else
h = h * 5.0_dp
endif
end subroutine GrowStepsize



subroutine ShrinkStepsize(errmax,h)
real(kind=dp) :: errmax, h
real(kind=dp) :: htemp

htemp = S*h*errmax**Pshrink
h = sign(max(abs(htemp),0.10_dp*abs(h)),h)

end subroutine ShrinkStepsize




end module numerical_methods
