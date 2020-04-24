module IO

use parameters
use constants
implicit none

private calculate_covariant_metric, calculate_contravariant_metric, zamo
public ToFile, ToFile_optimisation

contains

subroutine ToFile(array,n,alpha,beta,nu_obs,c, plot_status)
real(kind=dp), dimension(Nrows, Ncols),intent(in) :: array
real(kind=dp), dimension(4),intent(in) :: c
integer(kind=dp), intent(in) :: n
integer(kind=4) :: plot_status
real(kind=dp), intent(in) :: alpha,beta,nu_obs
character(len=200) :: Fname, aStr, bStr, ID, nuStr
integer(kind=dp) :: i
real(kind=dp) :: xC, yC,zC, mm,Lz
real(kind=dp) :: r,theta, pr,ptheta,phi,t
real(kind=dp) :: nu, g
real(kind=dp),dimension(4) :: pvector,zamo_vector
real(kind=dp),dimension(4,4) :: metric

!Define savefile
print *, 'IO=',alpha
write(aStr,'(F16.10)') alpha
write(bStr,'(F16.10)') beta
write(nuStr,'(F10.2)') nu_obs
ID = 'RT_alpha='//trim(adjustl(aStr))//'_beta='//trim(adjustl(bstr))//'_nu='//trim(adjustl(nuStr))//'.txt'
Fname = trim(adjustl(IOpath)) // trim(adjustl(ID))






!Write
open(unit=10, file = Fname, status = 'replace',form = 'formatted')

do i=1,n

r = array(1,i) ; theta = array(2,i); phi = array(3,i) ; t = array(4,i)
pr = array(5,i) ; ptheta = array(6,i)


!Define the general momentum vector - p0 has been normalised, so here we 'unnormalise'

pvector(1) = -Eobs_global
pvector(2) = pr * Eobs_global
pvector(3) = ptheta * Eobs_global
pvector(4) = c(1) * Eobs_global



if (i .eq. n) then



!call zamo(r,theta,zamo_vector)



call calculate_covariant_metric(r,theta,metric)



!print *, '---------'
!print *, 'Check points:'
!print *, r,theta
!print *, metric(1,1), metric(2,2),metric(3,3),metric(4,4),metric(4,1)
!print *, uvector(1), uvector(2), uvector(3),uvector(4)


!print *, 'check:', metric(1,1) * uvector(1)**2 &
 !                + metric(2,2) * uvector(2)**2 &
 !                + metric(3,3) * uvector(3)**2 &
 !                + metric(4,4) * uvector(4)**2 &
 !                +2.0_dp*metric(4,1)*uvector(1)*uvector(4)




g = -1.0 / ( pvector(1)*uvector(1) + pvector(2)*uvector(2) + pvector(3)*uvector(3) + pvector(4)*uvector(4))

!print *, 'u=',uvector
!print *, 'p=', pvector
!call gamma_shift(pvector,uvector,nu)
!print *, 'nu=', nu


endif


!Calculate the frequency shift as measured by the orbit
!call gamma_shift(pvector,uvector,g)


!and as measured by a cofalling ZAMO
!Create a zamo at r,theta
!call zamo(r,theta,zamo_vector)
!call gamma_shift(pvector,zamo_vector,nu)


mm = sqrt(array(1,i)**2 + a**2)
xC = mm * sin(array(2,i))*cos(array(3,i))
yC = mm * sin(array(2,i))*sin(array(3,i))
zC = mm * cos(array(2,i))

write(10,*) xC, yC ,zC,nu,&
            phiTarget, t, g, plot_status
enddo
close(10)


end subroutine ToFile




subroutine ToFile_optimisation(array)

real(kind=dp),intent(in),dimension(1000,2) :: array
character(len=200) :: Fname, ID
integer(kind=dp) :: i





ID = 'optimisation.txt'
Fname = trim(adjustl(IOpath)) // trim(adjustl(ID))

!Write
open(unit=10, file = Fname, status = 'replace',form = 'formatted')

do i=1,1000


write(10,*) array(i,1), array(i,2)

enddo




end subroutine ToFile_optimisation



subroutine calculate_contravariant_metric(r,theta,metric)
!Arguments
real(kind=dp), intent(IN) :: r, theta
real(kind=dp), intent(out), dimension(4,4) :: metric


!Internals
real(kind=dp) :: sigma, delta


sigma = r**2.0_dp + a**2*cos(theta)**2
delta = r**2.0_dp - 2.0_dp*r + a**2


metric(1,1) = - ((r**2.0_dp + a**2) + 2.0_dp*r*a**2*sin(theta)**2.0_dp/sigma)/delta
metric(2,2) = delta/sigma
metric(3,3) = 1.0_dp/sigma
metric(4,4) = (1.0_dp-2.0_dp*r/sigma)/(delta*sin(theta)**2.0_dp)

metric(1,4) = -2.0_dp*r*a/(sigma*delta)
metric(4,1) = metric(1,4)

!All other terms are zero
metric(1,2) = 0.0_dp
metric(1,3) = 0.0_dp

metric(2,1) = 0.0_dp
metric(2,3) = 0.0_dp
metric(2,4) = 0.0_dp

metric(3,1) = 0.0_dp
metric(3,2) = 0.0_dp
metric(3,4) = 0.0_dp

metric(4,2) = 0.0_dp
metric(4,3) = 0.0_dp


end subroutine calculate_contravariant_metric


subroutine calculate_covariant_metric(r,theta,metric)
!Arguments
real(kind=dp), intent(IN) :: r, theta
real(kind=dp), intent(out), dimension(4,4) :: metric


!Internals
real(kind=dp) :: sigma, delta


sigma = r**2.0_dp + a**2*cos(theta)**2
delta = r**2.0_dp - 2.0_dp*r + a**2



metric(1,1) = -(1.0_dp-2.0_dp*r/sigma)
metric(2,2) = sigma/delta
metric(3,3) = sigma
metric(4,4) = (r**2.0_dp + a**2 + 2.0_dp*r*a**2*sin(theta)**2.0_dp/sigma)*sin(theta)**2.0_dp


metric(1,4) = -2.0_dp*r*a*sin(theta)**2.0_dp/sigma
metric(4,1) = metric(1,4)

!All other terms are zero
metric(1,2) = 0.0_dp
metric(1,3) = 0.0_dp

metric(2,1) = 0.0_dp
metric(2,3) = 0.0_dp
metric(2,4) = 0.0_dp

metric(3,1) = 0.0_dp
metric(3,2) = 0.0_dp
metric(3,4) = 0.0_dp

metric(4,2) = 0.0_dp
metric(4,3) = 0.0_dp


end subroutine calculate_covariant_metric


subroutine zamo(r,theta,u)
!Create zamo at 4-position x
!Outputs contravariant 4-velocity u
real(kind=dp), intent(in) :: r, theta
real(kind=dp), intent(out),dimension(4) :: u

!Other
real(kind=dp),dimension(4,4) :: metric
real(kind=dp) :: omega, norm,delta

!Calculate the covariant metric
call calculate_covariant_metric(r,theta,metric)


!u(1) = metric(1,1) * uvector(1) + metric(3,1) * uvector(3)
!u(2) = metric(2,2) * uvector(2)
!u(3) = metric(3,3) * uvector(3)
!u(4) = metric(4,4) * uvector(4) + metric(3,1) * uvector(1)

!print *, uvector
!print *, u
!print *, u(1) * uvector(1) + u(2)*uvector(2) + u(3)*uvector(3) + u(4)*uvector(4)
!stop



!Construct the zamo
omega = -metric(1,4) / metric(4,4)
delta = r**2 -2.0_dp*r + a**2
norm = sqrt(metric(4,4)/ (delta*sin(theta)**2))

u(1) = 1.0_dp ; u(2) = 0.0_dp ; u(3) = 0.0_dp ; u(4) = omega

u = u*norm


end subroutine




!subroutine energy_shift(r,theta,Lz,freq)
!!Arguments
!real(kind=dp), intent(in) :: r,theta,Lz
!real(kind=dp), intent(out) :: freq
!!Other
!real(kind=dp) :: sigma, delta,gphph, AA, omega
!
!
!
!sigma = r**2 +a**2*cos(theta)**2
!delta = r**2 -2.0_dp*r + a**2
!
!gphph = sin(theta)**2*((r**2 +a**2)**2 - delta*a**2*sin(theta)**2)/sigma
!
!AA = sqrt(gphph / (delta*sin(theta)**2))
!
!omega = 2.0_dp*a*r / ((r**2+a**2)**2 - delta*a**2*sin(theta)**2)
!
!
!freq = -1.0_dp*AA + omega*AA*Lz
!
!end subroutine energy_shift
!
!
!
!

subroutine gamma_shift(p,u,g)
!arguments
real(kind=dp),dimension(4),intent(in) :: p,u
real(kind=dp), intent(out) :: g
!Other 

g = p(1)*u(1) + p(2)*u(2) + p(3)*u(3) + p(4)*u(4)

!g = g/p(1)

end subroutine gamma_shift

!

end module IO
