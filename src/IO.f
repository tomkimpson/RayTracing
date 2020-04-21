module IO

use parameters

implicit none

private
public ToFile, ToFile_optimisation

contains

subroutine ToFile(array,n,alpha,beta,nu_obs,Lz)
real(kind=dp), dimension(Nrows, Ncols),intent(in) :: array
integer(kind=dp), intent(in) :: n
real(kind=dp), intent(in) :: alpha,beta,nu_obs,Lz
character(len=200) :: Fname, aStr, bStr, ID, nuStr
integer(kind=dp) :: i
real(kind=dp) :: xC, yC,zC, mm
real(kind=dp) :: r,theta, pr,ptheta
real(kind=dp) :: nu


!Define savefile
write(aStr,'(F10.2)') alpha
write(bStr,'(F10.2)') beta
write(nuStr,'(F10.2)') nu_obs
ID = 'RT_alpha='//trim(adjustl(aStr))//'_beta='//trim(adjustl(bstr))//'_nu='//trim(adjustl(nuStr))//'.txt'
Fname = trim(adjustl(IOpath)) // trim(adjustl(ID))



!Write
open(unit=10, file = Fname, status = 'replace',form = 'formatted')

do i=1,n

r = array(1,i)
theta = array(2,i)

!Calculate the frequency shift as measured by a cofalling ZAMO
call energy_shift(r,theta,Lz,nu)

!Convert to Cartesian coordinates
mm = sqrt(array(1,i)**2 + a**2)
xC = mm * sin(array(2,i))*cos(array(3,i))
yC = mm * sin(array(2,i))*sin(array(3,i))
zC = mm * cos(array(2,i))

write(10,*) xC, yC ,zC,nu
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



subroutine energy_shift(r,theta,Lz,freq)
!Arguments
real(kind=dp), intent(in) :: r,theta,Lz
real(kind=dp), intent(out) :: freq
!Other
real(kind=dp) :: sigma, delta,gphph, AA, omega



sigma = r**2 +a**2*cos(theta)**2
delta = r**2 -2.0_dp*r + a**2

gphph = sin(theta)**2*((r**2 +a**2)**2 - delta*a**2*sin(theta)**2)/sigma

AA = sqrt(gphph / (delta*sin(theta)**2))

omega = 2.0_dp*a*r / ((r**2+a**2)**2 - delta*a**2*sin(theta)**2)


freq = -1.0_dp*AA + omega*AA*Lz

end subroutine energy_shift





end module IO
