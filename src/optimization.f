module optimization


use parameters 
use constants
use ray

implicit none

private 

public pattern_search

contains




subroutine pattern_search(alpha,beta,nu_obs,ds)

use parameters
use constants

implicit none

!Arguments
real(kind=dp), intent(inout) :: alpha,beta,ds
real(kind=dp),intent(in) :: nu_obs


integer(kind=dp) :: idx
real(kind=dp) :: dsR, dsL, dsU, dsD
real(kind=dp),dimension(4) :: ds_collection
real(kind=dp),dimension(4,2) :: AB_collection
real(kind=dp) :: aBest, bBest, dsBest
!Gradient alpha

!print *, 'r1'
!print *, alpha+dg,beta
call run(alpha+dg,beta,nu_obs,dsR,0)

!print *, 'r2'
call run(alpha-dg,beta,nu_obs,dsL,0)

!print *, 'r3'
call run(alpha, beta+dg, nu_obs, dsU,0)

!print *, 'r4'
call run(alpha, beta-dg, nu_obs, dsD,0)


ds_collection(1) = dsR 
ds_collection(2) = dsL 
ds_collection(3) = dsU 
ds_collection(4) = dsD 


AB_collection(1,1) = alpha+ dg
AB_collection(1,2) = beta


AB_collection(2,1) = alpha- dg
AB_collection(2,2) = beta

AB_collection(3,1) = alpha
AB_collection(3,2) = beta + dg


AB_collection(4,1) = alpha
AB_collection(4,2) = beta - dg



idx = minloc(ds_collection,1)


aBest = AB_collection(idx,1)
bBest = AB_collection(idx,2)

dsBest = ds_collection(idx)

!print *, 'out = ',aBest, bBest, dsBest,ds,dg,decay_factor

if (dsBest .LT. ds) then
alpha = aBest
beta = bBest
ds = dsBest

else

dg = dg / decay_factor

if (dg .LT. epsilon(dg)) then
    !Dont let stepsoie get too small
    print *, 'precision limit'
    !Set ds small so outer loop exits
    ds = 1e-21
    return
endif

endif






end subroutine pattern_search



! subroutine optimise_alpha_beta(alpha,beta,nu_obs,ds,globals)
! 
! use parameters
! use constants
! 
! implicit none
! 
! !Arguments
! real(kind=dp), intent(inout) :: alpha,beta,ds
! real(kind=dp),intent(in) :: nu_obs
! real(kind=dp), intent(inout),dimension(4) :: globals
! 
! !Other
! real(kind=dp) :: ds_alpha, ds_beta
! real(kind=dp) :: gA, gB
! real(kind=dp) :: zeta,hA,hB
! real(kind=dp) :: eta, t, ds_trial
! real(kind=dp) :: ds1, alpha1, beta1
! real(kind=dp) :: trial_alpha, trial_beta
! 
! 
! !Get the gradients in the alpha/beta directions
! 
! 
! 
! 
! 
! !Gradient alpha
! call run(alpha+dg,beta,nu_obs,ds_alpha,0)
! gA = - ((ds_alpha - ds)/dg)
! 
! !Gradient beta
! call run(alpha,beta+dg,nu_obs,ds_beta,0)
! gB = - ((ds_beta - ds)/dg)
! 
! 
! if (globals(1) .EQ. 0.0_dp) then
! !First go
! zeta = 0.0_dp
! else
! zeta = (gA*gA +gB*gB)/(globals(1)**2 + globals(2)**2)
! endif
! 
! hA = gA +zeta*globals(3)
! hB = gB +zeta*globals(4)
! 
! 
! !Normalise the direction. Just tells you higher or lower info now
! hB = gB / abs(gA)
! hA = gA / abs(gA)
! 
! 
! 
! 
! !Got the direction, now get the stepsizze by performing a line search
! eta = 2.0_dp !0.50_dp
! t = global_t !5e-9
! ds1 = 1e20
! 
! 
! print *, 'start iteration'
! print *, 'gradient = ', hA
! 11 do 
! 
!     !Trial values
!     trial_alpha = alpha + t*hA
!     trial_beta = beta + t*hB
! 
!     !Try a step in the direction
!     call run(trial_alpha,trial_beta,nu_obs,ds_trial,0)
! 
!     print *, 'Tr:', trial_alpha,ds_trial,t
!     if (ds_trial .LT. ds) then
! 
!     !Update the best values
!     alpha1 = alpha+t*hA
!     beta1 = beta+t*hB
!     ds1 = ds_trial
!     global_t = t 
! 
!     exit
! 
! 
!     else
!     !Has stopped improving.
!     !Use this stepsize going forwards
!     
!     stop
! 
!     !exit
!     if (t .LT. dg / 1e6) then
!     print *, 'no iszeable change', ds1,ds
!     dg = dg / 10.0_dp
!     exit    
!     endif
! 
! 
!     t = t/eta
!     goto 11
!     endif
! 
! 
! 
! 
! 
! 
! enddo
! 
! 
! 
! 
! !Update before exiting subroutine 
! 
! if (ds1 .LT. ds) then
! globals(1) = gA ; globals(2) = gB ; globals(3) = hA ; globals(4) = hB
! alpha = alpha1 ; beta = beta1 ; ds = ds1
! 
! else
! !Reset, gives it a kick
! print *, 'Failure- needs reset'
! print *, 'Gradients', gA, gB
! print *, 'adjusting gradient stepsie:', dg
! 
! 
! if ((alpha + dg) .EQ. alpha) then
! stop
! endif
! 
! globals = 0.0_dp
! !stop
! !global_t = 1e-3 !1.0_dp !global_t/10.0_dp
! endif
! 
! 
! 
! !print *, 'out = ', alpha,beta,ds
! print *, '---------------------'
! !stop
! 
! end subroutine optimise_alpha_beta





end module optimization




