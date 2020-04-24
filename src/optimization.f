module optimization


use parameters 
use constants
use ray

implicit none

private 

public pattern_search, optimise_alpha_beta

contains




subroutine pattern_search(alpha,beta,nu_obs,ds,stat)

use parameters
use constants

implicit none

!Arguments
real(kind=dp), intent(inout) :: alpha,beta,ds
real(kind=dp),intent(in) :: nu_obs
integer(kind=dp), intent(out) :: stat !status - did we exit ok, or else interrupt by error/precision limit

integer(kind=dp) :: idx
real(kind=dp) :: dsR, dsL, dsU, dsD
real(kind=dp),dimension(4) :: ds_collection
real(kind=dp),dimension(4,2) :: AB_collection
real(kind=dp) :: aBest, bBest, dsBest
!Gradient alpha

call run(alpha+dg,beta,nu_obs,dsR,0)

call run(alpha-dg,beta,nu_obs,dsL,0)

call run(alpha, beta+dg, nu_obs, dsU,0)

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


if (dsBest .lt. epsilon(dsBest)) then
! Has fallen ino BH likely
print *, 'Fallen into BH'
stat = 5 !Exit condition
return
endif






if (dsBest .LT. ds) then
alpha = aBest
beta = bBest
ds = dsBest

else

dg = dg / decay_factor

if (dg .LT. epsilon(dg) .or. dsBest .eq. 0.0_dp) then
    !Dont let stepsoie get too small

    !Also guard against BH fall in

    print *, 'REACHED PRECSIION LIMIT------------------------------------'
    print *, rTarget, thetaTarget, phiTarget
    print *, 'Best ds = ', ds
    stat = 1 !Exit condition

    return
endif

endif






end subroutine pattern_search



 subroutine optimise_alpha_beta(alpha,beta,nu_obs,ds,stat,globals)
 
 use parameters
 use constants
 
 implicit none
 
 !Arguments
 real(kind=dp), intent(inout) :: alpha,beta,ds
 real(kind=dp),intent(in) :: nu_obs
 integer(kind=dp), intent(out) :: stat   
 real(kind=dp), intent(inout),dimension(4) :: globals
 
 !Other
 real(kind=dp) :: ds_alpha, ds_beta
 real(kind=dp) :: gA, gB
 real(kind=dp) :: zeta,hA,hB
 real(kind=dp) :: eta, t, ds_trial
 real(kind=dp) :: ds1, alpha1, beta1
 real(kind=dp) :: trial_alpha, trial_beta
 
 
 !Get the gradients in the alpha/beta directions
 
 

!print *, 'dg = ', bit
 
 
 !Gradient alpha
 call run(alpha+bit,beta,nu_obs,ds_alpha,0)
 gA = - ((ds_alpha - ds)/bit)
! print *, ds_alpha, ds, ds_alpha -ds    


! call run(alpha-bit,beta,nu_obs,ds_alpha,1)
! gA = - ((ds_alpha - ds)/dg)
! print *, ds_alpha -ds    





 !Gradient beta
 call run(alpha,beta+bit,nu_obs,ds_beta,0)
 gB = - ((ds_beta - ds)/bit)
 !print *, ds_beta - ds    


 !call run(alpha,beta-bit,nu_obs,ds_beta,1)
 !gB = - ((ds_beta - ds)/dg)
 !print *, ds_beta - ds    

!stop



 if (globals(1) .EQ. 0.0_dp) then
 !First go
 zeta = 0.0_dp
 else
 zeta = (gA*gA +gB*gB)/(globals(1)**2 + globals(2)**2)
 endif
 
 hA = gA +zeta*globals(3)
 hB = gB +zeta*globals(4)
 
 
 
 !Got the direction, now get the stepsizze by performing a line search
 eta = 2.0_dp !0.50_dp
 t = global_t    
 ds1 = 1e20
 

!tracking line search         
 do 
 
     !Trial values
     trial_alpha = alpha + t*hA
     trial_beta = beta + t*hB
 
     !Try a step in the direction
     call run(trial_alpha,trial_beta,nu_obs,ds_trial,0)
 
     print *, 'Trial:', ds, ds_trial,t


     if (ds_trial .lt. ds1) then   
     !Update and exit
     alpha1 = trial_alpha        
     beta1 = trial_beta        
     ds1 = ds_trial      
  
  
    else

    exit    

    endif   

    t = t*eta
 
 enddo
 
 


if (ds1 .lt. ds) then
print *, 'Good'
alpha = alpha1 ; beta = beta1 ; ds = ds1
global_t = t/5.0_dp
!global_t = 5d-9
else
!globals = 0.0_dp
global_t = 5d-9
global_counter = global_counter + 1
print *, 'bad', global_counter
endif



globals(1) = gA ; globals(2) = gB ; globals(3) = hA ; globals(4) = hB
print *, '----------------------'


if (global_counter .GT. 10) then
globals=0.0_dp 
global_counter = 0
print *, 'resetting globals'
endif

 end subroutine optimise_alpha_beta





end module optimization




