MODULE transitmod

use params
use lunamod
use jasminemod
implicit none
      
contains
      
      
!=======================================================================
SUBROUTINE transit(Rin,nz,ressy,show,&
                   tobs,fobs,sigfobs,epoch_seq,fobswei,sigfobswei,&
                   Nresam,integ_big,seccy)
         
 implicit none
      
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 REAL(8), DIMENSION(nest_nPar) :: Rin           ! R in vector
 INTEGER, INTENT(IN) :: nz, Nresam, show
 LOGICAL, INTENT(IN) :: seccy           
 REAL(8), INTENT(IN) :: integ_big
 REAL(8) :: chi2                                ! merit function of fit
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 INTEGER :: i, j
 ! Data terms
 REAL(8), DIMENSION(nz), INTENT(IN) :: tobs, fobs, sigfobs, fobswei, sigfobswei
 INTEGER, DIMENSION(nz), INTENT(IN) :: epoch_seq
 REAL(8), DIMENSION(nz) :: tdiff, flux, ressy
 REAL(8), DIMENSION(OOTlen) :: OOTvec, weivec
 REAL(8), DIMENSION(taulen) :: tauvec
 !!REAL(8), DIMENSION(OOSlen) :: OOSvec
 ! Parameters
 REAL(8) :: p, rhomstar, bp, Pdays, gamglobal, samrecip, u1, tmid, w1, w2!, u1pu2
 REAL(8) :: Ps, phi_s, Rsp, asp_p, es, ws, is, Os, Msp
 REAL(8) :: rhoplan, rhosat
 ! Derived parameters
 REAL(8) :: u2, wrad, wavP, wavS, e, aR, as, s, rhoP, rhoS
 REAL(8) :: secpri, tsec, secphas, fpri
 REAL(8) :: tT, tF, DurP, t12, t34, rhostar, ideg, RVoffset, s14
 LOGICAL :: process
 ! Blending stuff
 REAL(8), DIMENSION(nz) :: gammy
 REAL(8), DIMENSION(nz*Nresam) :: gammy_big
 REAL(8) :: gamrecip, sam
 ! Explosion variables
 INTEGER :: k_big, k_max, nz_big
 REAL(8) :: Ndash, integ_small
 ! Explosion flag array
 INTEGER, DIMENSION(nz) :: explo
 ! Unflattened arrays
 REAL(8), DIMENSION(Nresam,nz) :: t_arr, flux_arr
 ! Flattened arrays
 REAL(8), DIMENSION(nz*Nresam) :: t_big
 REAL(8), DIMENSION(nz*Nresam) :: mulimb0!, b0
 INTEGER, DIMENSION(nz*Nresam) :: epoch_seq_big
 ! show=1 variables
 REAL(8), DIMENSION(nz) :: res
 REAL(8) :: time, epoch, epochmid, tmidfold, tsecfold
 LOGICAL :: reverse

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 1.0 DECLARATIONS ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - -

 ! PLANET TERMS
 p = Rin(1)   ! Planet's size
 rhomstar = Rin(2)  ! rho_{*}
 bp = Rin(3)  ! Barycentre's impact parameter
 Pdays = Rin(4)     ! Barycentre's period [days]
 gamglobal = 1.0D0       ! Blending factor
 samrecip = 1.0D0 !DABS(Rin(7)) ! Fp/F*
 w1 = Rin(6)    ! Limb darkening w1
 w2 = Rin(7) ! Limb darkening w2
 tmid = Rin(5)     ! Barycentre's transit time
 e = 2.0D-8 ! Barycentre's e
 wrad = 0.7853981633974483D0 ! Barycentre's w DATAN2(hb,kb)

 ! MOON TERMS
 Ps = Rin(8)			! Moon's period [days]
 asp_p = Rin(9)		! a_{SP}/R_P
 phi_s = Rin(10)		! Moon's phase angle [radians]
 is = Rin(11)		        ! Moon's inclination
 Os = 150.0D0*pi/180.0D0 !Rin(12) 			! Moon's longitude of ascending node [radians]
 es = 2.0D-8			! Moon's h
 ws = 0.7853981633974483D0	! Moon's k
 Msp = Rin(13) 			!(DSQRT((rhomsat/rhomplan)**3))*Rsp**3 ! Moon/Planet mass
 Rsp = Rin(14)			! Moon's radius
 IF( Rsp .LT. 0.0D0 ) THEN
   reverse = .TRUE.
   Rsp = -Rsp
 ELSE
   reverse = .FALSE.
 END IF

 ! Correct the moon's inclination
 IF( is .LT. 1.0D0 .AND. is .GE. -1.0D0 ) THEN
   is = DACOS(is)
 ELSE IF( is .LT. 3.0D0 .AND. is .GE. 1.0D0 ) THEN
   is = DACOS(is - 2.0D0) + pi
 ELSE
   is = 1984.0D0 ! error flag
 END IF

 ! tauarray
 IF( globalflag ) THEN
  ! donothing
 ELSE
   DO i=1,taulen
     tauvec(i) = Rin(nparamorig+i)
   END DO
 END IF

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 2.0 CONVERT FITTED PARAMETERS ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 nz_big = nz*Nresam

 ! Simple ones
 s = Rsp*p

 ! Calculate wavP
 wavP = fivehalfpi - wrad
 IF( wavP .GT. twopi ) THEN
   wavP = wavP - twopi
 END IF
 wavS = wavP + pi
 IF( wavS .GT. twopi ) THEN
   wavS = wavS - twopi
 END IF
 ! Calculate varrhoP and varrhoS
 rhoP = 1.0D0 - e**2
 rhoS = rhoP/(1.0D0 + e*DCOS(wavS))
 rhoP = rhoP/(1.0D0 + e*DCOS(wavP))

 ! Convert rhoplan to as
 as = (asp_p*p)/(1.0D0+Msp)

 ! Convert rhomstar to aR (keep units in days)
 aR = rhomstar*Grvx
 aR = aR + ( (as*(1.0D0+Msp))**3 )/( Ps**2 )
 aR = ( aR*Pdays**2 )**third

 ! Get u2
 u2 = DSQRT(w1) !u1+u2
 u1 = 2.0D0*u2*w2
 u2 = u2 - u1

 ! Override secondary eclipses to have no LD
 IF( seccy ) THEN
   u1 = 0.0D0
   u2 = 0.0D0
 END IF

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 2.1 REJECT TRIALS ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 process = .TRUE.

 ! Unphysical planet density
 rhoplan = Grvx*(1.0D0+Msp)*Ps**2
 rhoplan = asp_p**3/rhoplan
 IF( rhoplan .GT. 150000.0D0 ) THEN
   process = .FALSE.
 ELSE IF( rhoplan .LT. 30.0D0 ) THEN
   process = .FALSE.
 END IF
 
! ! Unphysical moon density
! IF( DABS(Rsp) .GT. 1.0D-9 ) THEN
!   rhosat = (rhoplan*Msp)/(DABS(Rsp)**3)
!   IF( rhosat .GT. 20000.0D0 )THEN
!     process = .FALSE.
!   ELSE IF( rhoplan .LT. 500.0D0 ) THEN
!	 process = .FALSE.
!   END IF
! END IF

 ! Contact planet-moon binaries
 IF( asp_p*(1.0D0-es) .LT. 2.0D0 ) THEN
   process = .FALSE.
 END IF

 ! Contact star-planet binaries                                                
 IF( aR*(1.0D0-e) .LT. 2.0D0 ) THEN
   process = .FALSE.
 END IF
 
! ! Exo-Hill-sphere moons
! IF( (3.0D0*(Ps/Pdays)**2)**third .GT. 0.9309D0*(1.0D0-1.0764D0*e-0.9812D0*es) ) THEN
!   process = .FALSE.
! END IF

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 3.0 TIME ARRAY ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 IF( process ) THEN

 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 ! 3.1 Offset time array for mid-transit time
 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

 DO i = 1,nz
   IF( globalflag .OR. seccy ) THEN
     ! Use a global tmid time
     tdiff(i) = tobs(i) - tmid
   ELSE
     ! Use individual transit times
     tdiff(i) = tobs(i) - tauvec(epoch_seq(i))
   END IF
 END DO

 ! Assign quarter-to-quarter blending factors
 gammy(:) = 1.0D0
 DO i=1,nz
   DO j=1,NQ
     IF( tobs(i) .GE. QXstart(j) .AND. tobs(i) .LT. QXend(j) ) THEN
       gammy(i) = gam(j)
     END IF
   END DO
 END DO

 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 ! 3.2 Jasmine Call
 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~    

 ! We call Jasmine here so that we may implement
 ! selecive resampling

 ! Jasmine calculates various durations in an exact manner
 IF( seccy ) THEN
   ! Use Jasmine's secpri value to get tsec & secphas
   call jasmine(bp,p,e,wrad,aR,Pdays,0,&
                tT,tF,DurP,t12,t34,rhostar,&
                ideg,secpri,s14,RVoffset,fpri)
 ELSE
   secpri = 0.5D0*86400.0D0*Pdays
   fpri = 1.5707963267948966D0 - wrad
 END IF
 tsec = (secpri/86400.0D0) + tmid
 secphas = secpri/(86400.0D0*Pdays)

 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 ! 3.3 Explode the time array
 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

 ! We now have to 'explode' the time array
 ! This is part of the process for accounting for the integration time
 ! of each time stamp.
 ! CAVEATS:
 ! * If m=0, no integration is required.
 ! * If we are in OOT, no integation is required (we define intransit times stamps
 !   as those < 1.1*0.5*(t_T + integration_time) (the 1.1 gives a 10% safety net)
 ! * Points which are exploded are assigned a logic flag 1 in the exploded(i) array

 IF( Nresam .GT. 1 ) THEN ! We only need to do this if Nresam > 1
  ! Stage one, create 2D time array
  k_big = 0
  Ndash = 0.5D0*(Nresam+1.0D0)     ! These two are commonly used
  integ_small = integ_big/Nresam   ! later, so easier to define here
  DO i=1,nz
   ! You add a 2nd condition here eg selective resampling, SC/LC mixed data, etc
   IF( sigfobs(i) .LE. flick ) THEN ! All data before this point is LC
    explo(i) = 1 ! Explosion is occuring, flag it
    DO j=1,Nresam
     ! Performing explosion
     k_big = k_big + 1
     t_arr(j,i) = tdiff(i) + (j-Ndash)*integ_small
     ! Stage two, flatten the 2D array into a 1D array
     t_big(k_big) = t_arr(j,i)
     gammy_big(k_big) = gammy(i)
     epoch_seq_big(k_big) = epoch_seq(i)
    END DO
   ELSE
    k_big = k_big + 1
    explo(i) = 0 ! No explosion occured for this timestamp
    t_big(k_big) = tdiff(i)
    gammy_big(k_big) = gammy(i)
    epoch_seq_big(k_big) = epoch_seq(i)
   END IF
  END DO
  k_max = k_big
 ELSE
  ! Infinitessimal integration time => go as normal
  DO i=1,nz
   t_big(i) = tdiff(i)
   gammy_big(i) = gammy(i)
   explo(i) = 0
   epoch_seq_big(i) = epoch_seq(i)
  END DO
  k_max = nz
 END IF

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 4.0 GENERATE LIGHTCURVE ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 ! 4.1 Main call to LUNA
 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

 call luna(t_big,Pdays,0.0D0,p,aR,e,wrad,bp,&
      Ps,phi_s,s,as,es,ws,is,Os,Msp,&
      u1,u2,fpri,reverse,&
      k_max,0,0,1,mulimb0)

 !  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 ! 4.2 Transformations
 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

 ! Define reciprocals of various parameters to speed up computation
 sam = 1.0D0/samrecip
 gamrecip = 1.0D0/gamglobal ! gamglobal = blending factor
 ! PERFORM MODEL -> OBSERVATIONS TRANSFORMATIONS
 IF( seccy ) THEN
   DO i=1,k_max
   ! i) Secondary Eclipse Transformation
     mulimb0(i) = mulimb0(i) + sam - 1.0D0
     mulimb0(i) = mulimb0(i)*samrecip
   END DO
 END IF
 DO i=1,k_max
   ! ii) Blending Transformation
   mulimb0(i) = mulimb0(i) + gammy_big(i) - 1.0D0
   mulimb0(i) = mulimb0(i)/gammy_big(i)
   mulimb0(i) = mulimb0(i) + gamglobal - 1.0D0
   mulimb0(i) = mulimb0(i)*gamrecip
 END DO

 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 ! 4.3 Implode the flux array
 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

 ! Implode the flux array, using standard binning
 ! First stage is to un-flatten the flux array from a
 ! a 1D vector to a 2D array

 IF( Nresam .GT. 1 ) THEN  ! Nresam>1 => implosions required
  k_big = 0
  DO i=1,nz
   ! For ith point, is an implosion needed?
   IF( explo(i) == 1 ) THEN
    ! Stage one, un-flatten the flux array
    DO j=1,Nresam
     k_big = k_big + 1
     flux_arr(j,i) = mulimb0(k_big)
    END DO
    ! Stage two, perform binning
    flux(i) = 0.0D0
    DO j=1,Nresam
     flux(i) = flux(i) + flux_arr(j,i)
    END DO
    flux(i) = flux(i)/Nresam
   ELSE
    k_big = k_big + 1
    flux(i) = mulimb0(k_big)
   END IF
  END DO
 ELSE
 ! Infinitessimal integration time => go as normal
  DO i=1,nz
   flux(i) = mulimb0(i)
  END DO
 END IF

 ! Linear optimization of OOTsfobswei, sigfobswei
 OOTvec(:) = 0.0D0; weivec(:) = 0.0D0
 DO j=1,OOTlen
   DO i=1,nz
     IF( epoch_seq(i) .EQ. j ) THEN
       OOTvec(j) = OOTvec(j) + flux(i)*fobswei(i)
       weivec(j) = weivec(j) + flux(i)**2*sigfobswei(i)
     END IF
   END DO
   IF( weivec(j) .NE. 0.0D0 ) THEN
     OOTvec(j) = OOTvec(j)/weivec(j)
   ELSE
     OOTvec(j) = 1984.0D0 ! error flag
   END IF
 END DO

 OOTvec(:) = 1.0D0
 ! Now normalize
 DO i=1,nz
   ! iii) OOT Transformation
   flux(i) = flux(i)*OOTvec(epoch_seq(i))
 END DO

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 5.0 PRINT RESULTS ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 res(:) = 0.0D0
 ! PRIMARY TRANSIT
 IF (show == 1) then
 ! First we do full lightcurve
   open(unit=91,file='PRI_full.jam')
   DO i=1,nz
     IF( globalflag .OR. seccy ) THEN
       time = tdiff(i) + tmid
     ELSE
       time = tdiff(i) + tauvec(epoch_seq(i))
     END IF
     res(i) = ( fobs(i) - flux(i) )*1.0D6
     write(91,*) time,fobs(i),sigfobs(i),flux(i)
   END DO
   close(91)
 ! Next we do the the folded LC
   epochmid = tmid/Pdays
   tmidfold = tmid - epochmid*Pdays
   open(unit=92,file='PRI_fold.jam')
   DO i=1,nz
     IF( globalflag .OR. seccy ) THEN
       time = tdiff(i) + tmid
     ELSE
       time = tdiff(i) + tauvec(epoch_seq(i))
     END IF
     epoch = time/(Pdays)
     time = time - epoch*Pdays
     time = time - tmidfold
     IF( time .GE. 0.5D0*Pdays) THEN
       time = time - Pdays
     ELSE IF( time .LE. (-0.5D0*Pdays)) THEN
       time = time + Pdays
     END IF
     write(92,*) time,fobs(i),sigfobs(i),flux(i)
   END DO
   close(92) 
 END IF
 ! SECONDARY ECLIPSE
 IF (show == 2) THEN
   ! First we do full lightcurve
   open(unit=93,file='SEC_full.jam')
   DO i=1,nz
     time = tdiff(i) + tmid
     res(i) = ( fobs(i) - flux(i) )*1.0D6
     write(93,*) time,fobs(i),sigfobs(i),flux(i)
   END DO
   close(93)
   ! Next we do the the folded LC
   epochmid = tsec/Pdays
   tsecfold = tsec - epochmid*Pdays
   open(unit=94,file='SEC_fold.jam')
   DO i=1,nz
     time = tdiff(i)+tmid
     epoch = time/(Pdays)
     time = time - epoch*Pdays
     time = time - tsecfold
     IF( time .GE. 0.5D0*Pdays) THEN
       time = time - Pdays
     ELSE IF( time .LE. (-0.5D0*Pdays)) THEN
       time = time + Pdays
     END IF
     write(94,*) time,fobs(i),sigfobs(i),flux(i)
   ENDDO
   close(94) 
 END IF

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 6.0 COMPUTE CHI^2 ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 ! 6.1 Basic chi^2
 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

 ! Calculated flux = flux(i)
 ! Observed flux = fobs(i)
 ! Observed error = sigfobs(i)
 chi2=0.0D0
 DO i=1,nz
   chi2 = chi2 + ((flux(i) - fobs(i))/sigfobs(i))**2
 END DO
! write(*,*) 'chi2 = ',chi2

 DO i=1,nz
   ressy(i) = flux(i) - fobs(i)
 END DO

 ELSE 
   ! process = .FALSE.
   ressy(:) = 1.0D0
 END IF

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 7.0 CLOSE PROGRAM ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

END SUBROUTINE transit
!=======================================================================

END MODULE transitmod
