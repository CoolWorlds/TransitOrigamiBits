! =====================================================================================
! LUNA version 1.4
!
! AUTHOR: David Kipping
!         Harvard-Smithsonian Center for Astrophysics
!         Email: dkipping@cfa.harvard.edu
!
! CITATION: If using this code, please cite Kipping, D., 2011, 
!           'LUNA: An algorithm for generating dynamic planet-moon transits', 
!           MNRAS, 416, 689
!
! DESCRIPTION: LUNA generates dynamically accurate lightcurves from a 
!              planet-moon pair
!         
!        File unit numbers...
!          1 - luna_in.txt
!          9 - luna.log
!         10 - full_output.dat
!         11 - lightcurve.dat
!
! CHANGES: * v1.0 is the first new version since v0.1. It is aimed at being a 
!            non-beta code. This entails removing a lot of inefficient code, 
!            comments, unnecessary checks. It also implements several more 
!            efficient versions of the LUNA eqns.
!          * v1.1 corrects some minor errors with gamma/beta
!          * v1.2 corrects psi1 = bb -> ab*DCOS(ib)
!          * v1.3 switches esinw/ecosw to e/w; removes Kstar & error call terms
!          * v1.4 includes some minor speed enhancements
!
! =====================================================================================

MODULE lunamod
use mandelmod
implicit none

CONTAINS

! ==============================================================
! === SUBROUTINE: LUNA ===
!
SUBROUTINE luna(t,Pb,T_b,p_in,ab,eb,wb,bb,Ps,phi_s,s_in,as,es,ws,is,Os,Msp_in,&
           u1,u2,fplan,reverse,nz,show,animate,limb,ItotL)

implicit none

 ! 1. LUNA controls
 INTEGER :: i, j, nz, error
 INTEGER :: show, animate, limb
 REAL(8) :: tstep
 ! Input parameters
 REAL(8), INTENT(IN) :: Pb, T_b, p_in, ab, eb, wb, bb
 REAL(8), INTENT(IN) :: Ps, phi_s, s_in, as, es, ws, is, Os, Msp_in
 REAL(8), INTENT(IN) :: u1, u2, fplan
 REAL(8) :: p, s, fmoon, Msp
 ! Preamble parameters
 REAL(8) :: T_s, rb_mid, rs_mid, ib
 REAL(8) :: cosib, sinib, cosis, sinis, coswbOs, sinwbOs
 REAL(8) :: coswb, sinwb, cosws, sinws
 REAL(8) :: pin1, pout1, sin1, sout1, psin1, psout1
 REAL(8) :: temp
 ! 2. Time/True Anomalies
 REAL(8), DIMENSION(nz) :: varrhob, varrhos
 REAL(8), DIMENSION(nz), INTENT(IN) :: t
 REAL(8), DIMENSION(nz) :: tb, ts, fb, fs
 REAL(8), DIMENSION(nz) :: cosfb, sinfb, cosfs, sinfs
 ! 3. Motion parameters
 REAL(8) :: gamma1, gamma2, beta1, beta2, psi1, epsilon2
 REAL(8), DIMENSION(nz) :: gammy, belly
 REAL(8), DIMENSION(nz) :: S_P, S_S, S_PS
 ! 4. Areal parameters
 INTEGER, DIMENSION(nz) :: prinny, subby        ! Case description
 REAL(8), DIMENSION(nz) :: x12, y12, x13, y13, x23, y23
 REAL(8), DIMENSION(nz) :: cost1, sint1, cost2, sint2
 REAL(8), DIMENSION(nz) :: c1, c2, c3
 REAL(8), DIMENSION(nz) :: A_p, A_s
 ! 5. Limb darkening
 REAL(8), DIMENSION(nz) :: IpL, IsL              ! Limb darkened output
 REAL(8), DIMENSION(nz), INTENT(OUT) :: ItotL    ! Limb darkened output
 REAL(8), DIMENSION(nz) :: IpN, IsN, ItotN       ! Non-LDed output
 ! 6. Animation
 REAL(8), DIMENSION(nz) :: Xp, Yp, Xs, Ys
 ! Constants
 REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
 REAL(8), PARAMETER :: twopi = 6.283185307179586D0
 REAL(8), PARAMETER :: halfpi = 1.5707963267948966D0
 REAL(8), PARAMETER :: piI = 0.318309886D0
 REAL(8), PARAMETER :: rad = 0.017453292519943295D0
 REAL(8), PARAMETER :: third = 0.3333333333333333D0
 REAL(8), PARAMETER :: Grv = 6.67428D-11
 LOGICAL :: verb, reverse

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! SECTION 1: COMMAND & CONTROL
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
! Open up luna.log
 verb = .FALSE. ! By default, verb 1 gives a LOT of output
  IF( verb ) THEN
   open(2,FILE='luna.log',FORM='FORMATTED',ACCESS='APPEND')
   write(2,*) ' - - - - - - - - - - - - - - - - - - - - - - - - -'
  END IF
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! SECTION 2: TIME CONVERSION
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ! SUBSECTION 2.1: Pre-amble

 error = 0
 p = p_in
 s = s_in
 Msp = Msp_in

 ! SUBSECTION 2.2: Pre-amble
 
 ! Get fmoon = true anomaly at time of inferior conjunction of moon
 fmoon = halfpi - ws
 IF( fmoon .LT. 0.0D0 ) THEN
  fmoon = fmoon + twopi
 END IF

 ! [Note, in all cases, ab&aw&as are in units of stellar radii]
 coswb = DCOS(wb); sinwb = DSIN(wb)
 cosws = DCOS(ws); sinws = DSIN(ws)
 rb_mid = ab*(1.0D0 - eb**2)/(1.0D0 + eb*sinwb)
 rs_mid = as*(1.0D0 - es**2)/(1.0D0 + es*sinws)
 cosib = bb/rb_mid
 ib = DACOS(cosib)      ! We use cosine here because ib=90 for bb=0
 sinib = DSIN(ib)
 cosis = DCOS(is); sinis = DSIN(is)
 coswbOs = DCOS(wb+Os); sinwbOs = DSIN(wb+Os)

 ! Contact points
 pin1   = (1.0D0 - p)
 pout1  = (1.0D0 + p)
 sin1   = (1.0D0 - s)
 sout1  = (1.0D0 + s)
 psin1  = (p - s)
 psout1 = (p + s)

 ! SUBSECTION 2.3: Call Kepler

 ! Offset times
 T_s = 0.5D0*piI*phi_s*Ps
 DO i=1,nz
  tb(i) = t(i) - T_b
  ts(i) = tb(i) - T_s
 END DO

 ! Kepler converts time array in f array
 call kepler(eb,wb,Pb,fplan,nz,tb,fb)
 call kepler(es,ws,Ps,fmoon,nz,ts,fs)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! SECTION 3: PLANET-MOON MOTION
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ! SUBSECTION 3.1: Pre-amble
 temp = 1.0D0 - eb**2
 DO i=1,nz
  cosfb(i) = DCOS(fb(i))
  sinfb(i) = DSIN(fb(i))
  cosfs(i) = DCOS(fs(i))
  sinfs(i) = DSIN(fs(i))
  varrhob(i) = temp/(1.0D0 + eb*cosfb(i))
 END DO
 temp = 1.0D0 - es**2
 DO i=1,nz
  varrhos(i) = temp/(1.0D0 + es*cosfs(i))
 END DO

 ! Gamma & beta definitions
 beta1  = -as*sinis*sinwbOs
 beta2  = as*coswbOs
 gamma1 = as*( cosib*sinis*coswbOs - sinib*cosis )
 gamma2 = as*cosib*sinwbOs
 ! epsilon and psi
 epsilon2 = ab
 psi1 = ab*cosib

 ! SUBSECTION 3.2: Planetary motion
 DO i=1,nz
  ! belly = beta in the guide
  belly(i) = beta1*varrhos(i)*( cosws*sinfs(i) + cosfs(i)*sinws ) + &
             beta2*varrhos(i)*( cosfs(i)*cosws - sinfs(i)*sinws )
  ! gammy = gamma in the guide
  gammy(i) = gamma1*varrhos(i)*( cosws*sinfs(i) + cosfs(i)*sinws ) + &
             gamma2*varrhos(i)*( cosfs(i)*cosws - sinfs(i)*sinws )
  ! Now S_P
  S_P(i) = ( epsilon2*varrhob(i)*( cosfb(i)*coswb - sinfb(i)*sinwb ) - &
             Msp*belly(i) )**2 + &
           ( psi1*varrhob(i)*( coswb*sinfb(i) + cosfb(i)*sinwb ) - &
             Msp*gammy(i) )**2
  S_P(i) = DSQRT(S_P(i))
  !!The following lines may be used if one wishes to fit for two-planets instead
  !!DO i=1,nz
  !! S_P(i) = rb(i)*DSQRT( 1.0D0 - ( DSIN(ib)*DSIN(fb(i)+wb) )**2 )
  !!END DO
 END DO

 ! SUBSECTION 3.3: Satellite motion
 DO i=1,nz
  S_S(i) = ( epsilon2*varrhob(i)*( cosfb(i)*coswb - sinfb(i)*sinwb ) + &
             belly(i) )**2 &
           + ( psi1*varrhob(i)*( coswb*sinfb(i) + cosfb(i)*sinwb ) + &
               gammy(i) )**2
  S_S(i) = DSQRT(S_S(i))
  !!The following lines may be used if one wishes to fit for two-planets instead
  !!DO i=1,nz
  !! S_S(i) = rs(i)*DSQRT( 1.0D0 - ( DSIN(is)*DSIN(fs(i)+ws) )**2 )
  !!END DO
 END DO

 ! SUBSECTION 3.4: PS motion
 DO i=1,nz
  S_PS(i) = (1.0D0+Msp)**2*(gammy(i)**2+belly(i)**2)
  S_PS(i) = DSQRT(S_PS(i))
  !!The following lines may be used if one wishes to fit for two-planets instead
  !!DO i=1,nz
  !!  S_PS(i) = ( rb(i)*DCOS(fb(i)+wb) - &
  !!              rs(i)*( DCOS(fs(i)+ws)*DCOS(Os) + &
  !!              DCOS(is)*DSIN(fs(i)+ws)*DSIN(Os) ) )**2 + &
  !!            ( rb(i)*DCOS(ib)*DSIN(fb(i)+wb) - &
  !!              rs(i)*DCOS(is)*DCOS(Os)*DSIN(fs(i)+ws) + &
  !!              rs(i)*DCOS(fs(i)+ws)*DSIN(Os) )**2
  !!  S_PS(i) = DSQRT(S_PS(i))
  !!END DO
 END DO

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! SECTION 4: AREA OF OCCULTATION
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ! SUBSECTION 4.1: Determine Principal Cases
 DO i=1,nz
  ! Defaults
  prinny(i) = 0; subby(i) = 0; temp = 0.0D0
  x12(i) = 0.0D0; y12(i) = 0.0D0
  x13(i) = 0.0D0; y13(i) = 0.0D0
  x23(i) = 0.0D0; y23(i) = 0.0D0
  cost1(i) = 0.0D0; sint1(i) = 0.0D0
  cost2(i) = 0.0D0; sint2(i) = 0.0D0
  c1(i) = 0.0D0 ; c2(i) = 0.0D0; c3(i) = 0.0D0
  ! Decision 1
  IF( S_P(i) .GE. pout1 ) THEN
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Tree: "p_out" => cases 1-9
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Decision 2
   IF ( S_S(i) .GE. sout1 ) THEN
    ! Sub-tree: "p_out-s_out" => cases 1,2,3 = case 1
    prinny(i) = 1
   ELSE IF ( S_S(i) .LE. sin1 ) THEN
    ! Sub-tree: "p_out-s_in" => cases 7,8,9 = case 7
    ! [cases 8,9 are unphysical => we must have case 7]
    prinny(i) = 7
   ELSE !IF(sin1.LT.S_S(i) .AND. S_S(i).LT.sout1) THEN !<-not needed since only 3 possibilities
    ! Sub-tree: "p_out-s_part" => cases 4,5,6 = 4
    ! [case 6 is unphysical, cases 4&5 give same A_S value]
    prinny(i) = 4
   END IF
  ELSE IF( S_P(i) .LE. pin1 ) THEN
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Tree: "p_in" => cases 19-27
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Decision 2
   IF ( S_S(i) .GE. sout1 ) THEN
    ! Sub-tree: "p_in-s_out" => cases 19,20,21 = case 19
    prinny(i) = 19
   ELSE IF ( S_S(i) .LE. sin1 ) THEN
    ! Sub-tree: "p_in-s_in" => cases 25,26,27
    ! Decision 3
    IF( S_PS(i) .LE. psin1 ) THEN
     prinny(i) = 27
    ELSE IF( S_PS(i) .GE. psout1 ) THEN
     prinny(i) = 25
    ELSE !IF(psin1.LT.S_PS(i) .AND. S_PS(i).LT.psout1) THEN
     prinny(i) = 26
    END IF
   ELSE !IF( sin1.LT.S_S(i) .AND. S_S(i).LT.sout1) THEN !<-not needed since only 3 possibilities
    ! Sub-tree: "p_in-s_part" => cases 22,23,24 = case 22 or 23
    ! [case 24 is unphysical]
    ! Decision 3
    IF( S_PS(i) .GE. psout1 ) THEN
     prinny(i) = 22
    ELSE !IF(psin1.LT.S_PS(i) .AND. S_PS(i).LT.psout1) THEN
     prinny(i) = 23
    END IF
   END IF
  ELSE !IF(pin1.LT.S_P(i) .AND. S_P(i).LT.pout1) THEN
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Tree: "p_part" => cases 10-19
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Decision 2
   IF ( S_S(i) .GE. sout1 ) THEN
    ! Sub-tree: "p_part-s_out" => cases 10,11,12 = case 10
    prinny(i) = 10
   ELSE IF ( S_S(i) .LE. sin1 ) THEN
    ! Sub-tree: "p_part-s_in" => cases 16,17,18
    ! Decision 3
    IF( S_PS(i) .LE. psin1 ) THEN
     prinny(i) = 18
    ELSE IF( S_PS(i) .GE. psout1 ) THEN
     prinny(i) = 16
    ELSE !IF(psin1.LT.S_PS(i) .AND. S_PS(i).LT.psout1) THEN
     prinny(i) = 17
    END IF
   ELSE !IF( sin1.LT.S_S(i) .AND. S_S(i).LT.sout1) THEN !<-not needed since only 3 possibilities
    ! Sub-tree: "p_part-s_part" => cases 13,14**,15
    ! Decision 3
    IF( S_PS(i) .LE. psin1 ) THEN
     prinny(i) = 15
    ELSE IF( S_PS(i) .GE. psout1 ) THEN
     prinny(i) = 13
    ELSE !IF(psin1.LT.S_PS(i) .AND. S_PS(i).LT.psout1) THEN
     prinny(i) = 14
     ! SUBSECTION 4.2: Determine Sub Cases
     ! Case 14 has 7 possible sub-case states.
     ! Decision 1: Evaluate condition B.
     ! Condition B requires new materials: x12, y12, cost1, sint1
     x12(i) = ( 1.0D0 - p**2 + S_P(i)**2 )/( 2.0D0*S_P(i) )
     y12(i) = ( DSQRT(2.0D0*S_P(i)**2*(1.0D0+p**2) &
              - (1.0D0-p**2)**2 - S_P(i)**4) )/( 2.0D0*S_P(i) )
     cost1(i) = ( S_P(i)**2 + S_S(i)**2 - S_PS(i)**2 ) /( 2.0D0*S_P(i)*S_S(i) )
     sint1(i) = DSQRT(1.0D0 - cost1(i)**2)
     ! Perform condition B evaluation
     IF( ( (x12(i)-S_S(i)*cost1(i))**2 & 
          + (y12(i)+S_S(i)*sint1(i))**2 ) .LT. s**2 ) THEN
      ! Condition B accepted => case 14.3a or 14.3b
      ! Perform condition F evaluation
      IF( S_P(i) .GT. 1.0D0 ) THEN
       subby(i) = 5  ! 14.3a
      ELSE !IF S_P(i) .LT. 1.0D0 ) THEN
       subby(i) = 6  ! 14.3b
      END IF
     ELSE
      ! Condition B rejected => case 14.1, 14.2 or 14.7
      ! The following material will be needed...
      y13(i) = -( DSQRT( 2.0D0*S_S(i)**2*(1.0D0+s**2) &
               - (1.0D0-s**2)**2 - S_S(i)**4 ) )/(2.0D0*S_S(i))
      temp = (1.0D0 - s**2 + S_S(i)**2)/(2.0D0*S_S(i))
      x13(i) = temp*cost1(i) - y13(i)*sint1(i)
      y13(i) = temp*sint1(i) + y13(i)*cost1(i)
      cost2(i) = -(S_P(i)**2 + S_PS(i)**2 - S_S(i)**2)/(2.0D0*S_P(i)*S_PS(i))
      sint2(i) = DSQRT(1.0D0 - cost2(i)**2)
      y23(i) = ( DSQRT(2.0D0*S_PS(i)**2*(p**2+s**2) - (p**2-s**2)**2 & 
               - S_PS(i)**4) )/(2.0D0*S_PS(i))
      temp = (p*p - s*s + S_PS(i)**2)/(2.0D0*S_PS(i))
      x23(i) = temp*cost2(i) - y23(i)*sint2(i) + S_P(i)
      y23(i) = temp*sint2(i) + y23(i)*cost2(i)
      ! Decision 2: Evaluate condition A
      ! Condition A does not require any new materials
      IF( ( (x12(i)-S_S(i)*cost1(i))**2 &
           + (y12(i)-S_S(i)*sint1(i))**2 ) .LT. s*s ) THEN
       ! Condition A accepted => case 14.1a or 14.1b
       ! Decision 3: Evaluate condition C
       ! Perform evaluation of condition C
       IF( (S_S(i)*sint1(i)) .GT. ( y13(i) &
          + ( (y23(i)-y13(i))/(x23(i)-x13(i)) )*( S_S(i)*cost1(i) - x13(i) ) ) ) THEN
        ! Condition C accepted => 14.1a
        subby(i) = 1
       ELSE
        ! Condition C rejected => 14.2b
        subby(i) = 2
       END IF
      ELSE
       ! Condition A rejected => 14.2 or 14.7
       ! Decision 3: Evaluate condition D
       IF( ((x13(i)-S_P(i))**2 + (y13(i))**2) .LT. p**2 ) THEN
        ! Condition D accepted => case 14.2a or 14.2b
        ! Decision 4: Evaluate condition E
        IF( (S_S(i)-s) .LT. (S_P(i)-p) ) THEN
         ! Condition E accepted => case 14.2a
         subby(i) = 3
        ELSE
         ! Condition E rejected => case 14.2b
         subby(i) = 4
        END IF
       ELSE
        ! Condition D rejected => case 14.7a or 14.7b
        ! Decision 4: Evaluate condition G
        IF( (x23(i)**2 + y23(i)**2) .LT. 1.0D0 ) THEN
         ! Condition G accepted => case 14.7b
         subby(i) = 8
        ELSE
         ! Condition G rejected => case 14.7a
         subby(i) = 7
        END IF
       END IF
      END IF
     END IF ! Condition B
     !
    END IF
   END IF
  END IF
 END DO

 ! SUBSECTION 4.3: Determine Area of Occultation of Planet
 ! [ *** This will eventually be deleted since MA02 will do this part ]

 DO i=1,nz
  IF( prinny(i) .LE. 9 ) THEN
   A_p(i) = 0.0D0
  ELSE IF( prinny(i) .GE. 19 ) THEN
   A_p(i) = pi*p**2
  ELSE
   call alpha(1.0D0,p,S_P(i),A_p(i))
  END IF
 END DO

 ! SUBSECTION 4.4: Determine Area of Occultation of Moon
 DO i=1,nz
 ! --------------------------------------------
 IF( prinny(i).EQ.1 .OR. prinny(i).EQ.10 .OR. &
 prinny(i).EQ.15 .OR. prinny(i).EQ.18 .OR. &
 prinny(i).EQ.19 .OR. prinny(i).EQ.27 ) THEN
  ! Zero cases... [1,2,3],[10,11,12],15,18,[19,20,21],27
  A_s(i) = 0.0D0
 ! --------------------------------------------
 ELSE IF( prinny(i).EQ.7 .OR. prinny(i).EQ.16 .OR. prinny(i).EQ.25 ) THEN
  ! Full cases... [7,8,9], 16, 25
  A_s(i) = pi*s**2
 ! --------------------------------------------
 ELSE IF( prinny(i).EQ.4 .OR. prinny(i).EQ.13 &
 .OR. prinny(i).EQ.22) THEN
  ! Alpha_S* cases... [4,5,6], 13, 22
  call alpha(1.0D0,s,S_S(i),A_s(i))
 ! --------------------------------------------
 ELSE IF( prinny(i).EQ.17 .OR. prinny(i).EQ.26 ) THEN
  ! Alpha_SP cases... 17, 26
  call alpha(p,s,S_PS(i),A_s(i))
  A_s(i) = pi*s**2 - A_s(i)
 ! --------------------------------------------
 ELSE IF( prinny(i).EQ.23 ) THEN
  ! Case 23 = alpha_S* - alpha_SP... 23
  call alpha(1.0D0,s,S_S(i),temp)
  call alpha(p,s,S_PS(i),A_s(i))
  A_s(i) = temp - A_s(i)
 ! --------------------------------------------
 ELSE IF( prinny(i).EQ.14 ) THEN
  ! Case 14... 14
  IF( subby(i).EQ.4 ) THEN  ! Case 14.2b
   A_s(i) = 0.0D0
  ELSE IF( subby(i).EQ.3 ) THEN  ! Case 14.2a
   call alpha(p,s,S_PS(i),A_s(i))
   A_s(i) = pi*s**2 - A_s(i)
  ELSE IF( subby(i).EQ.7 ) THEN ! Case 14.7a
   call alpha(1.0D0,s,S_S(i),A_s(i))
  ELSE IF( subby(i).EQ.5 ) THEN ! Case 14.3a
   call alpha(1.0D0,s,S_S(i),temp)
   call alpha(1.0D0,p,S_P(i),A_s(i))
   A_s(i) = temp - A_s(i)
  ELSE IF( subby(i).EQ.6 ) THEN ! Case 14.3b
   call alpha(p,s,S_PS(i),temp)      ! alpha_SP
   call alpha(1.0D0,p,S_P(i),A_s(i)) ! alpha_P*
   A_s(i) = A_s(i) + temp
   call alpha(1.0D0,s,S_S(i),temp) ! alpha_S*
   A_s(i) = pi*p**2 + temp - A_s(i)
  ELSE IF( subby(i).EQ.8 ) THEN ! Case 14.7b
   call alpha(1.0D0,s,S_S(i),temp)
   call alpha(p,s,S_PS(i),A_s(i))
   A_s(i) = temp - A_s(i)
  ELSE IF( subby(i).EQ.1 .OR. subby(i).EQ.2 ) THEN ! Case 14.1a&b
   ! Find chord lengths
   c3(i) = DSQRT( (x13(i)-x23(i))**2 + (y13(i)-y23(i))**2 )
   c1(i) = DSQRT( (x12(i)-x13(i))**2 + (y12(i)-y13(i))**2 )
   c2(i) = DSQRT( (x12(i)-x23(i))**2 + (y12(i)-y23(i))**2 )
   A_s(i) = 0.25D0*DSQRT( (c1(i)+c2(i)+c3(i))*(c2(i)+c3(i)-c1(i))*&
           (c1(i)+c3(i)-c2(i))*(c1(i)+c2(i)-c3(i)) ) + &
           DASIN(0.5D0*c1(i)) + p**2*DASIN(0.5D0*c2(i)/p) + &
           s**2*DASIN(0.5D0*c3(i)/s) - &
           0.25D0*c1(i)*DSQRT(4.0D0-c1(i)**2) - &
           0.25D0*c2(i)*DSQRT(4.0D0*p**2-c2(i)**2)
   IF( subby(i).EQ.1 ) THEN ! Case 14.1a
    A_s(i) = A_s(i) - 0.25D0*c3(i)*DSQRT(4.0D0*s**2-c3(i)**2)
   ELSE                     ! Case 14.1b
    A_s(i) = A_s(i) + 0.25D0*c3(i)*DSQRT(4.0D0*s**2-c3(i)**2)
   END IF
   call alpha(1.0D0,s,S_S(i),temp)
   A_s(i) = temp - A_s(i)
  ELSE
   !write(2,*) 'ERROR: no sub case identified for point ',i
  END IF
 ! --------------------------------------------
 ELSE
  !write(2,*) 'ERROR: no principal case identified for point ',i
 END IF
 ! Error catcher
 IF( A_s(i) .LT. 0.0D0 ) THEN
  !write(2,*) 'ERROR: Negative A_s for point ',i,' case ',prinny(i),' sub',subby(i)
  !write(2,*) 'A_s(i) = ',A_s(i)
  !write(2,*) 'S_P(i) = ',S_P(i),'; S_S(i) = ',S_S(i),';S_PS(i) = ',S_PS(i)
  !write(2,*) 'p = ',p,'; s = ',s,'; swapped = ',swapped
  A_s(i) = 0.0D0 ! Override implemented
 END IF
 END DO

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! SECTION 5: LIMB-DARKENING
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 ! Get non-LDed lightcurves for S and tot...
 DO i=1,nz
  IsN(i) = 1.0D0 - (A_s(i)*piI)
  ItotN(i) = 1.0D0 - ( (A_p(i)+A_s(i))*piI )
 END DO

 ! See if we want LD or not...
 IF( limb .EQ. 1 ) THEN
  ! SUBSECTION 5.1: Planetary limb darkening
  ! Call occultquad from mandelmod module (Mandel-Agol code)
  call occultquad(S_P,u1,u2,p,IpL,IpN,nz)
  ! IpL = LDed flux; IpN = non-LDed flux
  !
  ! SUBSECTION 5.2: Lunar limb darkening
  ! Call ldsmall subroutine for quadratic limb darkening
  call ldsmall(S_S,A_s,0.0D0,u1+2.0D0*u2,0.0D0,-u2,s,IsL,nz)
  !
 ELSE
  !'Limb darkening de-activated'
  DO i=1,nz
   IpN(i) = 1.0D0 - A_p(i)*piI  ! Non-LDed flux has not been calc'd
   IpL(i) = IpN(i)
   IsL(i) = 1.0D0 - A_s(i)*piI
  END DO
 END IF
 
 IF( s < 1.0D-6 ) THEN
   IsL(:) = 1.0D0
 END IF

 ! Combine
 IF( reverse ) THEN
   DO i=1,nz
    !ItotL(i) = 1.0D0 - ( (1.0D0-IpL(i)) + (1.0D0-IsL(i)) )
    ItotL(i) = IpL(i) + 1.0D0 - IsL(i)
   END DO
 ELSE
   DO i=1,nz
    !ItotL(i) = 1.0D0 - ( (1.0D0-IpL(i)) + (1.0D0-IsL(i)) )
    ItotL(i) = IpL(i) - 1.0D0 + IsL(i) 
   END DO
 END IF

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! SECTION 6: OUTPUT
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ! Test output: this is useful for debugging, but not normally used
! open(10,FILE='full_output.dat',FORM='FORMATTED',STATUS='UNKNOWN')
! DO i=1,nz
!  write(10,*) i,t(i),tb(i),ts(i),fb(i),fs(i),rb(i),rs(i),&
!              S_P(i),S_S(i),S_PS(i),&
!              prinny(i),subby(i),A_p(i),A_s(i)
! END DO
! close(10) 

 ! Test output: this is the outputted lightcurve
 ! The 'animate' option priduces the columns needed to generate an animation
 ! of the orbital motion

! ! REMOVE COMMENT SECTION FOR ANIMATION
! IF( show .EQ. 1 ) THEN
!  IF( animate .EQ. 1 ) THEN
!   ! Do animate
!   open(11,FILE='lightcurve.dat',FORM='FORMATTED',STATUS='UNKNOWN')
!   DO i=1,nz
!    ! Planet's motion
!    Xp(i) = rb(i)*DCOS(fb(i)+wb) + Msp*rs(i)*( -DCOS(wb+Os)*DCOS(fs(i)+ws) + &
!            DSIN(is)*DSIN(fs(i)+ws)*DSIN(wb+Os) )
!    Yp(i) = rb(i)*DCOS(ib)*DSIN(fb(i)+wb) - Msp*rs(i)*( &
!            DCOS(ib)*DCOS(wb+Os)*DCOS(fs(i)+ws) + &
!            DSIN(fs(i)+ws)*( DCOS(is)*DSIN(ib) + &
!            DCOS(ib)*DSIN(is)*DCOS(wb+Os) ) )
!  !! Use these lines if doing two-planet fit
!  !!Xp(i) = rb(i)*DCOS(fb(i)+wb)
!  !!Yp(i) = rb(i)*DCOS(ib)*DSIN(fb(i)+wb)
!    ! Moon's motion
!    Xs(i) = rb(i)*DCOS(fb(i)+wb) + rs(i)*(-DSIN(is)*DSIN(fs(i)+ws)*DSIN(wb+Os) + &
!            DCOS(fs(i)+ws)*DCOS(wb+Os) )
!    Ys(i) = rb(i)*DCOS(ib)*DSIN(fb(i)+wb) + rs(i)*( &
!            DSIN(fs(i)+ws)*( DCOS(is)*DSIN(ib) + &
!            DCOS(ib)*DSIN(is)*DCOS(wb+Os) ) + &
!            DCOS(fs(i)+ws)*DCOS(ib)*DSIN(wb+Os) )
!  !! Use these lines if doing two-planet fit
!  !!Xs(i) = rs(i)*( DCOS(fs(i)+ws)*DCOS(Os) + DCOS(is)*DSIN(fs(i)+ws)*DSIN(Os) )
!  !!Ys(i) = rs(i)*( DCOS(is)*DCOS(Os)*DSIN(fs(i)+ws) - DCOS(fs(i)+ws)*DSIN(Os) )
!    write(11,*) t(i),ItotL(i),IpL(i),IsL(i),&
!                ItotN(i),IpN(i),IsN(i),&
!                S_P(i),S_S(i),S_PS(i),&
!                Xp(i),Yp(i),Xs(i),Ys(i)
!   END DO
!   close(11)
!  ELSE
!   ! Do not animate
!   open(11,FILE='lightcurve.dat',FORM='FORMATTED',STATUS='UNKNOWN')
!   DO i=1,nz
!    write(11,*) t(i),&
!                ItotL(i),IpL(i),IsL(i),&
!                ItotN(i),IpN(i),IsN(i)
!   END DO
!   close(11)
!  END IF
! END IF
! ! END OF REMOVE COMMENT SECTION FOR ANIMATION

 ! Close the log file
 IF( verb ) THEN
  close(2)
 END IF

 END SUBROUTINE luna

! ==============================================================
! === SUBROUTINE: KEPLER ===
!
 SUBROUTINE kepler(e,wrad,Pdays,f_ref,n,t,f_)
 ! Solves Kepler's equation and calculates f(t)

 implicit none
 INTEGER :: i, j, n!, tt
 REAL(8), INTENT(IN) :: e, wrad, Pdays, f_ref
 REAL(8) :: E_ref, M_ref, ek, Pk, toler
 INTEGER :: ok
 REAL(8), DIMENSION(n), INTENT(IN) :: t
 REAL(8), DIMENSION(n) :: M_, E_, E_p
 REAL(8), DIMENSION(n), INTENT(OUT) :: f_
 REAL(8), PARAMETER :: halfpi = 1.5707963267948966D0
 REAL(8), PARAMETER :: twopi = 6.2831853071795864769D0
     
 ! Time is in days, so all references to P must be Pdays
 ! The following values are predefined for speed
 ek = DSQRT((1.0D0+e)/(1.0D0-e))
 Pk = twopi/Pdays
 IF( e .LT. 0.9D-4 ) THEN
   E_ref = f_ref
   M_ref = E_ref
   ! 2) Solve Kepler's equation
   DO i=1,n
     M_(i) = Pk*t(i) + M_ref
     E_(i) = M_(i)
     f_(i) = E_(i)
   END DO
 ELSE
   ! If e>0, we must first solve f(t).
   E_ref = 2.0D0*DATAN((1.0D0/ek)*DTAN(f_ref*0.5D0))
   IF(E_ref .LT. -halfpi) THEN
     E_ref = E_ref + twopi
   END IF
   M_ref = E_ref - e*DSIN(E_ref)
   ! 2) Solve Kepler's equation
   toler = 1.1574074074074074D-4*Pk ! Set to 10 second tolerance, divide by 10 for 1 second
   DO i=1,n
     M_(i) = Pk*t(i) + M_ref
   END DO
   E_(:) = M_(:)
   DO i=1,n
     ok = 0
     DO WHILE ( ok .EQ. 0 )
       E_p(i) = E_(i) - ((E_(i) - e*DSIN(E_(i)) &
                - M_(i))/(1.0D0-e*DCOS(E_(i))))
       IF( DABS(E_p(i)-E_(i)) .LT. toler ) THEN
         ok = 1
       END IF
       E_(i) = E_p(i)
     END DO
     f_(i) = 2.0D0*DATAN(ek*DTAN(E_(i)*0.5D0))
   END DO
 END IF

 END SUBROUTINE kepler
! =======================================================

! ==============================================================
! === SUBROUTINE: ALPHA ===
!
 SUBROUTINE alpha(Rb,Rs,S,alph)
 ! Computes area of intersection (alpha) between two circles 
 ! of radii R and r separated by distance S.

 implicit none
 REAL(8), INTENT(IN) :: Rb, Rs, S
 REAL(8) :: kap0, kap1, kap2
 REAL(8), INTENT(OUT) :: alph
     
 kap0 = DACOS( (S**2 + Rs**2 - Rb**2)/(2.0D0*S*Rs) )
 kap1 = DACOS( (S**2 + Rb**2 - Rs**2)/(2.0D0*S*Rb) )
 kap2 = DSQRT(0.25D0*(4.0D0*S**2*Rb**2 - (Rb**2+S**2-Rs**2)**2))
 
 alph = Rs**2*kap0 + Rb**2*kap1 - kap2

 END SUBROUTINE alpha
! =======================================================

! ==============================================================
! === SUBROUTINE: LDSMALL ===
!
 SUBROUTINE ldsmall(S,Ar_occ,c1,c2,c3,c4,r,I0,n)
 ! Computes limb darkened flux assuming small-planet approximation

 implicit none
 REAL(8), DIMENSION(n), INTENT(IN) :: S, Ar_occ
 REAL(8), INTENT(IN) :: c1, c2, c3, c4, r
 INTEGER, INTENT(IN) :: n
 INTEGER :: i
 REAL(8) :: Ftot
 REAL(8), DIMENSION(n) :: am, bm, amR, bmR
 REAL(8), DIMENSION(n) :: Ar_ann, Fl_ann
 REAL(8), DIMENSION(n), INTENT(OUT) :: I0
 REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
 REAL(8), PARAMETER :: third = 0.33333333333333333333D0
 REAL(8), PARAMETER :: seventh = 0.1428571428571428D0
 
 Ftot = 1.0D0 - 0.2D0*c1 - third*c2 - 3.0D0*seventh*c3 - 0.5D0*c4

 DO i=1,n
 am(i) = (S(i)-r)**2
 bm(i) = (S(i)+r)**2
 amR(i) = DSQRT(DSQRT(1.0D0-am(i)))
 bmR(i) = DSQRT(DSQRT(1.0D0-bm(i)))
 ! 4 CASES
  IF( S(i).GT.(1.0D0+r) ) THEN
   ! CASE I
   Ar_ann(i) = 0.0D0
   Fl_ann(i) = 0.0D0
   I0(i) = 1.0D0
  ELSE IF( S(i).GT.r .AND. S(i).LT.(1.0D0-r) ) THEN 
   ! CASE III
   ! Area of annulus is given by...
   Ar_ann(i) = pi*(bm(i) - am(i))
   ! Flux of annulus is given by...
   Fl_ann(i) = (am(i)-bm(i))*(c1+c2+c3+c4-1.0D0) + &
               0.8D0*c1*amr(i)**5 + 2.0D0*third*c2*amr(i)**6 + &
               4.0D0*seventh*c3*amr(i)**7 + 0.5D0*c4*amr(i)**8 - &
               0.8D0*c1*bmr(i)**5 - 2.0D0*third*c2*bmr(i)**6 - &
               4.0D0*seventh*c3*bmr(i)**7 - 0.5D0*c4*bmr(i)**8
    ! Intensity...
    I0(i) = 1.0D0 - (Ar_occ(i)/Ar_ann(i))*(Fl_ann(i)/Ftot)
  ELSE IF( S(i).LT.r ) THEN 
   ! CASE IX
   ! Area of annulus is given by...
   Ar_ann(i) = pi*bm(i)
   ! Flux of annulus is given by...
   Fl_ann(i) = -bm(i)*(c1+c2+c3+c4-1.0D0) + &
               0.8D0*c1 + 2.0D0*third*c2 + 4.0D0*seventh*c3 + 0.5D0*c4 - &
               0.8D0*c1*bmr(i)**5 - 2.0D0*third*c2*bmr(i)**6 - &
               4.0D0*seventh*c3*bmr(i)**7 - 0.5D0*c4*bmr(i)**8
   ! Intensity...    
   I0(i) = 1.0D0 - (Ar_occ(i)/Ar_ann(i))*(Fl_ann(i)/Ftot)
  ELSE !IF( S(i).LT.(1.0D0+r) .AND. S(i).LT.(1.0D0-r) ) THEN
   ! CASE II
   ! Area of annulus is given by...
   Ar_ann(i) = pi*(1.0D0-am(i))
   ! Flux of annulus is given by...
   Fl_ann(i) = (am(i)-1.0D0)*(c1+c2+c3+c4-1.0D0) + &
               0.8D0*c1*amr(i)**5 + 2.0D0*third*c2*amr(i)**6 + &
               4.0D0*seventh*c3*amr(i)**7 + 0.5D0*c4*amr(i)**8
   ! Intensity...    
   I0(i) = 1.0D0 - (Ar_occ(i)/Ar_ann(i))*(Fl_ann(i)/Ftot)
  END IF
 END DO

 END SUBROUTINE ldsmall
! =======================================================

END MODULE lunamod
