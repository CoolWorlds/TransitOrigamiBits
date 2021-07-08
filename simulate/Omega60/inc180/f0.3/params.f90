! Include file for example MultiNest program 'Gaussian' (see arXiv:1001.0719)

module params
implicit none


! ==============================================================================
! Model Parameters

      	INTEGER, PARAMETER :: nplen = 70488
        REAL(8), DIMENSION(nplen) :: tp, fp, sigfp, fpwei, sigfpwei, logsigfp
        INTEGER, DIMENSION(nplen) :: epochp_seq
        INTEGER, PARAMETER :: NresamP = 30
        REAL(8), PARAMETER :: integ_bigP = 0.02043361111111111D0
        REAL(8), PARAMETER :: flick = 0.0D0
        INTEGER, PARAMETER :: NQ = 17			! Number of quarters
		
		REAL(8), DIMENSION(NQ) :: gam
        LOGICAL, PARAMETER :: globalflag = .TRUE.
        INTEGER, PARAMETER :: nest_nPar = 14
      	INTEGER, PARAMETER :: nparamorig = 14
      	INTEGER, PARAMETER :: OOTlen = 24
      	INTEGER, PARAMETER :: taulen = 0
		
        REAL(8), PARAMETER :: pivotRV = 55784.7922045D0
        REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
        REAL(8), PARAMETER :: fivehalfpi = 7.853981633974483D0
        REAL(8), PARAMETER :: twopi = 6.283185307179586D0
        REAL(8), PARAMETER :: third = 0.333333333333333333D0
        REAL(8), PARAMETER :: Grvx = 5.2864092327892624D-2 ! Grvx = (86400^2*Grv)/(3pi)
 
 ! Quarter-to-quarter start/end points
 REAL(8), DIMENSION(NQ), PARAMETER :: QXstart = (/ 54953.52862037D0, &		! 1
                                                55002.50923614D0, &		! 2
						55093.2141461D0, &		! 3
						55185.36716152D0, &		! 4
						55276.48039581D0, &		! 5
						55372.43908262D0, &		! 6
						55463.1644707D0, &		! 7
						55568.35386174D0, &		! 8
						55641.50607599D0, &		! 9
						55739.83509146D0, &		! 10
						55834.1977742D0, &		! 11
						55932.39901141D0, &		! 12
						56015.72702732D0, &		! 13
						56107.12890281D0, &		! 14
						56206.47750326D0, &		! 15
						56305.08740097D0, &		! 16
						56390.97969746D0 /)		! 17

 REAL(8), DIMENSION(NQ), PARAMETER :: QXend = (/ 54997.99334075D0, &		! 1
						55091.47732973D0, &		! 2
						55182.50650574D0, &		! 3
						55275.21348059D0, &		! 4
						55371.17221239D0, &		! 5
						55462.30628212D0, &		! 6
						55552.55903682D0, &		! 7
						55635.35547074D0, &		! 8
						55738.93604628D0, &		! 9
						55833.27831544D0, &		! 10
						55931.33647373D0, &		! 11
						56015.03229749D0, &		! 12
						56106.0663811D0, &		! 13
						56204.33204943D0, &		! 14
						56304.14747395D0, &		! 15
						56390.96969746D0, &		! 16
						66390.96969746D0 /)		! 17
 
end module params
