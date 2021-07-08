PROGRAM plotit

use params
use transitmod

implicit none

 INTEGER :: i, j
 REAL(8), DIMENSION(nest_nPar) :: Rvec   ! Fitted-parameter vector
 REAL(8), DIMENSION(nplen) :: resP

 ! Open up jammi_in.txt and read file
 open(2784,FILE='blends.dat',FORM='FORMATTED',STATUS='UNKNOWN')
 DO j=1,NQ
   read(2784,*) gam(j)
   gam(j) = 1.0D0/gam(j)
 END DO
 close(2784)

 ! Read-in the seriesP.dat
 open(unit=30301,file='seriesP.jam')
 DO i = 1,nplen
   read(30301,*) tp(i),fp(i),sigfp(i),epochp_seq(i)
   fpwei(i) = fp(i)/(sigfp(i)**2)
   sigfpwei(i) = 1.0D0/(sigfp(i)**2)
 END DO
 close(30301)
 
 ! Read-in the bestmode
 open(unit=30302,file='bestmode.jam')
 DO i = 1,nest_nPar
   read(30302,*) Rvec(i)
 END DO
 close(30302)

 ! === Call transit to primary transit ===
 call transit(Rvec,nplen,resP,1,&
              tp,fp,sigfp,epochp_seq,fpwei,sigfpwei,&
              NresamP,integ_bigP,.FALSE.)

END PROGRAM plotit

