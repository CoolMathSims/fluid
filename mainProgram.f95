PROGRAM shalWater
USE param
USE sub
!using wet/dry algorithm  and topography from Ocean Modeling for Beginners
REAL :: time,dtmax
INTEGER :: n,ntot,nout

!CALL init
CALL initB

! output of initial eta distribution
OPEN(10,file ='eta0.dat',form='formatted')
  DO j = 1,ny
    WRITE(10,'(101F12.6)')(eta(j,i),i=1,nx)
  END DO
CLOSE(10)

! output of initial layer thickness distribution
OPEN(10,file ='h0.dat',form='formatted')
  DO j = 1,ny
    WRITE(10,'(101F12.6)')(h0(j,i),i=1,nx)
  END DO
CLOSE(10)


ntot = 5000
nout = 5



OPEN(10,file ='eta.dat',form='formatted')
OPEN(20,file ='h.dat',form='formatted')
OPEN(30,file ='u.dat',form='formatted')
OPEN(40,file ='v.dat',form='formatted')


DO n = 1,ntot

time = REAL(n)*dt


CALL fluid

DO j = 1,ny
  etan(j,0) = etan(j,nx)
  etan(j,nx+1) = etan(j,1)
END DO

CALL filter

DO j = 0,ny+1
DO i = 0,nx+1
  h(j,i) = h0(j,i) + eta(j,i)
  wet(j,i) = 1
  IF(h(j,i)<hmin)wet(j,i) = 0
  u(j,i) = un(j,i)
  v(j,i) = vn(j,i)
END DO
END DO

! write to file
IF(MOD(n,nout)==0)THEN
  DO j = 1,ny
    WRITE(10,'(101F12.6)')(eta(j,i),i=1,nx)
    WRITE(20,'(101F12.6)')(h(j,i),i=1,nx)
    WRITE(30,'(101F12.6)')(u(j,i),i=1,nx)
    WRITE(40,'(101F12.6)')(v(j,i),i=1,nx)
  END DO
  WRITE(6,*)"Data output at time = ",time
ENDIF

END DO 

END PROGRAM shalWater
