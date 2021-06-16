MODULE param

INTEGER, PARAMETER :: nx = 100
INTEGER, PARAMETER :: ny = 50

REAL :: h0(0:ny+1,0:nx+1), h(0:ny+1,0:nx+1)
REAL :: eta(0:ny+1,0:nx+1),etan(0:ny+1,0:nx+1)
REAL :: u(0:ny+1,0:nx+1), un(0:ny+1,0:nx+1)
REAL :: v(0:ny+1,0:nx+1), vn(0:ny+1,0:nx+1)
REAL :: dt, dx, dy, g, ah, rho
REAL :: tx, ty, eps

INTEGER :: j,i

INTEGER :: wet(0:ny+1,0:nx+1)
REAL :: hmin



END MODULE param
