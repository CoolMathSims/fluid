MODULE sub
USE param

CONTAINS
!=======================
SUBROUTINE init !wind forced model

  hmin = 0.10

  ! grid parameters
  dx = 100
  dy = 100
  dt = 1.
  g = 9.81
  rho = 1028.0
  ah = 1e-5
  tx = 0.8



  OPEN(10,file ='topo.dat',form='formatted')
  do j = 0,ny+1
    READ(10,'(103F12.6)')(h0(j,i),i=0,nx+1)
  end do

  do j = 0,ny+1
  do i = 0,nx+1
    eta(j,i) = -MIN(0.0, h0(j,i))
    etan(j,i) = eta(j,i)
  end do
  end do

  do j = 0,ny+1
  do i = 0,nx+1
    h(j,i) = h0(j,i) + eta(j,i)
  ! wet == 1 for wet and 0 for dry
    wet(j,i) = 1
    if(h(j,i) < hmin) wet(j,i) = 0
    u(j,i) = 0.
    un(j,i) = 0.
    v(j,i) = 0.
    vn(j,i) = 0.
  end do
  end do
end SUBROUTINE init

subroutine initB
  hmin = 0.05

dx = 10.0
dy = 10.0
dt = 0.1

g = 9.81
tx = 0.!0.2
ah = 1e-5


do j = 0,ny+1
do i = 0,nx+1
  h0(j,i) = 10.0
end do
end do

do i = 0,nx+1
 h0(0,i) = -10.0
 h0(ny+1,i) = -10.0
end do



do j = 0,ny+1
 h0(j,0) = -10.0
 h0(j,nx+1) = -10.0
end do

do j = 0,ny+1
do i = 0,nx+1
  eta(j,i) = -MIN(0.0,h0(j,i))
  etan(j,i) = eta(j,i)
end do
end do

do j = 0,ny+1
do i = 0,nx+1
  h(j,i) = h0(j,i)+eta(j,i)
  wet(j,i) = 1
  if(h(j,i) < hmin) wet(j,i) = 0
  u(j,i) = 0.
  un(j,i) = 0.
  v(j,i) = 0.
  vn(j,i) = 0.
end do
end do

do j = 24,26
do i = 49,51
eta(j,i) = 1.0
end do
end do

end subroutine initB

subroutine fluid

  REAL :: du(0:ny+1,0:nx+1), dv(0:ny+1,0:nx+1)
  REAL :: advx(0:ny+1,0:nx+1),advy(0:ny+1,0:nx+1)
  REAL :: hux,huy,hvx,hvy
  REAL :: diffu(0:ny+1,0:nx+1), diffv(0:ny+1,0:nx+1)
  REAL :: term1,term2,term3,term4,term5,term6,term7,term8
  REAL :: plu,miu,plv,miv
  REAL :: upl,umi,vpl,vmi
  REAL :: us, vs, dus, dvs
  REAL :: plh,mih

  REAL :: visc
  REAL :: k1,k2,k3,k4
  REAL :: etx,ety

!start with spatial terms, then time integration
!pgrad
  do j = 1,ny
    do i = 1,nx
      du(j,i) = -g*(eta(j,i+1) - eta(j,i))/dx
      dv(j,i) = -g*(eta(j+1,i) - eta(j,i))/dy
    end do
  end do

!advective
  do j = 1,ny
    do i = 1,nx
      plu = max(u(j,i),0.0)
      miu = min(u(j,i),0.0)
      plv = max(v(j,i),0.0)
      miv = min(v(j,i),0.0)

      upl = (u(j,i+1) - u(j,i))/dx
      umi = (u(j,i) - u(j,i-1))/dx
      vpl = (u(j+1,i) - u(j,i))/dy
      vmi = (u(j,i) - u(j-1,i))/dy

      advx(j,i) = -1*((plu*umi + miu*upl) + (plv*umi + miv*vpl))

      upl = (v(j,i+1) - v(j,i))/dx
      umi = (v(j,i) - v(j,i-1))/dx
      vpl = (v(j+1,i) - v(j,i))/dy
      vmi = (v(j,i) - v(j-1,i))/dy

      advy(j,i) = -1*((plv*vmi + miv*vpl) + (plu*umi + miu*upl))
    end do
  end do

!Diffusive, no slip
  do j = 1,ny
    do i = 1,nx
      visc = ah
      upl = (u(j,i+1) - u(j,i))/dx
      umi = (u(j,i) - u(j,i-1))/dx
      !hux = (plh*umi - mih*upl)/dx
      hux = (upl - umi)
!
      upl = (u(j+1,i) - u(j,i))/dy
      umi = (u(j,i) - u(j-1,i))/dy
      if (h(j,i)<hmin)then    !no slip condition
        upl = -(u(j+1,i) - u(j,i))/dy
        umi = -(u(j,i) - u(j-1,i))/dy
      end if
      huy = (upl - umi)
!
      vpl = (v(j+1,i) - v(j,i))/dy
      vmi = (v(j,i) - v(j-1,i))/dy
      hvy = (vpl - vmi)
!
      vpl = (v(j,i+1) - v(j,i))/dx
      vmi = (v(j,i) - v(j,i-1))/dx
      if (h(j,i)<hmin)then
        vpl = -(v(j,i+1) - v(j,i))/dx
        vmi = -(v(j,i) - v(j,i-1))/dx
      end if
      hvx = (vpl - vmi)

      diffu(j,i) = (hux*visc/dx) + (huy*visc/dy)
      diffv(j,i) = (hvx*visc/dx) + (hvy*visc/dy)

      !diffu(j,i) = visc*(dudx+dudy)
      !diffv(j,i) = visc*(dvdx*dvdy)
    end do
  end do


  !prediction time integration
  do j = 1,ny
    do i = 1,nx
      un(j,i) = 0.0
      us = u(j,i)
      dus = du(j,i)
      term1 = tx/(rho*h(j,i))

      vn(j,i) = 0.0
      vs = v(j,i)
      dvs = dv(j,i)
      if(wet(j,i)==1)then
        if((wet(j,i+1)==1).or.(dus>0.0)) then !uncomment term1 for forcing
          k1 = dus + diffu(j,i) + advx(j,i)! + term1
          k2 = dus+(dt*0.5*k1) + diffu(j,i)+(0.5*dt*k1) + advx(j,i)+(dt*0.5*k1)! + term1+(dt*0.5*k1)
          k3 = dus+(dt*0.5*k2) + diffu(j,i)+(0.5*dt*k2) + advx(j,i)+(dt*0.5*k2)! + term1+(dt*0.5*k2)
          k4 = dus+(dt*k3) + diffu(j,i)+(dt*k3) + advx(j,i)+(dt*k3)! + term1+(dt*k3)
          un(j,i) = us + dt*(k1+k2+k2+k3+k3+k4)/6.0
        end if
      else
        if((wet(j,i+1)==1).and.(dus<0.0)) then
          k1 = dus + diffu(j,i) + advx(j,i)! + term1
          k2 = dus+(dt*0.5*k1) + diffu(j,i)+(0.5*dt*k1) + advx(j,i)+(dt*0.5*k1)! + term1+(dt*0.5*k1)
          k3 = dus+(dt*0.5*k2) + diffu(j,i)+(0.5*dt*k2) + advx(j,i)+(dt*0.5*k2)! + term1+(dt*0.5*k2)
          k4 = dus+(dt*k3) + diffu(j,i)+(dt*k3) + advx(j,i)+(dt*k3)! + term1+(dt*k3)
          un(j,i) = us + dt*(k1+k2+k2+k3+k3+k4)/6.0
        end if
      end if
      if(wet(j,i)==1)then
        if((wet(j+1,i)==1).or.(dvs>0.0)) then
          k1 = dvs + diffv(j,i) + advy(j,i)
          k2 = dvs+(dt*0.5*k1) + diffv(j,i)+(dt*0.5*k1) + advy(j,i)+(dt*0.5*k1)
          k3 = dvs+(dt*0.5*k2) + diffv(j,i)+(dt*0.5*k2) + advy(j,i)+(dt*0.5*k2)
          k4 = dvs+(dt*k3) + diffv(j,i)+(dt*k3) + advy(j,i)+(dt*k3)
          vn(j,i) = vs + dt*(k1+k2+k2+k3+k3+k4)/6.0
        end if
      else
        if((wet(j+1,i)==1).and.(dvs<0.0)) then
          k1 = dvs + diffv(j,i) + advy(j,i)
          k2 = dvs+(dt*0.5*k1) + diffv(j,i)+(dt*0.5*k1) + advy(j,i)+(dt*0.5*k1)
          k3 = dvs+(dt*0.5*k2) + diffv(j,i)+(dt*0.5*k2) + advy(j,i)+(dt*0.5*k2)
          k4 = dvs+(dt*k3) + diffv(j,i)+(dt*k3) + advy(j,i)+(dt*k3)
          vn(j,i) = vs + dt*(k1+k2+k2+k3+k3+k4)/6.0
        end if
      end if
    end do
  end do


  do j = 1,ny
    un(j,0) = un(j,nx)
    un(j,nx+1) = un(j,1)
    vn(j,0) = vn(j,nx)
    vn(j,nx+1) = vn(j,1)
  end do




  do j = 1,ny
    do i = 1,nx

      plh = max(h(j,i),0.0)
      mih = min(h(j,i),0.0)

      upl = (un(j,i+1) - un(j,i))/dx
      umi = (un(j,i) - un(j,i-1))/dx
      etx = (plh*umi + mih*upl)

      !plv = max(vn(j,i),0.0)
      !miv = min(vn(j,i),0.0)

      vpl = (vn(j+1,i) - vn(j,i))/dy
      vmi = (vn(j,i) - vn(j-1,i))/dy
      ety = (plh*vmi + mih*vpl)

      k1 = etx + ety
      k2 = (etx+(dt*0.5*k1)) + (ety+(dt*0.5*k1))
      k3 = (etx+(dt*0.5*k2)) + (ety+(dt*0.5*k2))
      k4 = (etx+(dt*k3)) + (ety+(dt*k3))

      etan(j,i) = eta(j,i) - dt*(k1+k2+k2+k3+k3+k4)/6.0
    end do
  end do

  do j = 1,ny
    etan(j,0) = etan(j,nx)
    etan(j,nx+1) = etan(j,1)
  end do

  !do j = 0,ny+1
  !  do i = 0,nx+1
  !    eta(j,i) = etan(j,i)
  !  end do
  !end do

end subroutine fluid

subroutine filter
  REAL :: eps
  REAL :: etNew
  REAL :: t1,t2,t3
  eps = 0.1
  do j = 0,ny+1
    do i = 0,nx+1
      if (wet(j,i)==1)then
        t1 = (1.0 - 0.25*eps*(wet(j,i+1)+wet(j,i-1)+wet(j+1,i)+wet(j-1,i)))*etan(j,i)
        t2 = 0.25*eps*(wet(j+1,i)*etan(j+1,i) + wet(j-1,i)*etan(j-1,i))
        t3 = 0.25*eps*(wet(j,i+1)*etan(j,i+1) + wet(j,i-1)*etan(j,i-1))
        eta(j,i) = etan(j,i)
      else
        eta(j,i) = etan(j,i)
      end if
!
    end do
  end do


end subroutine filter

end module sub
