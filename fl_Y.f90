! вычисление потоков по направлению y 
subroutine fl_Y(b1,b2,b3,b4,r,u,v,p,unn,cos,sin)
use commonArrays
use commonVariables

implicit none

real*8, intent(in) :: r,u,v,p,cos,sin,unn
real*8, intent(out) :: b1,b2,b3,b4
real*8 :: e

!unn=v*cos-u*sin           !unn=v*xksi(i,j)-u*yksi(i,j)
e=p/gamm+r*(u**2+v**2)/2.  

if (label==0) then
! unn=0.
endif

b1=r*unn
b2=r*u*unn-p*sin
b3=r*v*unn+p*cos
b4=(e+p)*unn

end subroutine fl_Y