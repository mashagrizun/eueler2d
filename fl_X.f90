! вычисление потоков по направлению х
subroutine fl_X(b1,b2,b3,b4,r,u,v,p,unn,cos,sin)
use commonArrays
use commonVariables

implicit none

real*8, intent(in) :: r,u,v,p,cos,sin,unn
real*8, intent(out) :: b1,b2,b3,b4
real*8 :: e

!unn=u*cos-v*sin           !unn=u*yeta(i,j)-v*xeta(i,j)
e=p/gamm+r*(u**2+v**2)/2.

b1=r*unn
b2=r*u*unn+p*cos
b3=r*v*unn-p*sin
b4=(e+p)*unn

end subroutine fl_X
