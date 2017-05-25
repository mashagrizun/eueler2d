subroutine multAeta(b1,b2,b3,b4,a1,a2,a3,a4,r,u,v,p,coefy,coefx)

use commonVariables
implicit none
real*8, intent(in)  :: a1,a2,a3,a4,r,u,v,p,coefy,coefx 
real*8, intent(out) :: b1,b2,b3,b4 
real*8 :: uneta

uneta=-u*coefy+v*coefx

b1=uneta*a1-r*coefy*a2+r*coefx*a3
b2=uneta*a2-coefy*a4/r
b3=uneta*a3+coefx*a4/r
b4=-gam*p*coefy*a2+gam*p*coefx*a3+uneta*a4

end subroutine multAeta