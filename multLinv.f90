! умножение на матрицу перехода от характеристических переменных к примитивым

subroutine multLinv(b1,b2,b3,b4,a1,a2,a3,a4,r,p,coefy,coefx,direction)

 use commonVariables

 implicit NONE

 real*8, intent(in) :: a1,a2,a3,a4,r,p,coefy,coefx
 real*8, intent(out) :: b1,b2,b3,b4
 real*8 :: c,kron
 character*3 :: direction

 c=dsqrt(gam*p/r)

 select case (direction)

 case ('ksi') ! направление кси

	kron=1.

 case ('eta') ! направление эта

	kron=-1.

 end select

 b1= a1*r/c/2.+a3/c+a4*r/c/2.
 b2=-a1*kron*coefy/2.+a2*coefx+a4*kron*coefy/2.
 b3= a1*kron*coefx/2.+a2*coefy-a4*kron*coefx/2.
 b4= a1*r*c/2.+a4*r*c/2.

end subroutine multLinv
