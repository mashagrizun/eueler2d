! умножение на матрицу перехода от примитивных переменных к характеристическим

subroutine multL(b1,b2,b3,b4,a1,a2,a3,a4,r,p,coefy,coefx,direction)

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

 b1=-a2*kron*coefy+a3*kron*coefx+a4/r/c
 b2= a2*coefx+a3*coefy
 b3= a1*c-a4/c
 b4= a2*kron*coefy-a3*kron*coefx+a4/r/c

end subroutine multL