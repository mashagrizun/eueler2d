! умножение на матрицу переход от консервативных переменных к примитивным
subroutine multT(b1,b2,b3,b4,a1,a2,a3,a4,r,u,v,p,markmark)

 use commonVariables

 implicit NONE

 real*8, intent(in) :: a1,a2,a3,a4,r,u,v,p
 real*8, intent(out) :: b1,b2,b3,b4
 character*4, intent(in) :: markmark
 real*8 :: rn,un,vn

 select case (markmark)

 case ('appr') ! приближенная матрица перехода

	b1= a1
	b2=-a1*u/r+a2/r
	b3=-a1*v/r+a3/r
	b4= gamm*(a1*(u**2+v**2)/2.-a2*u-a3*v+a4)

 case ('accu') ! точная матрица перехода

	b1= a1
	rn= r+b1

	b2=-a1*u/rn+a2/rn
	un= u+b2

	b3=-a1*v/rn+a3/rn
	vn= v+b3

	b4= gamm*(a1*(u*un+v*vn)/2.-a2*(u+un)/2.-a3*(v+vn)/2.+a4)

 end select

end subroutine multT