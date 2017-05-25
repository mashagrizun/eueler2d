real*8 function isnas(a,b)

implicit NONE

real*8, intent(in) :: a,b
real*8 resultat

if (b*a>0.) then
 resultat=(a**2+3.*b*a)/(b+a)**2
else
 resultat=0.
endif

isnas=resultat

end function isnas