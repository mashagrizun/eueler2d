real*8 function minmod(a,b)

implicit NONE

real*8, intent(in) :: a,b
real*8 resultat

if (a*b > 0.) then
  if (dabs(a) < dabs(b)) then
	resultat=a
  else
	resultat=b
  endif 
else
  resultat=0.
endif	

minmod=resultat

end function minmod