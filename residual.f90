subroutine residual
use commonArrays
use commonVariables

implicit none

integer :: i,j,ii=0


RoMaxMax=1.		 
RoMaximum=0.

do i=1,NNx
 do j=1,NNy

  if(abs(dqn(i,j,1))>=RoMaximum) then

	RoMaximum=abs(dqn(i,j,1))

  endif	

 enddo
enddo 

if (Ntime==1.or.RoMaximum>RoMaxMax) RoMaxMax=RoMaximum

RoMaximum=RoMaximum/RoMaxMax

open (7,file='RoMaximum.dat')
open (88,file='RoMaximum10.dat')

if (mark==1) then
	do
		read (7,*,end=10) RoMaximum
	enddo
10	continue		
	do
		read (88,*,end=20) RoMaximum10
	enddo
20 continue
endif

write(7,*) RoMaximum

if (mod(Ntime,ResPeriod)==0) write(88,*) RoMaximum

! вычисление максимальной разности значений плотности через заданное количество временных шагов (Ncontrol)

ii=ii+1

! запись в массив значенний на каждом 1 слое

if (ii==1) then
 do i=1,NNx
  do j=1,NNy

   ro1(i,j)=q(i,j,1)

  enddo
 enddo
endif

! запись в массив значений через заданное количество шагов Ncontrol

if (ii==Ncontrol) then
 
 ii=0 ! счетчик обнуляется
 
 do i=1,NNx
  do j=1,NNy
 
	ro2(i,j)=q(i,j,1)

  enddo
 enddo

 romax=0. ! максимальная разность

 do i=1,NNx
  do j=1,NNy

   if(abs(ro1(i,j)-ro2(i,j))>=romax) then

		romax=abs(ro1(i,j)-ro2(i,j))	

   endif

  enddo
 enddo 

 write(*,*)

write(*,'(a22,xxxxxxxxxxx,d15.8)') 'Max Residual NControl', romax


write(*,*)
write(*,'(a6,xxxxxxxxxxxxxx,f5.3)') 'Gamma', gam
write(*,'(a5,xxxxxxxxxxxxxx,f8.3)') 'Rg  ', rg
write(*,*)
write(*,'(a6,xxxxxxxxxxxxx,d13.6)') 'T_in ', t_in
write(*,'(a6,xxxxxxxxxxxxx,d13.6)') 'P_in ', p_in
write(*,'(a10,xxxxxxxxxx,f5.2)') 'Alpha_in ', alpha_in/pi*180.
write(*,'(a7,xxxxxxxxxxxx,d13.6)') 'P_out ', p_out
write(*,*)
write(*,'(a4,xxxxxxxxxxxxxxx,f5.2)') 'CFL', cfl
write(*,*)
write(*,'(a12,xxxxxxxx,i1)') 'Minimzation', charmin
write(*,'(a15,xxxxx,f5.3)') 'Preconditioner', s_prec
write(*,*)
write(*,'(a19,x,i1)') 'Implicit In/Outlet', implicitInOutlet		
write(*,'(a17,xxx,i1)') 'Implicit Periods', implicitPeriod		 
write(*,'(a15,xxxxx,i1)') 'Implicit Walls', implicitWalls		 
write(*,*)
write(*,'(a13,xxxxxxx,f5.3)') 'Square Coeff', coeff
write(*,'(a8,xxxxxxxxxxx,d13.6)') 'Epsilon', epsilon
write(*,*)
write(*,'(a10,xxxxxxxxxx,i9)') 'Ntime Max', NtimeMax
write(*,'(a9,xxxxxxxxxxx,i5)') 'NControl', Ncontrol
write(*,'(a19,x,i5)') 'NDerivatives(Newt)', Nderivatives
write(*,'(a11,xxxxxxxxx,i5)') 'Res Period', ResPeriod
write(*,'(a15,xxxxx,i1)') 'SubitMax(Newt)', SubitMax
write(*,*)
write(*,'(a17,xx,d13.6)') 'EpsilonSub(Newt)', epsilonSub
write(*,*)
write(*,'(a15,xxxxx,i1)') 'TimeStepChoice', timeStepChoice
write(*,'(a11,xxxxxxxxx,a5)') 'GridChoice', GridChoice
write(*,'(a13,xxxxxxx,a7)') 'RRTypeChoice', RRTypeChoice
write(*,*)
write(*,'(a19,x,a3)') 'ExplicitSchemeType', ExplicitSchemeType
write(*,'(a19,x,a5)') 'ImplicitSchemeType', ImplicitSchemeType
write(*,'(a14,xxxxxx,a7)') 'TimeDerivType',timeDerivType
write(*,*)

 call PrintFromF

endif

end subroutine residual