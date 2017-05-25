! расчет

subroutine calculation2D
use commonArrays
use commonVariables

implicit none

integer :: i,j,k,subit=0,ii=0
real*8 :: p_dim,r_dim,c,re1,re2,re3,re4,un,ut 
real*8 :: q_cons(4),RHS(4),LHS(4)
real*8 :: cc,cfl_check,kapa
real*8 :: qnp(3,200),delro(0:10000) 
real*4 :: vo

select case (ImplicitSchemeType)

case ('BW') ! Beam-Warming*********************************************************************************

 call derivatives				! вычисление производных

 call conditionsInlet           ! граничные условия на входе

 call conditionsOutlet          ! граничные условия на выходе

 call conditionsBound           ! граничные условия периодичноcти, на стенке

 call fluxCalc                  ! вычисление потоков

 call eigenvalues				! вычисление собственных значений

 ! вычисление слагаемого потоков в правой части

 do i=1,NNx
  do j=1,NNy
   do k=1,4

		delta_god(i,j,k)=-tau(i,j)*(flx(i+1,j,k)*deta(i+1,j)-flx(i,j,k)*deta(i,j)+fly(i,j+1,k)*dksi(i,j+1)-fly(i,j,k)*dksi(i,j))/S(i,j)
   
   enddo
  enddo
 enddo

 ! проход по направлению кси вперед
 
 do j=1,NNy  

 	dfi(0,j,1,:)=0. ! обнуление "нулевого" приращения

	do i=1,NNx


		! переход от консервативных переменных к примитивным - умножение на матрицу Т (точную - 'accu')

		call MultT(re1,re2,re3,re4,delta_god(i,j,1),delta_god(i,j,2),delta_god(i,j,3),delta_god(i,j,4),q(i,j,1),q(i,j,2),q(i,j,3),q(i,j,4),'accu')

		! вычисление правой части схемы

		RHS(1)=dqn(i,j,1)/3.+re1*2./3.
		RHS(2)=dqn(i,j,2)/3.+re2*2./3.
		RHS(3)=dqn(i,j,3)/3.+re3*2./3.
		RHS(4)=dqn(i,j,4)/3.+re4*2./3.

		! переход от примитивных переменных к характеристическим относительно направления кси

		re1=0.; re2=0.; re3=0.; re4=0.

		call MultL(re1,re2,re3,re4,RHS(1),RHS(2),RHS(3),RHS(4),q(i,j,1),q(i,j,4),yetac(i,j),xetac(i,j),'ksi')

		! вычисление правой части схемы после расщепления
   
		RHS(1)=re1+2./3.*tau(i,j)/dnoh(i,j)/Jac(i,j)*max(0.,lamksi(i-1,j,1))*dfi(i-1,j,1,1)
		RHS(2)=re2+2./3.*tau(i,j)/dnoh(i,j)/Jac(i,j)*max(0.,lamksi(i-1,j,2))*dfi(i-1,j,1,2)
		RHS(3)=re3+2./3.*tau(i,j)/dnoh(i,j)/Jac(i,j)*max(0.,lamksi(i-1,j,3))*dfi(i-1,j,1,3)
		RHS(4)=re4+2./3.*tau(i,j)/dnoh(i,j)/Jac(i,j)*max(0.,lamksi(i-1,j,4))*dfi(i-1,j,1,4)

		! вычисление левой части схемы после расщепления

		LHS(1)=1.+2./3.*tau(i,j)/dnoh(i,j)/Jac(i,j)*max(0.,lamksi(i,j,1))
		LHS(2)=1.+2./3.*tau(i,j)/dnoh(i,j)/Jac(i,j)*max(0.,lamksi(i,j,2))
		LHS(3)=1.+2./3.*tau(i,j)/dnoh(i,j)/Jac(i,j)*max(0.,lamksi(i,j,3))
		LHS(4)=1.+2./3.*tau(i,j)/dnoh(i,j)/Jac(i,j)*max(0.,lamksi(i,j,4))

		do k=1,4
			dfi(i,j,1,k)=RHS(k)/LHS(k)  ! приращения
		enddo
   
	enddo
 enddo

 ! проход по направлению кси назад

 do j=1,NNy

	dfi(NNx+1,j,2,:)=0. ! обнуление граничных приращений
 
	do i=NNx,1,-1

		! вычисление правой части схемы после расщепления

		RHS(1)=dfi(i,j,1,1)-2./3.*tau(i,j)/dnoh(i+1,j)/Jac(i,j)*min(0.,lamksi(i+1,j,1))*dfi(i+1,j,2,1)
		RHS(2)=dfi(i,j,1,2)-2./3.*tau(i,j)/dnoh(i+1,j)/Jac(i,j)*min(0.,lamksi(i+1,j,2))*dfi(i+1,j,2,2)
		RHS(3)=dfi(i,j,1,3)-2./3.*tau(i,j)/dnoh(i+1,j)/Jac(i,j)*min(0.,lamksi(i+1,j,3))*dfi(i+1,j,2,3)
		RHS(4)=dfi(i,j,1,4)-2./3.*tau(i,j)/dnoh(i+1,j)/Jac(i,j)*min(0.,lamksi(i+1,j,4))*dfi(i+1,j,2,4)
		
		! вычисление левой части схемы после расщепления

		LHS(1)=1.-2./3.*tau(i,j)/dnoh(i+1,j)/Jac(i,j)*min(0.,lamksi(i,j,1))
		LHS(2)=1.-2./3.*tau(i,j)/dnoh(i+1,j)/Jac(i,j)*min(0.,lamksi(i,j,2))
		LHS(3)=1.-2./3.*tau(i,j)/dnoh(i+1,j)/Jac(i,j)*min(0.,lamksi(i,j,3))
		LHS(4)=1.-2./3.*tau(i,j)/dnoh(i+1,j)/Jac(i,j)*min(0.,lamksi(i,j,4))

	    do k=1,4
		   dfi(i,j,2,k)=RHS(k)/LHS(k)   ! приращения
		enddo

	enddo 
 enddo

 ! переход от характеристических переменных по направлению кси к примитивным

 do i=1,NNx
  do j=1,NNy
	call multLinv(delq(i,j,1),delq(i,j,2),delq(i,j,3),delq(i,j,4),dfi(i,j,2,1),dfi(i,j,2,2),dfi(i,j,2,3),dfi(i,j,2,4),q(i,j,1),q(i,j,4),yetac(i,j),xetac(i,j),'ksi')
  enddo
 enddo

 ! проход вверх по направлению эта

 do i=1,NNx
  
	dfi(i,0,3,:)=0. ! обнуление граничных приращений

	do j=1,NNy

		! переход от примитивных переменных к характеристическим по направлению эта

		call multL(re1,re2,re3,re4,delq(i,j,1),delq(i,j,2),delq(i,j,3),delq(i,j,4),q(i,j,1),q(i,j,4),yksic(i,j),xksic(i,j),'eta')

		! вычисление правой части схемы после расщепления

		RHS(1)=re1+2./3.*tau(i,j)/dnov(i,j)/Jac(i,j)*max(0.,lameta(i,j-1,1))*dfi(i,j-1,3,1)
		RHS(2)=re2+2./3.*tau(i,j)/dnov(i,j)/Jac(i,j)*max(0.,lameta(i,j-1,2))*dfi(i,j-1,3,2)
		RHS(3)=re3+2./3.*tau(i,j)/dnov(i,j)/Jac(i,j)*max(0.,lameta(i,j-1,3))*dfi(i,j-1,3,3)
		RHS(4)=re4+2./3.*tau(i,j)/dnov(i,j)/Jac(i,j)*max(0.,lameta(i,j-1,4))*dfi(i,j-1,3,4)

		! вычисление левой части схемы после расщепления

		LHS(1)=1.+2./3.*tau(i,j)/dnov(i,j)/Jac(i,j)*max(0.,lameta(i,j,1))
		LHS(2)=1.+2./3.*tau(i,j)/dnov(i,j)/Jac(i,j)*max(0.,lameta(i,j,2))
		LHS(3)=1.+2./3.*tau(i,j)/dnov(i,j)/Jac(i,j)*max(0.,lameta(i,j,3))
		LHS(4)=1.+2./3.*tau(i,j)/dnov(i,j)/Jac(i,j)*max(0.,lameta(i,j,4))

		do k=1,4
			dfi(i,j,3,k)=RHS(k)/LHS(k)  ! приращения
		enddo

	enddo
 enddo

 ! проход вниз по направлению эта

 do i=1,NNx
  
	dfi(i,NNy+1,4,:)=0. ! обнуление граничных приращений
	
	do j=NNy,1,-1

		! вычисление правой части схемы после расщепления

		RHS(1)=dfi(i,j,3,1)-2./3.*tau(i,j)/dnov(i,j+1)/Jac(i,j)*min(0.,lameta(i,j+1,1))*dfi(i,j+1,4,1)
		RHS(2)=dfi(i,j,3,2)-2./3.*tau(i,j)/dnov(i,j+1)/Jac(i,j)*min(0.,lameta(i,j+1,2))*dfi(i,j+1,4,2)
		RHS(3)=dfi(i,j,3,3)-2./3.*tau(i,j)/dnov(i,j+1)/Jac(i,j)*min(0.,lameta(i,j+1,3))*dfi(i,j+1,4,3)
		RHS(4)=dfi(i,j,3,4)-2./3.*tau(i,j)/dnov(i,j+1)/Jac(i,j)*min(0.,lameta(i,j+1,4))*dfi(i,j+1,4,4)

		! вычисление левой части схемы после расщепления

		LHS(1)=1.-2./3.*tau(i,j)/dnov(i,j+1)/Jac(i,j)*min(0.,lameta(i,j,1))
		LHS(2)=1.-2./3.*tau(i,j)/dnov(i,j+1)/Jac(i,j)*min(0.,lameta(i,j,2))
		LHS(3)=1.-2./3.*tau(i,j)/dnov(i,j+1)/Jac(i,j)*min(0.,lameta(i,j,3))
		LHS(4)=1.-2./3.*tau(i,j)/dnov(i,j+1)/Jac(i,j)*min(0.,lameta(i,j,4))

		do k=1,4
			dfi(i,j,4,k)=RHS(k)/LHS(k) ! приращения
		enddo

	enddo
 enddo

 ! переход от характеристических переменных по направлению эта к примитивным

 do i=1,NNx
  do j=1,NNy
	call multLinv(delq(i,j,1),delq(i,j,2),delq(i,j,3),delq(i,j,4),dfi(i,j,4,1),dfi(i,j,4,2),dfi(i,j,4,3),dfi(i,j,4,4),q(i,j,1),q(i,j,4),yksic(i,j),xksic(i,j),'eta')
  enddo
 enddo

 q=q+delq   ! значение параметров на новом слое

 do i=1,NNx
  do j=1,NNy
   if (delq(i,j,1)<0.) then
    q(i,j,1)=q(i,j,1)+coeff*delq(i,j,1)**2/q(i,j,1)
   endif
   if (delq(i,j,4)<0.) then
    q(i,j,4)=q(i,j,4)+coeff*delq(i,j,4)**2/q(i,j,4)
   endif
  enddo
 enddo

 dqn=q-qn   ! временн'ое приращение

 qn=q       ! значению на старом слое присваиваем значение на новом слое

! энтропийная функция

do i=1,NNx
 do j=1,NNy
!	s_dim(i,j)=q(i,j,4)/q(i,j,1)**gam
 enddo
enddo

total=Ntime

case ('NONE')   ! явная схема **************************************************************************************


 call derivatives

 call conditionsInlet           ! вход
 call conditionsOutlet          ! выход
 call conditionsBound           ! периодичноcти, на стенке

 call fluxCalc                  ! вычисление потоков

do i=1,NNx
 do j=1,NNy
  
  ! вычисление приращений потоков (правая часть)

  do k=1,4
   delta_god(i,j,k)=-tau(i,j)*(flx(i+1,j,k)*deta(i+1,j)-flx(i,j,k)*deta(i,j)+fly(i,j+1,k)*dksi(i,j+1)-fly(i,j,k)*dksi(i,j))/S(i,j)
  enddo 

    
  ! консервативные переменные на новом слое = на старом + приращения

  q_cons(1)=q(i,j,1)+delta_god(i,j,1)
  q_cons(2)=q(i,j,2)*q(i,j,1)+delta_god(i,j,2)
  q_cons(3)=q(i,j,3)*q(i,j,1)+delta_god(i,j,3)
  q_cons(4)=q(i,j,4)/gamm+q(i,j,1)*(q(i,j,2)**2+q(i,j,3)**2)/2.+delta_god(i,j,4)

  ! переход от консервативных переменных к примитивным (без умножения на матрицу)

  q(i,j,1)=q_cons(1)
  q(i,j,2)=q_cons(2)/q_cons(1)
  q(i,j,3)=q_cons(3)/q_cons(1)
  q(i,j,4)=gamm*(q_cons(4)-(q_cons(2)**2+q_cons(3)**2)/2./q_cons(1))

  ! энтропия

!  s_dim(i,j)=q(i,j,4)/q(i,j,1)**gam

  
 enddo
enddo

dqn=q-qn ! приращение по времени

total=Ntime

case ('Newton')  ! неявная схема Ньютона  **********************************************************************************

 if ((Ntime==1.or.mod(Ntime,Nderivatives)==0.).and.(Nderivatives.ne.0)) then
  call derivatives
 endif

 10 continue

 delro(0)=0.

 subit=subit+1
 total=total+1

!write(*,*)'subiteration',subit

 if (Nderivatives==0) then
  call derivatives
 endif
 
 call conditionsInlet           ! вход
 call conditionsOutlet          ! выход
 call conditionsBound           ! периодичноcти, на стенке

 call fluxCalc                  ! вычисление потоков

!if (subit==1) then 

 call eigenvalues 

!endif
 
 do i=1,NNx
  do j=1,NNy
   do k=1,4
    delta_god(i,j,k)=-tau(i,j)*(flx(i+1,j,k)*deta(i+1,j)-flx(i,j,k)*deta(i,j)+fly(i,j+1,k)*dksi(i,j+1)-fly(i,j,k)*dksi(i,j))/S(i,j)
   enddo
  enddo
 enddo

 do i=1,NNx
  do j=1,NNy

!    if (subitMax==1) then
!		call MultT(re1,re2,re3,re4,delta_god(i,j,1),delta_god(i,j,2),delta_god(i,j,3),delta_god(i,j,4),q(i,j,1),q(i,j,2),q(i,j,3),q(i,j,4),'accu')
!	else
		call MultT(re1,re2,re3,re4,delta_god(i,j,1),delta_god(i,j,2),delta_god(i,j,3),delta_god(i,j,4),q(i,j,1),q(i,j,2),q(i,j,3),q(i,j,4),'appr')
!	endif

    delta_god(i,j,1)=re1
	delta_god(i,j,2)=re2
	delta_god(i,j,3)=re3
	delta_god(i,j,4)=re4

  enddo
 enddo

 do i=1,NNx
  do j=1,NNy
   if (timeDerivType=='3layers'.or.timeDerivType=='combine') then
    teta(i,j)=1./2.
   else if (((dqn(i,j,1)*delta_god(i,j,1)>=0.).and.timeDerivType=='combine').or.(timeDerivType=='2layers')) then
    teta(i,j)=0.
   endif
  enddo
 enddo

 
 do j=1,NNy  
   
	if (implicitInOutlet.ne.1) then 
		 dfi(0,j,1,:)=0.
	endif

	do i=1,NNx

!		if (timeDerivType=='3layers'.or.timeDerivType=='combine') then
 
			RHS(1)=s_prec*(-dq(i,j,1)+teta(i,j)*dqn(i,j,1)/(1.+teta(i,j))+delta_god(i,j,1)/(1.+teta(i,j)))
			RHS(2)=s_prec*(-dq(i,j,2)+teta(i,j)*dqn(i,j,2)/(1.+teta(i,j))+delta_god(i,j,2)/(1.+teta(i,j)))
			RHS(3)=s_prec*(-dq(i,j,3)+teta(i,j)*dqn(i,j,3)/(1.+teta(i,j))+delta_god(i,j,3)/(1.+teta(i,j)))
			RHS(4)=s_prec*(-dq(i,j,4)+teta(i,j)*dqn(i,j,4)/(1.+teta(i,j))+delta_god(i,j,4)/(1.+teta(i,j)))
   
!		else if (((dqn(i,j,1)*delta_god(i,j,1)>=0.).and.timeDerivType=='combine').or.(timeDerivType=='2layers')) then
!
!			RHS(1)=s_prec*(-dq(i,j,1)+(delta_god(i,j,1)))!+delta_godn(i,j,1))*1./2.)
!			RHS(2)=s_prec*(-dq(i,j,2)+(delta_god(i,j,2)))!+delta_godn(i,j,2))*1./2.)
!			RHS(3)=s_prec*(-dq(i,j,3)+(delta_god(i,j,3)))!+delta_godn(i,j,3))*1./2.)
!			RHS(4)=s_prec*(-dq(i,j,4)+(delta_god(i,j,4)))!+delta_godn(i,j,4))*1./2.)
   
!		endif
    
		call MultL(re1,re2,re3,re4,RHS(1),RHS(2),RHS(3),RHS(4),q(i,j,1),q(i,j,4),yetac(i,j),xetac(i,j),'ksi')
         
		if (implicitInOutlet==1.and.i==1) then
			
			dfi(0,j,1,1)=re1
			dfi(0,j,1,2)=re2
			dfi(0,j,1,3)=re3
			dfi(0,j,1,4)=re4
		
		endif

		RHS(1)=re1+tau(i,j)/dnoh(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lamksi(i-1,j,1))*dfi(i-1,j,1,1)
		RHS(2)=re2+tau(i,j)/dnoh(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lamksi(i-1,j,2))*dfi(i-1,j,1,2)
		RHS(3)=re3+tau(i,j)/dnoh(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lamksi(i-1,j,3))*dfi(i-1,j,1,3)
		RHS(4)=re4+tau(i,j)/dnoh(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lamksi(i-1,j,4))*dfi(i-1,j,1,4)
   
		LHS(1)=1.+tau(i,j)/dnoh(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lamksi(i,j,1))
		LHS(2)=1.+tau(i,j)/dnoh(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lamksi(i,j,2))
		LHS(3)=1.+tau(i,j)/dnoh(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lamksi(i,j,3))
		LHS(4)=1.+tau(i,j)/dnoh(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lamksi(i,j,4))

		do k=1,4
			dfi(i,j,1,k)=RHS(k)/LHS(k)
		enddo
		
		if (implicitInOutlet==1.and.i==1) then 
			do k=1,4
				dfi(1,j,1,k)=dfi(0,j,1,k)
			enddo
		endif
   
	enddo
 enddo

 do j=1,NNy
  
	if (implicitInOutlet.ne.1) then 
		dfi(NNx+1,j,2,:)=0.
	else
		dfi(NNx+1,j,2,:)=dfi(NNx,j,1,:)
	endif

	do i=NNx,1,-1

		RHS(1)=dfi(i,j,1,1)-tau(i,j)/dnoh(i+1,j)/Jac(i,j)/(1.+teta(i,j))*min(0.,lamksi(i+1,j,1))*dfi(i+1,j,2,1)
		RHS(2)=dfi(i,j,1,2)-tau(i,j)/dnoh(i+1,j)/Jac(i,j)/(1.+teta(i,j))*min(0.,lamksi(i+1,j,2))*dfi(i+1,j,2,2)
		RHS(3)=dfi(i,j,1,3)-tau(i,j)/dnoh(i+1,j)/Jac(i,j)/(1.+teta(i,j))*min(0.,lamksi(i+1,j,3))*dfi(i+1,j,2,3)
		RHS(4)=dfi(i,j,1,4)-tau(i,j)/dnoh(i+1,j)/Jac(i,j)/(1.+teta(i,j))*min(0.,lamksi(i+1,j,4))*dfi(i+1,j,2,4)

		LHS(1)=1.-tau(i,j)/dnoh(i+1,j)/Jac(i,j)/(1.+teta(i,j))*min(0.,lamksi(i,j,1))
		LHS(2)=1.-tau(i,j)/dnoh(i+1,j)/Jac(i,j)/(1.+teta(i,j))*min(0.,lamksi(i,j,2))
		LHS(3)=1.-tau(i,j)/dnoh(i+1,j)/Jac(i,j)/(1.+teta(i,j))*min(0.,lamksi(i,j,3))
		LHS(4)=1.-tau(i,j)/dnoh(i+1,j)/Jac(i,j)/(1.+teta(i,j))*min(0.,lamksi(i,j,4))	
		
		do k=1,4
			dfi(i,j,2,k)=RHS(k)/LHS(k)
		enddo

		if (implicitInOutlet==1.and.i==NNx) then 
  			do k=1,4
				dfi(NNx,j,2,k)=dfi(NNx+1,j,2,k)
			enddo
		endif

	enddo 
 enddo

 do i=1,NNx
  do j=1,NNy
   call multLinv(delq(i,j,1),delq(i,j,2),delq(i,j,3),delq(i,j,4),dfi(i,j,2,1),dfi(i,j,2,2),dfi(i,j,2,3),dfi(i,j,2,4),q(i,j,1),q(i,j,4),yetac(i,j),xetac(i,j),'ksi')
  enddo
 enddo

 do i=1,NNx
 
	if ((implicitWalls.ne.1).and.(implicitPeriod.ne.1)) then 
		dfi(i,0,3,:)=0.
	else
		if (i<=N1.or.i>N2) then
		    if ((Ntime==1.and.subit==1).or.(implicitPeriod.ne.1)) then
				dfi(i,0,3,:)=0.
			else
				dfi(i,0,3,:)=dfi(i,NNy,3,:)
			endif
		else
		    if (implicitWalls==1) then

				call multL(re1,re2,re3,re4,delq(i,1,1),delq(i,1,2),delq(i,1,3),delq(i,1,4),q(i,1,1),q(i,1,4),yksic(i,1),xksic(i,1),'eta')

				kapa=tau(i,1)/dnov(i,1)/Jac(i,1)/(1.+teta(i,1))

				RHS(1)=re1*(1.+kapa*max(0.,lameta(i,1,4)))+kapa*max(0.,lameta(i,0,1))*re4
				RHS(2)=re2
				RHS(3)=re3
				RHS(4)=re4*(1.+kapa*max(0.,lameta(i,1,1)))+kapa*max(0.,lameta(i,0,4))*re1
   
				LHS(1)=(1.+kapa*max(0.,lameta(i,1,4)))*(1.+kapa*max(0.,lameta(i,1,1)))-kapa**2*max(0.,lameta(i,0,1))*max(0.,lameta(i,0,4))
				LHS(2)=1.+kapa*(max(0.,lameta(i,1,2))-max(0.,lameta(i,0,2)))
				LHS(3)=LHS(2)
				LHS(4)=LHS(1)

				do k=1,4
					dfi(i,1,3,k)=RHS(k)/LHS(k)
				enddo

!				dfi(i,0,3,1)=dfi(i,1,3,4)
!				dfi(i,0,3,2)=dfi(i,1,3,2)
!				dfi(i,0,3,3)=dfi(i,1,3,3)
!				dfi(i,0,3,4)=dfi(i,1,3,1)

				dfi(i,0,3,1)=(dfi(i,1,3,1)+dfi(i,1,3,4))/2.
				dfi(i,0,3,2)=dfi(i,1,3,2)
				dfi(i,0,3,3)=dfi(i,1,3,3)
				dfi(i,0,3,4)=(dfi(i,1,3,1)+dfi(i,1,3,4))/2.


			else
				dfi(i,0,3,:)=0.
			endif

		endif
	endif

	do j=1,NNy

		call multL(re1,re2,re3,re4,delq(i,j,1),delq(i,j,2),delq(i,j,3),delq(i,j,4),q(i,j,1),q(i,j,4),yksic(i,j),xksic(i,j),'eta')

		RHS(1)=re1+tau(i,j)/dnov(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lameta(i,j-1,1))*dfi(i,j-1,3,1)
		RHS(2)=re2+tau(i,j)/dnov(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lameta(i,j-1,2))*dfi(i,j-1,3,2)
		RHS(3)=re3+tau(i,j)/dnov(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lameta(i,j-1,3))*dfi(i,j-1,3,3)
		RHS(4)=re4+tau(i,j)/dnov(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lameta(i,j-1,4))*dfi(i,j-1,3,4)

		LHS(1)=1.+tau(i,j)/dnov(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lameta(i,j,1))
		LHS(2)=1.+tau(i,j)/dnov(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lameta(i,j,2))
		LHS(3)=1.+tau(i,j)/dnov(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lameta(i,j,3))
		LHS(4)=1.+tau(i,j)/dnov(i,j)/Jac(i,j)/(1.+teta(i,j))*max(0.,lameta(i,j,4))
		
	    do k=1,4
			dfi(i,j,3,k)=RHS(k)/LHS(k)
	    enddo
       
	    if (implicitWalls==1.and.j==1.and.(i>N1.and.i<=N2)) then 
!	 		dfi(i,1,3,1)=dfi(i,0,3,4)
!			dfi(i,1,3,2)=dfi(i,0,3,2)
!			dfi(i,1,3,3)=dfi(i,0,3,3)
!			dfi(i,1,3,4)=dfi(i,0,3,1)
	    endif

	enddo
 enddo


 do i=1,NNx

	if ((implicitWalls.ne.1).and.(implicitPeriod.ne.1)) then 
		dfi(i,NNy+1,4,:)=0.
	else
		if (i<=N1.or.i>N2) then
			
			if ((Ntime==1.and.subit==1).or.(implicitPeriod.ne.1)) then
				dfi(i,NNy+1,4,:)=0.
			else
				dfi(i,NNy+1,4,:)=dfi(i,1,4,:)
			endif
		
		else
			
			if (implicitWalls==1) then

				kapa=tau(i,NNy)/dnov(i,NNy+1)/Jac(i,NNy)/(1.+teta(i,NNy))

				RHS(1)=dfi(i,NNy,3,1)*(1.-kapa*min(0.,lameta(i,NNy,4)))-kapa*min(0.,lameta(i,NNy+1,1))*dfi(i,NNy,3,4)
				RHS(2)=dfi(i,NNy,3,2)
				RHS(3)=dfi(i,NNy,3,3)
				RHS(4)=dfi(i,NNy,3,4)*(1.-kapa*min(0.,lameta(i,NNy,1)))-kapa*min(0.,lameta(i,NNy+1,4))*dfi(i,NNy,3,1)
   
				LHS(1)=(1.-kapa*min(0.,lameta(i,NNy,4)))*(1.-kapa*min(0.,lameta(i,NNy,1)))-kapa**2*min(0.,lameta(i,NNy+1,1))*min(0.,lameta(i,NNy+1,4))
				LHS(2)=1+kapa*(min(0.,lameta(i,NNy+1,2))-min(0.,lameta(i,NNy,2)))
				LHS(3)=LHS(2)
				LHS(4)=LHS(1)

				do k=1,4
					dfi(i,NNy,4,k)=RHS(k)/LHS(k)
				enddo   
			
!				dfi(i,NNy+1,4,1)=dfi(i,NNy,4,4)
!				dfi(i,NNy+1,4,2)=dfi(i,NNy,4,2)
!				dfi(i,NNy+1,4,3)=dfi(i,NNy,4,3)
!				dfi(i,NNy+1,4,4)=dfi(i,NNy,4,1)

				dfi(i,NNy+1,4,1)=(dfi(i,NNy,4,1)+dfi(i,NNy,4,4))/2.
				dfi(i,NNy+1,4,2)=dfi(i,NNy,4,2)
				dfi(i,NNy+1,4,3)=dfi(i,NNy,4,3)
				dfi(i,NNy+1,4,4)=(dfi(i,NNy,4,1)+dfi(i,NNy,4,4))/2.
			
			else
				dfi(i,NNy+1,4,:)=0.
			endif

		endif
	endif 

	do j=NNy,1,-1

		RHS(1)=dfi(i,j,3,1)-tau(i,j)/dnov(i,j+1)/Jac(i,j)/(1.+teta(i,j))*min(0.,lameta(i,j+1,1))*dfi(i,j+1,4,1)
		RHS(2)=dfi(i,j,3,2)-tau(i,j)/dnov(i,j+1)/Jac(i,j)/(1.+teta(i,j))*min(0.,lameta(i,j+1,2))*dfi(i,j+1,4,2)
		RHS(3)=dfi(i,j,3,3)-tau(i,j)/dnov(i,j+1)/Jac(i,j)/(1.+teta(i,j))*min(0.,lameta(i,j+1,3))*dfi(i,j+1,4,3)
		RHS(4)=dfi(i,j,3,4)-tau(i,j)/dnov(i,j+1)/Jac(i,j)/(1.+teta(i,j))*min(0.,lameta(i,j+1,4))*dfi(i,j+1,4,4)

		LHS(1)=1.-tau(i,j)/dnov(i,j+1)/Jac(i,j)/(1.+teta(i,j))*min(0.,lameta(i,j,1))
		LHS(2)=1.-tau(i,j)/dnov(i,j+1)/Jac(i,j)/(1.+teta(i,j))*min(0.,lameta(i,j,2))
		LHS(3)=1.-tau(i,j)/dnov(i,j+1)/Jac(i,j)/(1.+teta(i,j))*min(0.,lameta(i,j,3))
		LHS(4)=1.-tau(i,j)/dnov(i,j+1)/Jac(i,j)/(1.+teta(i,j))*min(0.,lameta(i,j,4))   

		do k=1,4
			dfi(i,j,4,k)=RHS(k)/LHS(k)
		enddo
		
	    if (implicitWalls==1.and.j==NNy.and.(i>N1.and.i<=N2)) then 
!			dfi(i,NNy,4,1)=dfi(i,NNy+1,4,4)
!			dfi(i,NNy,4,2)=dfi(i,NNy+1,4,2)
!			dfi(i,NNy,4,3)=dfi(i,NNy+1,4,3)
!			dfi(i,NNy,4,4)=dfi(i,NNy+1,4,1)
		 endif

	enddo
 enddo


 do i=1,NNx
  do j=1,NNy
   call multLinv(delq(i,j,1),delq(i,j,2),delq(i,j,3),delq(i,j,4),dfi(i,j,4,1),dfi(i,j,4,2),dfi(i,j,4,3),dfi(i,j,4,4),q(i,j,1),q(i,j,4),yksic(i,j),xksic(i,j),'eta')
  enddo
 enddo
 
 q=q+delq

 do i=1,NNx
  do j=1,NNy
   if (delq(i,j,1)<0.) then
    q(i,j,1)=q(i,j,1)+coeff*delq(i,j,1)**2/q(i,j,1)
   endif
   if (delq(i,j,4)<0.) then
    q(i,j,4)=q(i,j,4)+coeff*delq(i,j,4)**2/q(i,j,4)
   endif
  enddo
 enddo

 dq=q-qn 

delro(subit)=0. 

do i=1,NNx
 do j=1,NNy
  if(abs(dq(i,j,1))>=delro(subit)) then
   delro(subit)=abs(dq(i,j,1))
  endif
 enddo
enddo 

!if (romax>epsilon.or.Ntime<Ncontrol) then
 if (abs(delro(subit)-delro(subit-1))>epsilonSub.and.subit<SubitMax) then 
 ! if (abs(delro(subit)-delro(subit-1))>epsilonSub) then 
 ! if (subit<1) then 
   goto 10
  else
   write(*,'(a16,i2)') "subiterations= ",subit
   subit=0
  endif

qn=q

dqn=dq

!delta_godn=delta_god

end select   ! окончание расчета на одном слое

call residual

end subroutine calculation2D



