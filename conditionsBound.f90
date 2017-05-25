! граничные условия (периодичности, на стенке)

subroutine conditionsBound
use commonArrays 
use commonVariables

implicit none

integer :: i,j,k

real*8 :: un1,un2,ut1,ut2
real*8 :: q11,q12,q41,q42 
real*8 :: re1,re2,re3,re4
real*8 :: left(4),right(4)

select case (RRTypeChoice)
case ('classic')

do i=1,NNx

 if (i<=N1.or.i>N2) then ! условие периодичности ==================================================
!if (i==-1) then
   
   do k=1,4
	if (ExplicitSchemeType=='Zij') then
    	left(k)=(dqy_minus(i,NNy,k)*dvert(i,NNy)+dqt(i,NNy,k)*tau(i,NNy))/2.
		right(k)=(-dqy_plus(i,1,k)*dvert(i,1)+dqt(i,1,k)*tau(i,1))/2.
	else
    	left(k)=(dqy(i,NNy,k)*dvert(i,NNy)+dqt(i,NNy,k)*tau(i,NNy))/2.
		right(k)=(-dqy(i,1,k)*dvert(i,1)+dqt(i,1,k)*tau(i,1))/2.
	endif
   enddo

	! нормальная и касательная составляющие скорости (нормаль-внешняя) в "нулевой (ghost)" ячейке
	! (берем значения параметров из последней ячейки)

	q11= q(i,NNy,1)+left(1) 
	un1=(q(i,NNy,3)+left(3))*xksi(i,NNy+1)-(q(i,NNy,2)+left(2))*yksi(i,NNy+1) 
	ut1=(q(i,NNy,3)+left(3))*yksi(i,NNy+1)+(q(i,NNy,2)+left(2))*xksi(i,NNy+1) 
	q41= q(i,NNy,4)+left(4) 

	! нормальная и касательная составляющие скорости в первой ячейке

	q12= q(i,1,1)+right(1) 
	un2=(q(i,1,3)+right(3))*xksi(i,1)-(q(i,1,2)+right(2))*yksi(i,1)  
	ut2=(q(i,1,3)+right(3))*yksi(i,1)+(q(i,1,2)+right(2))*xksi(i,1)
	q42= q(i,1,4)+right(4) 

	call crash(3,q11,q41,i,NNy) ! проверка на возможную отрицательность
	call crash(3,q12,q42,i,1)

	! выполняется распад разрыва

	call RR1d(q11,un1,ut1,q41,q12,un2,ut2,q42)

	! запись параметров после вычисления распада разрыва на границах ячеек
   
	qy(i,1,1)= q11
	qy(i,1,2)=-un1*yksi(i,1)+ut1*xksi(i,1)   ! скорости преобразованы обратно      
	qy(i,1,3)= un1*xksi(i,1)+ut1*yksi(i,1)
	qy(i,1,4)= q41

	uny(i,1)=qy(i,1,3)*xksi(i,1)-qy(i,1,2)*yksi(i,1) ! нормальная компонента скорости

	! перенос параметров на верхнюю границу (из первой ячейки в последнюю)

	qy(i,NNy+1,1)=qy(i,1,1)
	qy(i,NNy+1,2)=qy(i,1,2)
	qy(i,NNy+1,3)=qy(i,1,3)
	qy(i,NNy+1,4)=qy(i,1,4)

	uny(i,NNy+1)=qy(i,NNy+1,3)*xksi(i,NNy+1)-qy(i,NNy+1,2)*yksi(i,NNy+1)
!if (i==-1) then
 else ! стенка ====================================================================================

  ! нижняя граница_________________________________________________________________________________
 
  ! различные случаи производных явной аппроксимации
   do k=1,4
	if (ExplicitSchemeType=='Zij') then   
		right(k)=(-dqy_plus(i,1,k)*dvert(i,1)+dqt(i,1,k)*tau(i,1))/2.
	else
    	right(k)=(-dqy(i,1,k)*dvert(i,1)+dqt(i,1,k)*tau(i,1))/2.
	endif
   enddo

  ! плотность, давление, нормальная и касательная составляющие скорости в первой ячейке

  q12= q(i,1,1)+right(1) 
  un2=(q(i,1,3)+right(3))*xksi(i,1)-(q(i,1,2)+right(2))*yksi(i,1)
  ut2=(q(i,1,3)+right(3))*yksi(i,1)+(q(i,1,2)+right(2))*xksi(i,1)
  q42= q(i,1,4)+right(4) 

  ! в "нулевой" ячейке нормальная составляющая направлена в противополжную сторону
  ! по сравнению с первой ячейкой

  q11= q12
  un1=-un2
  ut1= ut2
  q41= q42

  ! проверка на возможную отрицательность

	call crash(4,q11,q41,i,1)
	call crash(4,q12,q42,i,1)
  
  ! выполняется решение задачи распада разрыва

  call RR1d(q11,un1,ut1,q41,q12,un2,ut2,q42)

	! после распада разрыва, если поданы одинаковые параметры и противоположные нормальные компоненты скорости, 
	! нормальная компонента будет равна нулю. un_dn=0. 

  qy(i,1,1)=q11
  qy(i,1,2)=ut1*xksi(i,1)  ! скорости преобразованы обратно
  qy(i,1,3)=ut1*yksi(i,1)
  qy(i,1,4)=q41

    uny(i,1)=qy(i,1,3)*xksi(i,1)-qy(i,1,2)*yksi(i,1)  ! нормальная компонента скорости

  ! верхняя граница _______________________________________________________________________________
   
  ! различные случаи производных явной аппроксимации

   do k=1,4
  
	if (ExplicitSchemeType=='Zij') then
		 left(k)=(dqy_minus(i,NNy,k)*dvert(i,NNy)+dqt(i,NNy,k)*tau(i,NNy))/2.
	else
		 left(k)=(dqy(i,NNy,k)*dvert(i,NNy)+dqt(i,NNy,k)*tau(i,NNy))/2.
	endif
   
   enddo
  
  ! плотность, давление, нормальная и касательная составляющие скорости в последней ячейке

  q11= q(i,NNy,1)+left(1)
  un1=(q(i,NNy,3)+left(3))*xksi(i,NNy+1)-(q(i,NNy,2)+left(2))*yksi(i,NNy+1)
  ut1=(q(i,NNy,3)+left(3))*yksi(i,NNy+1)+(q(i,NNy,2)+left(2))*xksi(i,NNy+1)
  q41= q(i,NNy,4)+left(4)

  ! в "Ny+1" ячейке нормальная составляющая направлена в противополжную сторону
  ! по сравнению с последней ячейкой

  q12= q11
  un2=-un1
  ut2= ut1
  q42= q41
  		
  ! проверка на возможную отрицательность

	call crash(5,q11,q41,i,NNy)
	call crash(5,q12,q42,i,NNy)

  ! решение задачи распада разрыва

  call RR1d(q11,un1,ut1,q41,q12,un2,ut2,q42)

	! после распада разрыва, если поданы одинаковые параметры и противоположные нормальные компоненты скорости, 
	! нормальная компонента будет равна нулю. un_dn=0. 


  qy(i,NNy+1,1)=q11
  qy(i,NNy+1,2)=ut1*xksi(i,NNy+1)  ! скорости преобразованы обратно 
  qy(i,NNy+1,3)=ut1*yksi(i,NNy+1)
  qy(i,NNy+1,4)=q41

    uny(i,NNy+1)=qy(i,NNy+1,3)*xksi(i,NNy+1)-qy(i,NNy+1,2)*yksi(i,NNy+1)  ! нормальная компонента скорости

 endif
enddo

end select

end subroutine conditionsBound


