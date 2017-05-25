! вычисление производных 

subroutine derivatives
use commonVariables
use commonArrays

implicit none

integer :: i,j,k
real*8 :: re1,re2,re3,re4,minmod,isnas,rre1,rre2,rre3,rre4,hdzl
real*8 :: un1,un2,ut1,ut2,c2
real*8 :: q11,q12,q41,q42 ,right(4),left(4)
real*8 :: koNy1,koNy2,ko1,ko2,ko

rq=0.
rrq=0.
dqx=0.
dqx_minus=0.
dqx_plus=0.
dqtx=0.
dqt=0.

! вычисление произодных по направлению кси===============================================================================

! первые разности

if (ExplicitSchemeType=='God') then ! Годунов------------------------------------------------------------------------------------------------------
 
 dqx=0.

else

! первые разности

do j=1,NNy
  
	if (Ntime==1) then  ! если слой по времени - первый, значения на границах берем из значений в ближайших ячейках
	do k=1,4
		qx(1,j,k)=q(1,j,k)
		qx(NNx+1,j,k)=q(NNx,j,k)
	enddo

	endif

	do k=1,4

		rq(1,j,k)=2.*(q(1,j,k)-qx(1,j,k))/dh(1,j)
		rq(NNx+1,j,k)=2.*(qx(NNx+1,j,k)-q(NNx,j,k))/dh(NNx,j)

	enddo	
	
	c2=(gam*(qx(1,j,4)+q(1,j,4))/(qx(1,j,1)+q(1,j,1)))**2
	rq(1,j,1)=rq(1,j,4)-c2*rq(1,j,1)

	c2=(gam*(qx(NNx+1,j,4)+q(NNx,j,4))/(qx(NNx+1,j,1)+q(NNx,j,1)))**2
	rq(NNx+1,j,1)=rq(NNx+1,j,4)-c2*rq(NNx+1,j,1)

		
	do i=2,NNx
	 do k=1,4
		rq(i,j,k)=2.*(q(i,j,k)-q(i-1,j,k))/(dh(i,j)+dh(i-1,j))
	 enddo

	c2=(gam*(q(i,j,4)+q(i-1,j,4))/(q(i,j,1)+q(i-1,j,1)))**2
	rq(i,j,1)=rq(i,j,4)-c2*rq(i,j,1)

	enddo
enddo


! переход к характеристическим переменным, если задано-------------------------------------------------------------------

if (charmin==1) then
do i=1,NNx
 do j=1,NNy

	call multL(re1,re2,re3,re4,rq(i,j,1),rq(i,j,2),rq(i,j,3),rq(i,j,4),q(i,j,1),q(i,j,4),yetac(i,j),xetac(i,j),'ksi')

	rq(i,j,1)=re1
	rq(i,j,2)=re2
	rq(i,j,3)=re3
	rq(i,j,4)=re4

  enddo
 enddo
endif

select case (ExplicitSchemeType)  ! выбор типа явной аппроксимации=============================================================

case ('TVD') ! Колган ------------------------------------------------------------------------------------------------------

do i=1,NNx
 do j=1,NNy

	do k=1,4

		dqx(i,j,k)=minmod(rq(i,j,k),rq(i+1,j,k))

	enddo

	c2=(gam*q(i,j,4)/q(i,j,1))**2
	dqx(i,j,1)=(dqx(i,j,4)-dqx(i,j,1))/c2
  
 enddo
enddo

case ('ENO') ! ENO-схема ------------------------------------------------------------------------------------------------------

! вторые разности (производные)
 do j=1,NNy
  do i=2,NNx-1
	do k=1,4

		rrq(i,j,k)=4./(dh(i+1,j)+2.*dh(i,j)+dh(i-1,j))*(rq(i+1,j,k)-rq(i,j,k))

	enddo
  enddo
  do k=1,4
	
	rrq(1,j,k)=4./(dh(2,j)+2.*dh(1,j))*(rq(2,j,k)-rq(1,j,k))
	rrq(0,j,k)=rrq(1,j,k)

	rrq(NNx,j,k)=4./(2.*dh(NNx,j)+dh(NNx-1,j))*(rq(NNx+1,j,k)-rq(NNx,j,k))
	rrq(NNx+1,j,k)=rrq(NNx,j,k)

  enddo
 enddo

! вычисление производных

do j=1,NNy
 do i=2,NNx-1

  do k=1,4

	dqx(i,j,k)=minmod(rq(i,j,k)+1./2.*(dh(i,j)+dh(i-1,j))/2.*minmod(rrq(i,j,k),rrq(i-1,j,k)),rq(i+1,j,k)-1./2.*(dh(i+1,j)+dh(i,j))/2.*minmod(rrq(i+1,j,k),rrq(i,j,k)))

  enddo

  do k=1,4

	dqx(1,j,k)=minmod(rq(1,j,k)+1./2.*dh(1,j)/2.*minmod(rrq(1,j,k),rrq(0,j,k)),rq(2,j,k)-1./2.*(dh(2,j)+dh(1,j))/2.*minmod(rrq(2,j,k),rrq(1,j,k)))

	dqx(NNx,j,k)=minmod(rq(NNx,j,k)+1./2.*(dh(NNx,j)+dh(NNx-1,j))/2.*minmod(rrq(NNx,j,k),rrq(NNx-1,j,k)),rq(NNx+1,j,k)-1./2.*dh(NNx,j)/2.*minmod(rrq(NNx+1,j,k),rrq(NNx,j,k)))
	
  enddo

 enddo
enddo


case ('Zij') ! ISNAS Зийлемы------------------------------------------------------------------------------------------------------

do j=1,NNy

 do k=1,4
 
  dqx_minus(NNx,j,k)=0. 

 enddo

  do i=1,NNx-1
   do k=1,4

     dqx_minus(i,j,k)=isnas(rq(i,j,k),rq(i+1,j,k))*rq(i+1,j,k)

   enddo
 enddo
enddo   

do j=1,NNy
  
 do k=1,4
  dqx_plus(NNx,j,k)=0.
 enddo

  do i=1,NNx-1	 
   do k=1,4

      dqx_plus(i,j,k)=isnas(rq(i,j,k),rq(i+1,j,k))*rq(i,j,k)

   enddo
 enddo
enddo

end select

! если задана минимизация характеристических производных, переход обратно примитивным переменным

if (charmin==1) then
 do i=1,NNx
  do j=1,NNy

   ! производные ENO,TVD

   call multLinv(re1,re2,re3,re4,dqx(i,j,1),dqx(i,j,2),dqx(i,j,3),dqx(i,j,4),q(i,j,1),q(i,j,4),yetac(i,j),xetac(i,j),'ksi')

   dqx(i,j,1)=re1
   dqx(i,j,2)=re2
   dqx(i,j,3)=re3
   dqx(i,j,4)=re4

   ! "правая" и "левая" производные ISNAS
   
   call multLinv(re1,re2,re3,re4,dqx_minus(i,j,1),dqx_minus(i,j,2),dqx_minus(i,j,3),dqx_minus(i,j,4),q(i,j,1),q(i,j,4),yetac(i,j),xetac(i,j),'ksi')

   dqx_minus(i,j,1)=re1
   dqx_minus(i,j,2)=re2
   dqx_minus(i,j,3)=re3
   dqx_minus(i,j,4)=re4

   call multLinv(re1,re2,re3,re4,dqx_plus(i,j,1),dqx_plus(i,j,2),dqx_plus(i,j,3),dqx_plus(i,j,4),q(i,j,1),q(i,j,4),yetac(i,j),xetac(i,j),'ksi')

   dqx_plus(i,j,1)=re1
   dqx_plus(i,j,2)=re2
   dqx_plus(i,j,3)=re3
   dqx_plus(i,j,4)=re4

  enddo
 enddo
endif

endif

! вычисление производных по времени ========================================================================================================================

if ((ExplicitSchemeType=='God').or.(ImplicitSchemeType.ne.'NONE')) then 

	dqtx=0.

else

do i=1,NNx
 do j=1,NNy

  if (ImplicitSchemeType=='NONE') then

	if (ExplicitSchemeType.ne.'Zij') then

      call multAksi(dqtx(i,j,1),dqtx(i,j,2),dqtx(i,j,3),dqtx(i,j,4),dqx(i,j,1),dqx(i,j,2),dqx(i,j,3),dqx(i,j,4),q(i,j,1),q(i,j,2),q(i,j,3),q(i,j,4),yetac(i,j),xetac(i,j))

      do k=1,4
        dqtx(i,j,k)=-dqtx(i,j,k)
      enddo

	else if (ExplicitSchemeType=='Zij') then
	  
	  re1=(dqx_minus(i,j,1)+dqx_plus(i,j,1))/2.
	  re2=(dqx_minus(i,j,2)+dqx_plus(i,j,2))/2.
	  re3=(dqx_minus(i,j,3)+dqx_plus(i,j,3))/2.
	  re4=(dqx_minus(i,j,4)+dqx_plus(i,j,4))/2.
      
	  call multAksi(dqtx(i,j,1),dqtx(i,j,2),dqtx(i,j,3),dqtx(i,j,4),re1,re2,re3,re4,q(i,j,1),q(i,j,2),q(i,j,3),q(i,j,4),yetac(i,j),xetac(i,j))
      
	  do k=1,4
        dqtx(i,j,k)=-dqtx(i,j,k)
      enddo
	
	endif
  else
   do k=1,4

	dqtx(i,j,k)=0.  ! если схема неявная, берем равными нулю
   
   enddo

  endif

 enddo
enddo
endif

! вычисление произодных по направлению эта====================================================================================================

rq=0.
rrq=0.
dqy=0.
dqy_minus=0.
dqy_plus=0.
dqty=0.


if (ExplicitSchemeType=='God') then !Годунов------------------------------------------------------------------------------------------------------

 dqy=0.

else

! первые разности

do i=1,NNx

 if (i<=N1.or.i>N2) then ! на границе периодичности
  
  do k=1,4
	rq(i,1,k)=2.*(q(i,1,k)-q(i,NNy,k))/(dvert(i,1)+dvert(i,NNy))
	rq(i,NNy+1,k)=rq(i,1,k)
  enddo

  	c2=(gam*(q(i,1,4)+q(i,NNy,4))/(q(i,1,1)+q(i,NNy,1)))**2
	rq(i,1,1)=rq(i,1,4)-c2*rq(i,1,1)

	rq(i,NNy+1,1)=rq(i,1,1)

 else ! на стенке
  
	do k=1,4
		right(k) = 2.*(q(i,2,k)-q(i,1,k))/(dvert(i,2)+dvert(i,1))
	enddo

	c2=(gam*(qy(i,1,4)+q(i,1,4))/(qy(i,1,1)+q(i,1,1)))**2
	if ((right(4)-c2*right(1))<=0.) then 
		do k=1,4
			right(k)=0.
		enddo
	endif

	c2=(gam*q(i,1,4)/q(i,1,1))**2
	right(1)=(right(4)-right(1))/c2

  q12= q(i,1,1)-right(1)*dvert(i,1)/2.
  un2=(q(i,1,3)-right(3)*dvert(i,1)/2.)*xksi(i,1)-(q(i,1,2)-right(2)*dvert(i,1)/2.)*yksi(i,1)
  ut2=(q(i,1,3)-right(3)*dvert(i,1)/2.)*yksi(i,1)+(q(i,1,2)-right(2)*dvert(i,1)/2.)*xksi(i,1)
  q42= q(i,1,4)-right(4)*dvert(i,1)/2. 

  
  q11= q12
  un1=-un2
  ut1= ut2
  q41= q42

  call RR1d(q11,un1,ut1,q41,q12,un2,ut2,q42)

  qy(i,1,1)=q11
  qy(i,1,2)=ut1*xksi(i,1)  ! скорости преобразованы обратно
  qy(i,1,3)=ut1*yksi(i,1)
  qy(i,1,4)=q41


	do k=1,4
		left(k) = 2.*(q(i,NNy,k)-q(i,NNy-1,k))/(dvert(i,NNy)+dvert(i,NNy-1))
	enddo

	c2=(gam*(qy(i,NNy+1,4)+q(i,NNy,4))/(qy(i,NNy+1,1)+q(i,NNy,1)))**2
	if ((left(4)-c2*left(1))>=0.) then
		do k=1,4
		   left(k)=0.
		enddo
	endif

	c2=(gam*q(i,NNy,4)/q(i,NNy,1))**2
	left(1)=(left(4)-left(1))/c2

  q11= q(i,NNy,1)+left(1)*dvert(i,NNy)/2.
  un1=(q(i,NNy,3)+left(3)*dvert(i,NNy)/2.)*xksi(i,NNy+1)-(q(i,NNy,2)+left(2)*dvert(i,NNy)/2.)*yksi(i,NNy+1)
  ut1=(q(i,NNy,3)+left(3)*dvert(i,NNy)/2.)*yksi(i,NNy+1)+(q(i,NNy,2)+left(2)*dvert(i,NNy)/2.)*xksi(i,NNy+1)
  q41= q(i,NNy,4)+left(4)*dvert(i,NNy)/2.

  q12= q11
  un2=-un1
  ut2= ut1
  q42= q41
  		
  call RR1d(q11,un1,ut1,q41,q12,un2,ut2,q42)


  qy(i,NNy+1,1)=q11
  qy(i,NNy+1,2)=ut1*xksi(i,NNy+1)  ! скорости преобразованы обратно 
  qy(i,NNy+1,3)=ut1*yksi(i,NNy+1)
  qy(i,NNy+1,4)=q41


  do k=1,4
	rq(i,1,k)=2.*(q(i,1,k)-qy(i,1,k))/dvert(i,1)
	rq(i,NNy+1,k)=2.*(qy(i,NNy+1,k)-q(i,NNy,k))/dvert(i,NNy)
  enddo

    c2=(gam*(q(i,1,4)+qy(i,1,4))/(q(i,1,1)+qy(i,1,1)))**2
	rq(i,1,1)=rq(i,1,4)-c2*rq(i,1,1)

	c2=(gam*(q(i,NNy,4)+qy(i,NNy+1,4))/(q(i,NNy,1)+qy(i,NNy+1,1)))**2
	rq(i,NNy+1,1)=rq(i,NNy+1,4)-c2*rq(i,NNy+1,1)

 endif

 do j=2,NNy
  do k=1,4
	rq(i,j,k)=2.*(q(i,j,k)-q(i,j-1,k))/(dvert(i,j)+dvert(i,j-1))
  enddo

   	c2=(gam*(q(i,j,4)+q(i,j-1,4))/(q(i,j,1)+q(i,j-1,1)))**2
	rq(i,j,1)=rq(i,j,4)-c2*rq(i,j,1)

 enddo
enddo


! если задана минимизация характеристических производных, переход к характ. переменным

if (charmin==1) then
do i=1,NNx
 do j=1,NNy

	call multL(re1,re2,re3,re4,rq(i,j,1),rq(i,j,2),rq(i,j,3),rq(i,j,4),q(i,j,1),q(i,j,4),yksic(i,j),xksic(i,j),'eta')

	rq(i,j,1)=re1
	rq(i,j,2)=re2
	rq(i,j,3)=re3
	rq(i,j,4)=re4

  enddo

	call multL(re1,re2,re3,re4,rq(i,NNy+1,1),rq(i,NNy+1,2),rq(i,NNy+1,3),rq(i,NNy+1,4),q(i,NNy,1),q(i,NNy,4),yksic(i,NNy),xksic(i,NNy),'eta')

	rq(i,NNy+1,1)=re1
	rq(i,NNy+1,2)=re2
	rq(i,NNy+1,3)=re3
	rq(i,NNy+1,4)=re4

 enddo
endif

select case (ExplicitSchemeType) ! выбор типа аппроксимации========================================================================================

case ('TVD') ! Колган------------------------------------------------------------------------------------------------------------------


do i=1,NNx
 do j=1,NNy
  do k=1,4

  	dqy(i,j,k)=minmod(rq(i,j,k),rq(i,j+1,k))

  enddo

  	c2=(gam*q(i,j,4)/q(i,j,1))**2
	dqy(i,j,1)=(dqy(i,j,4)-dqy(i,j,1))/c2

 enddo
enddo



case ('ENO') ! ENO-схема---------------------------------------------------------------------------------------------------------------

! вторые разности (производные)
 do i=1,NNx
  
  if (i<=N1.or.i>N2) then ! граница периодичности

  	ko1=  4./(dvert(i,2)+2.*dvert(i,1)+dvert(i,NNy))*(dvert(i,1)+dvert(i,NNy))/2.
	koNy1=4./(dvert(i,1)+2.*dvert(i,NNy)+dvert(i,NNy-1))*(dvert(i,NNy)+dvert(i,NNy-1))/2.
	  
   do k=1,4
	
	rrq(i,1,k)=ko1*(rq(i,2,k)-rq(i,1,k))
	rrq(i,NNy,k)=koNy1*(rq(i,NNy+1,k)-rq(i,NNy,k))

!	rrq(i,NNy+1,k)=4./(dvert(i,1)+2.*dvert(i,NNy)+dvert(i,NNy-1))*(rq(i,NNy+1,k)-rq(i,NNy,k))
	rrq(i,0,k)=rrq(i,NNy,k)
	rrq(i,NNy+1,k)=rrq(i,1,k)

   enddo

  else       ! стенка

	ko1=  4./(dvert(i,2)+2.*dvert(i,1))*dvert(i,1)/2.
	koNy1=4./(2.*dvert(i,NNy)+dvert(i,NNy-1))*(dvert(i,NNy)+dvert(i,NNy-1))/2.

   do k=1,4
	
	rrq(i,1,k)=ko1*(rq(i,2,k)-rq(i,1,k))
	rrq(i,0,k)=rrq(i,1,k)

	rrq(i,NNy,k)=koNy1*(rq(i,NNy+1,k)-rq(i,NNy,k))
	rrq(i,NNy+1,k)=rrq(i,NNy,k)

   enddo
  
  endif

  do j=2,NNy-1

	ko=4./(dvert(i,j+1)+2.*dvert(i,j)+dvert(i,j-1))*(dvert(i,j)+dvert(i,j-1))/2.
	do k=1,4

		rrq(i,j,k)=ko*(rq(i,j+1,k)-rq(i,j,k))

	enddo
  enddo

 enddo
 
! вычисление производных --------------------------------------------------------------------------------------------------------

 do i=1,NNx
  do j=1,NNy
   do k=1,4

	dqy(i,j,k)=minmod(rq(i,j,k)+1./2.*minmod(rrq(i,j,k),rrq(i,j-1,k)),rq(i,j+1,k)-1./2.*minmod(rrq(i,j+1,k),rrq(i,j,k)))
	
   enddo
  enddo
 enddo


case ('Zij')  ! ISNAS Зийлемы------------------------------------------------------------------------------------------------------


do i=1,NNx
 do j=1,NNy
  do k=1,4

	dqy_minus(i,j,k)=isnas(rq(i,j,k),rq(i,j+1,k))*rq(i,j+1,k)

  enddo
 enddo
enddo


do i=1,NNx
 do j=1,NNy
  do k=1,4

	dqy_plus(i,j,k)=isnas(rq(i,j,k),rq(i,j+1,k))*rq(i,j,k)
   
  enddo
 enddo
enddo

end select !==================================================================================================================

! если задана минимизация характеристических производных, переход обратно к примитивным переменным

if (charmin==1) then
 do i=1,NNx
  do j=1,NNy

	call multLinv(re1,re2,re3,re4,dqy(i,j,1),dqy(i,j,2),dqy(i,j,3),dqy(i,j,4),q(i,j,1),q(i,j,4),yksic(i,j),xksic(i,j),'eta')
	
	dqy(i,j,1)=re1
	dqy(i,j,2)=re2
	dqy(i,j,3)=re3
	dqy(i,j,4)=re4
	
   ! "правая" и "левая" производные ISNAS

	call multLinv(re1,re2,re3,re4,dqy_minus(i,j,1),dqy_minus(i,j,2),dqy_minus(i,j,3),dqy_minus(i,j,4),q(i,j,1),q(i,j,4),yksic(i,j),xksic(i,j),'eta')
   
	dqy_minus(i,j,1)=re1
	dqy_minus(i,j,2)=re2
	dqy_minus(i,j,3)=re3
	dqy_minus(i,j,4)=re4

	call multLinv(re1,re2,re3,re4,dqy_plus(i,j,1),dqy_plus(i,j,2),dqy_plus(i,j,3),dqy_plus(i,j,4),q(i,j,1),q(i,j,4),yksic(i,j),xksic(i,j),'eta')
	
	dqy_plus(i,j,1)=re1
	dqy_plus(i,j,2)=re2
	dqy_plus(i,j,3)=re3
	dqy_plus(i,j,4)=re4

  enddo
 enddo
endif

end if

! вычисление производных по времени==================================================================================================================

if ((ExplicitSchemeType=='God').or.(ImplicitSchemeType.ne.'NONE')) then
	
	 dqty=0.
else	


do i=1,NNx
 do j=1,NNy

  if (ImplicitSchemeType=='NONE') then

   if (ExplicitSchemeType.ne.'Zij') then

      call multAeta(dqty(i,j,1),dqty(i,j,2),dqty(i,j,3),dqty(i,j,4),dqy(i,j,1),dqy(i,j,2),dqy(i,j,3),dqy(i,j,4),q(i,j,1),q(i,j,2),q(i,j,3),q(i,j,4),yksic(i,j),xksic(i,j))
    
	  do k=1,4
        dqty(i,j,k)=-dqty(i,j,k) 
      enddo
	
	else if (ExplicitSchemeType=='Zij') then
	
	  re1=(dqy_minus(i,j,1)+dqy_plus(i,j,1))/2.
	  re2=(dqy_minus(i,j,2)+dqy_plus(i,j,2))/2.
	  re3=(dqy_minus(i,j,3)+dqy_plus(i,j,3))/2.
	  re4=(dqy_minus(i,j,4)+dqy_plus(i,j,4))/2.
	
	  call multAeta(dqty(i,j,1),dqty(i,j,2),dqty(i,j,3),dqty(i,j,4),re1,re2,re3,re4,q(i,j,1),q(i,j,2),q(i,j,3),q(i,j,4),yksic(i,j),xksic(i,j))
    
	  do k=1,4
        dqty(i,j,k)=-dqty(i,j,k) 
      enddo
	
	endif
  
  else

  do k=1,4

   dqty(i,j,k)=0.  ! если схема неявная, берем равными нулю
     
  enddo

  endif

 enddo
enddo
endif

dqt=dqtx+dqty


end subroutine derivatives
