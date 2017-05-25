! граничные услови€ на выходе 

subroutine conditionsOutlet
use commonArrays
use commonVariables

implicit none
real*8 :: dp,du,dv,dr,c0
integer :: j
real*8 :: un1,ut1,un2,ut2
real*8 :: q11,q12,q41,q42
real*8 :: left(4)=0.,right(4)=0.

do j=1,NNy 

	c0=dsqrt(gam*q(NNx,j,4)/q(NNx,j,1))  ! скорость звука в ближайщей к границе €чейке

!	if ((flowType == 0).or.(q(NNx,j,2)*dcos(datan2(q(NNx,j,3),q(NNx,j,2)))<dsqrt(gam*q(NNx,j,1)/q(NNx,j,4)))) then 
	
	if ((q(NNx,j,2)*yeta(NNx+1,j)-q(NNx,j,3)*xeta(NNx+1,j))<c0) then 

 
		dp=p_out-q(NNx,j,4)                
		du=-dp/q(NNx,j,1)/c0
		dv=0.
		dr=dp/c0**2

 
	 q11=q(NNx,j,1)
     un1=q(NNx,j,2)*yeta(NNx+1,j)-q(NNx,j,3)*xeta(NNx+1,j) 
     ut1=q(NNx,j,2)*xeta(NNx+1,j)+q(NNx,j,3)*yeta(NNx+1,j) 
     q41=q(NNx,j,4)  

	 !q12= q(NNx,j,1)+dr  
     !un2=(q(NNx,j,2)+du)*yeta(NNx+1,j)-(q(NNx,j,3)+dv)*xeta(NNx+1,j) 
     !ut2=(q(NNx,j,2)+du)*xeta(NNx+1,j)+(q(NNx,j,3)+dv)*yeta(NNx+1,j) 
     !q42= q(NNx,j,4)+dp 

	 q12=q11
     un2=un1
     ut2=ut1
     q42=p_out

	 call RR1d(q11,un1,ut1,q41,q12,un2,ut2,q42)
	 
     qx(NNx+1,j,1)= q11
     qx(NNx+1,j,2)= un1*yeta(NNx+1,j)+ut1*xeta(NNx+1,j)   ! обратное преобразование скоростей
     qx(NNx+1,j,3)=-un1*xeta(NNx+1,j)+ut1*yeta(NNx+1,j)
     qx(NNx+1,j,4)= q41
	
	else if ((q(NNx,j,2)*yeta(NNx+1,j)-q(NNx,j,3)*xeta(NNx+1,j))>=c0) then

		qx(NNx+1,j,1)=q(NNx,j,1)
		qx(NNx+1,j,2)=q(NNx,j,2)
		qx(NNx+1,j,3)=q(NNx,j,3)
		qx(NNx+1,j,4)=q(NNx,j,4)

	endif
enddo

do j=1,NNy
	unx(NNx+1,j)=qx(NNx+1,j,2)*yeta(NNx+1,j)-qx(NNx+1,j,3)*xeta(NNx+1,j) ! нормальна€ компонента скорости
enddo


do j=1,NNy
  s_dim(NNx+1,j)=qx(NNx+1,j,4)/qx(NNx+1,j,1)**gam
enddo

if (Ntime==1.or.mod(Ntime,100)==0) then
 open (19, file = 'check_entropy.dat',status='old')
 write (19,'(a15,i5)') 'Time iteration=',Ntime
 write(19,*) 's_dim'
 do j=1,NNy
  write(19,*) s_dim(NNx+1,j)
 enddo
 write(19,*)
endif


end subroutine conditionsOutlet