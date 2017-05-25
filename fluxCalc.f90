! вычисление потоков

subroutine fluxCalc

use commonArrays
use commonVariables

implicit none

integer :: i,j,k

real*8 :: un1,ut1,un2,ut2
real*8 :: q11,q12,q41,q42
real*8 :: left(4)=0.,right(4)=0.,qnp(3,200)
real*4 :: p1,r1,u11,p2,r2,u22,u11RR,r1RR, p1RR
real*8 :: un1RR,q11RR,q41RR,hdzl

select case(RRTypeChoice)

case('classic')!************************************************************************************************************

 ! задача распада разрыва решаетс€ через грань по каждому направлению


 do i=2,NNx
  do j=1,NNy

    do k=1,4
 
	  if (ExplicitSchemeType=='Zij') then    
		left(k)=(dqx_minus(i-1,j,k)*dh(i-1,j)+dqt(i-1,j,k)*tau(i-1,j))/2.
		right(k)=(-dqx_plus(i,j,k)*dh(i,j)+dqt(i,j,k)*tau(i,j))/2.
	  else
		left(k)=(dqx(i-1,j,k)*dh(i-1,j)+dqt(i-1,j,k)*tau(i-1,j))/2.
		right(k)=(-dqx(i,j,k)*dh(i,j)+dqt(i,j,k)*tau(i,j))/2.  
	  endif

	enddo

    q11= q(i-1,j,1)+left(1)  
    un1=(q(i-1,j,2)+left(2))*yeta(i,j)-(q(i-1,j,3)+left(3))*xeta(i,j) ! нормальна€ составл€юща€ скорости в "левой" €чейке
    ut1=(q(i-1,j,2)+left(2))*xeta(i,j)+(q(i-1,j,3)+left(3))*yeta(i,j) ! касательна€ составл€юща€ скорости в "левой" €чейке
    q41= q(i-1,j,4)+left(4)  

    q12= q(i,j,1)+right(1)  
    un2=(q(i,j,2)+right(2))*yeta(i,j)-(q(i,j,3)+right(3))*xeta(i,j) ! нормальна€ составл€юща€ скорости в "правой" €чейке
    ut2=(q(i,j,2)+right(2))*xeta(i,j)+(q(i,j,3)+right(3))*yeta(i,j) ! касательна€ составл€юща€ скорости в "правой" €чейке
    q42= q(i,j,4)+right(4)  

		 if (Ntime==NtimeMax) then
!	if (i==25) then

!		  open (1000,file='p1D.dat')
!		  write (1000,'(d13.7)') dqx(i-1,j,1)


!		  open (2000,file='r1D.dat')
!		  write (2000,'(d13.7)') dqx(i-1,j,2)


!		  open (3000,file='u11D.dat')
!		  write (3000,'(d13.7)') dqx(i-1,j,3)
!endif
	endif

	call crash(1,q11,q41,i-1,j)
	call crash(1,q12,q42,i,j)

		 if (Ntime==NtimeMax) then
!if (i==25) then
!		  open (7000,file='p2D.dat')
!		  write (7000,'(d13.7)') dqx(i-1,j,4)

!		  open (8000,file='r2D.dat')
!		  write (8000,'(d13.7)') dqx(i,j,1)

!		  open (9000,file='u22D.dat')
!		  write (9000,'(d13.7)') dqx(i,j,2)
!endif
	endif

	call RR1d(q11,un1,ut1,q41,q12,un2,ut2,q42)

		 if (Ntime==NtimeMax) then
!if (i==25) then
!		  open (14000,file='u11DRR.dat')
!		  write (14000,'(d13.7)') dqx(i,j,3)

!		  open (15000,file='r1DRR.dat')
!		  write (15000,'(d13.7)') dqx(i,j,4)

!		  open (13000,file='p1DRR.dat')
!		  write (13000,'(d13.7)') q41
!endif
	endif

    qx(i,j,1)= q11
    qx(i,j,2)= un1*yeta(i,j)+ut1*xeta(i,j)   ! обратное преобразование скоростей
    qx(i,j,3)=-un1*xeta(i,j)+ut1*yeta(i,j)
    qx(i,j,4)= q41

!	 if (Ntime==NtimeMax) then

!		  open (14000,file='u11DRR.dat')
!		  write (14000,'(d13.7)') 0.0

!		  open (15000,file='r1DRR.dat')
!		  write (15000,'(d13.7)') qx(i,j,3)

!		  open (13000,file='p1DRR.dat')
!		  write (13000,'(d13.7)') qx(i,j,2)

!	endif

    unx(i,j)=qx(i,j,2)*yeta(i,j)-qx(i,j,3)*xeta(i,j)

  enddo
 enddo


 do i=1,NNx
  do j=2,NNy

	do k=1,4

	  if (ExplicitSchemeType=='Zij') then
		left(k)=(dqy_minus(i,j-1,k)*dvert(i,j-1)+dqt(i,j-1,k)*tau(i,j-1))/2.
		right(k)=(-dqy_plus(i,j,k)*dvert(i,j)+dqt(i,j,k)*tau(i,j))/2.
	  else
		left(k)=(dqy(i,j-1,k)*dvert(i,j-1)+dqt(i,j-1,k)*tau(i,j-1))/2.
		right(k)=(-dqy(i,j,k)*dvert(i,j)+dqt(i,j,k)*tau(i,j))/2.
      endif
	
	enddo

    q11= q(i,j-1,1)+left(1) 
    un1=(q(i,j-1,3)+left(3))*xksi(i,j)-(q(i,j-1,2)+left(2))*yksi(i,j) ! нормальна€ составл€юща€ скорости в "нижней" €чейке
    ut1=(q(i,j-1,3)+left(3))*yksi(i,j)+(q(i,j-1,2)+left(2))*xksi(i,j) ! касательна€ составл€юща€ скорости в "нижней" €чейке
    q41= q(i,j-1,4)+left(4) 

	 if ((Ntime==NtimeMax).and.(i>N1.and.i<=N2)) then
		open (1000,file='p1D.dat')
		
		if (j==2) then
		  
!		  write (1000,'(i3,i3,a1,d13.7)')i,j-1,' ', rq(i,j-1,1)
!		  write (1000,'(i3,i3,a1,d13.7)')i,j,' ', rq(i,j,1)

		endif
		
		if (j==NNy) then
		
!		  write (1000,'(i3,i3,a1,d13.7)')i,j,' ', rq(i,j,1)
!		  write (1000,'(i3,i3,a1,d13.7)')i,j+1,' ', rq(i,j+1,1)

		endif
		!if (dqy(i,j,4)>1.0e-3) then
		  open (2000,file='r1D.dat')
		  if (j==2) then

!			write (2000,'(i3,i3,a1,d13.7)') i,j-1, ' ', rq(i,j-1,4)
!			write (2000,'(i3,i3,a1,d13.7)') i,j, ' ', rq(i,j,4)
		  endif

		  if (j==NNy) then
		  	
!			write (2000,'(i3,i3,a1,d13.7)') i,j, ' ', rq(i,j,4)

!			write (2000,'(i3,i3,a1,d13.7)') i,j+1, ' ', rq(i,j+1,4)

		  endif
		!endif

!		  open (3000,file='u11D.dat')
!		  write (3000,'(d13.7)') un1

	endif

    q12= q(i,j,1)+right(1) 
    un2=(q(i,j,3)+right(3))*xksi(i,j)-(q(i,j,2)+right(2))*yksi(i,j) ! нормальна€ составл€юща€ скорости в "верхней" €чейке
    ut2=(q(i,j,3)+right(3))*yksi(i,j)+(q(i,j,2)+right(2))*xksi(i,j) ! касательна€ составл€юша€ скорости в "верхней" €чейке
    q42= q(i,j,4)+right(4) 

	call crash(2,q11,q41,i,j-1)
	call crash(2,q12,q42,i,j)

	 if (Ntime==NtimeMax) then

!		  open (7000,file='p2D.dat')
!		  write (7000,'(d13.7)') q42

!		  open (8000,file='r2D.dat')
!		  write (8000,'(d13.7)') q12

!		  open (9000,file='u22D.dat')
!		  write (9000,'(d13.7)') un2

	endif


	call RR1d(q11,un1,ut1,q41,q12,un2,ut2,q42)

	 if (Ntime==NtimeMax) then

!		  open (14000,file='u11DRR.dat')
!		  write (14000,'(d13.7)') un1

!		  open (15000,file='r1DRR.dat')
!		  write (15000,'(d13.7)') q11

!		  open (13000,file='p1DRR.dat')
!		  write (13000,'(d13.7)') q41

	endif
  
    qy(i,j,1)= q11
    qy(i,j,2)=-un1*yksi(i,j)+ut1*xksi(i,j)  ! обратное преобразование скоростей
    qy(i,j,3)= un1*xksi(i,j)+ut1*yksi(i,j)
    qy(i,j,4)= q41

!		 if (Ntime==NtimeMax) then

!		  open (14000,file='u11DRR.dat')
!		  write (14000,'(d13.7)') 0.0

!		  open (15000,file='r1DRR.dat')
!		  write (15000,'(d13.7)') qy(i,j,3)

!		  open (13000,file='p1DRR.dat')
!		  write (13000,'(d13.7)') qy(i,j,2)

! 	endif

    uny(i,j)=qy(i,j,3)*xksi(i,j)-qy(i,j,2)*yksi(i,j)

  enddo
 enddo

! 	close(1000);close(2000);close(3000)
!	close(7000);close(8000);close(9000)
!	close(14000);close(15000);close(13000)


 ! вычисление вектора потока по направлению х

 do i=1,NNx+1
  do j=1,NNy

    call fl_X(flx(i,j,1),flx(i,j,2),flx(i,j,3),flx(i,j,4),qx(i,j,1),qx(i,j,2),qx(i,j,3),qx(i,j,4),unx(i,j),yeta(i,j),xeta(i,j))

  enddo
 enddo



 ! вычисление вектора потока по направлению у

 do i=1,NNx
  do j=1,NNy+1
  
    call fl_Y(fly(i,j,1),fly(i,j,2),fly(i,j,3),fly(i,j,4),qy(i,j,1),qy(i,j,2),qy(i,j,3),qy(i,j,4),uny(i,j),xksi(i,j),yksi(i,j))
  
  enddo
 enddo

end select

!call checkPrint('qx    ')
!call checkPrint('qy    ')

!call checkPrint('flx   ')
!call checkPrint('fly   ')


!open (1010101,file='results1.dat')
!open (1010102,file='results2.dat')
!open (1010103,file='results3.dat')
!open (1010104,file='results4.dat')
!open (1010105,file='results5.dat')
!open (1010106,file='results6.dat')
!open (1010107,file='results7.dat')
!open (1010108,file='results8.dat')
!open (1010109,file='results9.dat')

if (Ntime==NtimeMax) then
! do i=1,NNx
i=1
  do j=1,NNy
!   		  open (81000,file='p1D.dat')
!		  read (81000,*) q41

!		  open (82000,file='r1D.dat')
!		  read (82000,*) q11

!		  open (83000,file='u11D.dat')
!		  read (83000,*) un1


!		  open (87000,file='p2D.dat')
!		  read (87000,*) q42

!		  open (88000,file='r2D.dat')
!		  read (88000,*) q12

!		  open (89000,file='u22D.dat')
!		  read (89000,*) un2

		  

!		  open (814000,file='u11DRR.dat')
!		  read (814000,*) un1RR

!		  open (815000,file='r1DRR.dat')
!		  read (815000,*) q11RR

!		  open (813000,file='p1DRR.dat')
!		  read (813000,*) q41RR


		  
!		  open (810001,file='p1Dzeta.dat')
!
!		  open (820001,file='r1Dzeta.dat')

!		  open (830001,file='u11Dzeta.dat')
!
!		  read(810001,*) p1

!		  read(820001,*) r1
!
!!		  read(830001,*) u11


		  
!		  open (870001,file='p2Dzeta.dat')

!		  open (880001,file='r2Dzeta.dat')
!
!		  open (890001,file='u22Dzeta.dat')


!		  read(870001,*) p2

!		  read(880001,*) r2

!		  read(890001,*) u22


		  
!		  open (8140001,file='u11DzetaRR.dat')

!		  open (8150001,file='r1DzetaRR.dat')

!		  open (8130001,file='p1DzetaRR.dat')

 		  
!		  read(8140001,*) u11RR

!		  read(8150001,*) r1RR
!
!		  read(8130001,*) p1RR

	
    if (abs(p1-q41)>1.e-10) then
!		write(1010101,*) 'p1',i,j, p1, q41,abs(p1-q41)
	endif

	if (abs(r1-q11)>1.e-10) then
!		write(1010102,*) 'r1',i,j, r1, q11,abs(r1-q11)
	endif

    if (abs(u11-un1)>1.e-10) then
!		write(1010103,*) 'u11',i,j, u11, un1,abs(u11-un1)
	endif




    if (abs(p2-q42)>1.e-10) then
!		write(1010104,*) 'p2',i,j, p2, q42,abs(p2-q42)
	endif

	if (abs(r2-q12)>1.e-10) then
!		write(1010105,*) 'r2',i,j, r2, q12,abs(r2-q12)
	endif

    if (abs(u22-un2)>1.e-10) then
!		write(1010106,*) 'u22',i,j, u22, un2,abs(u22-un2)
	endif




    if (abs(u11RR-un1RR)>1.e-10) then
!		write(1010108,*) 'u11RR',i,j, u11RR, un1RR,abs(u11RR-un1RR)
	endif

	if (abs(r1RR-q11RR)>1.e-10) then
!		write(1010107,*)'r1RR', i,j, r1RR, q11RR, abs(r1RR-q11RR)
	endif

    if (abs(p1RR-q41RR)>1.e-10) then
!		write(1010109,*) 'p1RR',i,j, p1RR, q41RR,abs(p1RR-q41RR)
	endif





    if (abs(r2-qx(i,j,1))>1.e-7) then
!		write(1010105,*) 'r',i,j, r2, qx(i,j,1),abs(r2-qx(i,j,1))
	endif

	
    if (abs(u22-0.)>1.e-7) then
!		write(1010106,*) 'u',i,j, u22, 0.,abs(u22-0.)
	endif


	if (abs(u11RR-qx(i,j,3))>1.e-7) then
!		write(1010107,*)'v', i,j, u11RR, qx(i,j,3), abs(u11RR-qx(i,j,3))
	endif

    if (abs(r1RR-qx(i,j,2))>1.e-7) then
!		write(1010108,*) 'w',i,j, r1RR, qx(i,j,2),abs(r1RR-qx(i,j,2))
	endif


    if (abs(p1RR-qx(i,j,4))>1.e-7) then
!		write(1010109,*) 'p',i,j, p1RR, qx(i,j,4),abs(p1RR-qx(i,j,4))
	endif

  enddo
! enddo   
endif

!close(1010101)
!close(1010102)
!close(1010103)
!close(1010104)
!close(1010105)
!close(1010106)
!close(1010107)
!close(1010108)
!close(1010109)


end subroutine fluxCalc
