! ¬ычисление метрических коэффициентов, параметров экстрапол€ции

subroutine geom
use commonArrays
use commonVariables
implicit none
integer i,j,ii,k,m,mid,numb,NNx1,NNy1
integer r

real*8 kkk,a1,b1,c1,dminus
real*4 :: lks1(1,251,301),let1(2,250,301),ldz1(2,251,300)
real*4 :: ndzx(1,250,301),ndzz1(1,250,301),ndzy1(1,250,301),netx1(1,251,300),nety1(1,251,300),netz1(1,251,300)
real*4 :: ndzz(250,301),jac1(1,250,300)
!=======================================================================================

! выбор венца, проверка, существует ли выбранный венец

write (*,*) 'There is(are)',rmax,' row(s).'

10 continue

write(*,*) 'Enter number of row:'
read (*,*) r

if ((r>rmax).or.(r<=0)) then
 write (*,*) 'Row with this number does not exist.'
 goto 10
endif

!=======================================================================================

mid=Nx1/2  ! среднее сечение по трехмерному направлению икс

!==================================================================================================
!  чтение файлов с координатами сетки

open (1, file = 'dkx.dat', form='unformatted', access='direct', recl=8*NN)
read (1, rec=1)  (((xold(ii,k,m), k=1,Ny+1), ii=1,Nx+1), m=1,Nz+1) 

open (2, file = 'dky.dat', form='unformatted', access='direct', recl=8*NN)
read (2, rec=1)  (((yold(ii,k,m), k=1,Ny+1), ii=1,Nx+1), m=1,Nz+1)

open (3, file = 'dkz.dat', form='unformatted', access='direct', recl=8*NN)
read (3, rec=1)  (((zold(ii,k,m), k=1,Ny+1), ii=1,Nx+1), m=1,Nz+1)

close(1)
close(2)
close(3)

ii=mid

!==================================================================================================    
! запись координат узлов сетки                  

do k=1,Ny+1 
 do m=1,Nz+1
	xn(m,k)=zold(ii,k,m)
	yn(m,k)=yold(ii,k,m)
 enddo
enddo


! координаты узлов фиктивных €чеек (без учета выполнени€ условий отражени€ на стенке - не нужны в программе)

do i=1,NNx+1

	xn(i,0)=xn(i,1)-(xn(i,NNy+1)-xn(i,NNy))
	yn(i,0)=yn(i,1)-(yn(i,NNy+1)-yn(i,NNy))

	xn(i,NNy+2)=xn(i,NNy+1)+(xn(i,2)-xn(i,1))
	yn(i,NNy+2)=yn(i,NNy+1)+(yn(i,2)-yn(i,1))

enddo

do j=1,NNy+1

	xn(0,j)=xn(1,j)-(xn(2,j)-xn(1,j))
	yn(0,j)=yn(1,j)-(yn(2,j)-yn(1,j))

	xn(NNx+2,j)=xn(NNx+1,j)+(xn(NNx+1,j)-xn(NNx,j))
	yn(NNx+2,j)=yn(NNx+1,j)+(yn(NNx+1,j)-yn(NNx,j))

enddo

xn(0,0)=xn(0,1)-(xn(0,NNy+1)-xn(0,NNy))
yn(0,0)=yn(0,1)-(yn(0,NNy+1)-yn(0,NNy))

xn(0,NNy+2)=xn(0,NNy+1)+(xn(0,2)-xn(0,1))
yn(0,NNy+2)=yn(0,NNy+1)+(yn(0,2)-yn(0,1))

xn(NNx+2,0)=xn(NNx+2,1)-(xn(NNx+2,NNy+1)-xn(NNx+2,NNy))
yn(NNx+2,0)=yn(NNx+2,1)-(yn(NNx+2,NNy+1)-yn(NNx+2,NNy))

xn(NNx+2,NNy+2)=xn(NNx+2,NNy+1)+(xn(NNx+2,2)-xn(NNx+2,1))
yn(NNx+2,NNy+2)=yn(NNx+2,NNy+1)+(yn(NNx+2,2)-yn(NNx+2,1))


!==================================================================================================    

NNx1=NNx+1; NNy1=NNy+1

!==================================================================================================

! 'координаты центров €чеек'
 
do i=0,NNx+1                     
 do j=0,NNy+1

	xc(i,j)=(xn(i,j)+xn(i+1,j)+xn(i+1,j+1)+xn(i,j+1))/4.
	yc(i,j)=(yn(i,j)+yn(i+1,j)+yn(i+1,j+1)+yn(i,j+1))/4.

 enddo
enddo

!==================================================================================================
!'"горизонтальные" грани €чеек '

do i=1,NNx
 do j=1,NNy1

	!координаты центра грани

	xh(i,j)=(xn(i+1,j)+xn(i,j))/2.
	yh(i,j)=(yn(i+1,j)+yn(i,j))/2.

	!приращени€

	dxh(i,j)=xn(i+1,j)-xn(i,j)                      
	dyh(i,j)=yn(i+1,j)-yn(i,j)
	
	!длина грани

	dksi(i,j)=dsqrt(dxh(i,j)**2+dyh(i,j)**2) 

	!метрические коэффициенты
 
	xksi(i,j)=dxh(i,j)/dksi(i,j)
	yksi(i,j)=dyh(i,j)/dksi(i,j)

 
 enddo
enddo


!==================================================================================================
! '"вертикальные" грани €чеек'

do j=1,NNy
 do i=1,NNx1
 
	!координаты центра грани

	xv(i,j)=(xn(i,j+1)+xn(i,j))/2.
	yv(i,j)=(yn(i,j+1)+yn(i,j))/2.

	!приращени€

	dxv(i,j)=xn(i,j+1)-xn(i,j)                       
	dyv(i,j)=yn(i,j+1)-yn(i,j)
 
	!длина грани

	deta(i,j)=dsqrt(dxv(i,j)**2+dyv(i,j)**2)    
 
	!метрические коэффициенты
 
	xeta(i,j)=dxv(i,j)/deta(i,j)
	yeta(i,j)=dyv(i,j)/deta(i,j)
 
 enddo
enddo

!==================================================================================================

do i=1,NNx
 do j=1,NNy
  
	!рассто€ние от центра до середины грани €чейки 

	dh_half_lt(i,j)=dsqrt((xv(i,j)-xc(i,j))**2+(yv(i,j)-yc(i,j))**2)
	dh_half_rt(i,j)=dsqrt((xv(i+1,j)-xc(i,j))**2+(yv(i+1,j)-yc(i,j))**2)
  
	dv_half_dn(i,j)=dsqrt((xh(i,j)-xc(i,j))**2+(yh(i,j)-yc(i,j))**2)
	dv_half_up(i,j)=dsqrt((xh(i,j+1)-xc(i,j))**2+(yh(i,j+1)-yc(i,j))**2)

	!рассто€ние между серединами противолежащих граней
  
	dh(i,j)=dsqrt((xv(i+1,j)-xv(i,j))**2+(yv(i+1,j)-yv(i,j))**2)
	dvert(i,j)=dsqrt((xh(i,j+1)-xh(i,j))**2+(yh(i,j+1)-yh(i,j))**2)

	!приращени€ по средним лини€м

	dxh_mid(i,j)=xv(i+1,j)-xv(i,j)
	dyh_mid(i,j)=yv(i+1,j)-yv(i,j)

	dxv_mid(i,j)=xh(i,j+1)-xh(i,j)
	dyv_mid(i,j)=yh(i,j+1)-yh(i,j)

	! "центральные" метрические коэффициенты

	xksic(i,j)=dxh_mid(i,j)/dh(i,j)
	yksic(i,j)=dyh_mid(i,j)/dh(i,j)

	xetac(i,j)=dxv_mid(i,j)/dvert(i,j)
	yetac(i,j)=dyv_mid(i,j)/dvert(i,j)

 enddo
enddo

! "центральные" метрические коэффициенты фиктивных €чеек

do j=1,NNy 

  xetac(0,j)=xeta(1,j)
  yetac(0,j)=yeta(1,j)  

  xetac(NNx+1,j)=xeta(NNx+1,j)   
  yetac(NNx+1,j)=yeta(NNx+1,j)
   
enddo

do i=1,NNx

  xksic(i,0)=xksi(i,1)
  yksic(i,0)=yksi(i,1)
  
  xksic(i,NNy+1)=xksi(i,NNy+1)
  yksic(i,NNy+1)=yksi(i,NNy+1)
  
enddo

!==================================================================================================
! 'рассто€ние между центрами €чеек'

do j=1,NNy

  dnoh(1,j)=dh_half_lt(1,j)      ! в первой €чейке берем рассто€ние от центра грани до центра €чейки
  dnoh(NNx+1,j)=dh_half_rt(NNx,j)! в последней €чейке - аналогично

  do i=2,NNx
	 dnoh(i,j)=dsqrt((xc(i-1,j)-xc(i,j))**2+(yc(i-1,j)-yc(i,j))**2) ! по-горизонтали
  enddo

enddo

do i=1,NNx

 dnov(i,1)=dv_half_dn(i,1)	    ! в первой €чейке берем рассто€ние от центра грани до центра €чейки
 dnov(i,NNy+1)=dv_half_up(i,NNy)! в последней - аналогично

 do j=2,NNy 
	dnov(i,j)=dsqrt((xc(i,j-1)-xc(i,j))**2+(yc(i,j-1)-yc(i,j))**2) ! по-вертикали
 enddo

enddo


!==================================================================================================
!'площадь €чейки (как 1/2 векторного произведени€ диагоналей)'

do i=1,NNx
 do j=1,NNy

	S(i,j)=abs((xn(i+1,j+1)-xn(i,j))*(yn(i,j+1)-yn(i+1,j))-(yn(i+1,j+1)-yn(i,j))*(xn(i,j+1)-xn(i+1,j)))/2.

 enddo
enddo

!==================================================================================================
! 'якобиан в центре €чейки'

do i=1,NNx
 do j=1,NNy

	Jac(i,j)=xksic(i,j)*yetac(i,j)-xetac(i,j)*yksic(i,j)

 enddo
enddo

!==================================================================================================
! приращени€ по х дл€ вычислени€ шага по времени (2*рассто€ние от центра до "вертикальной" стороны €чейки)

do i=1,NNx
 do j=1,NNy

  a1=yn(i,j+1)-yn(i,j)
  b1=xn(i,j)-xn(i,j+1)
  c1=xn(i,j+1)*yn(i,j)-xn(i,j)*yn(i,j+1)

  dminus=abs(a1*xc(i,j)+b1*yc(i,j)+c1)/dsqrt(a1**2+b1**2)

  h_x(i,j)=2.*dminus
  
 enddo
enddo

!open (11000, file = 'lks1.dat', form='unformatted', access='direct', recl=4*Nx*Ny1*Nz1)
!read (11000, rec=1)  (((lks1(ii,k,m), k=1,Ny1), ii=1,Nx), m=1,Nz1) 

!open (33000, file = 'let1.dat', form='unformatted', access='direct', recl=4*Nx1*Ny*Nz1)
!read (33000, rec=1)  (((let1(ii,k,m), k=1,Ny), ii=1,Nx1), m=1,Nz1)

!open (44000, file = 'ldz1.dat', form='unformatted', access='direct', recl=4*Nx1*Ny1*Nz)
!read (44000, rec=1)  (((ldz1(ii,k,m), k=1,Ny1), ii=1,Nx1), m=1,Nz)


!open (1000, file = 'jac1.dat', form='unformatted', access='direct', recl=4*Nx*Ny*Nz)
!read (1000, rec=1)  (((jac1(ii,k,m), k=1,Ny), ii=1,Nx), m=1,Nz) 

!open (2000, file = 'ndzy1.dat', form='unformatted', access='direct', recl=4*Nx*Ny*Nz1)
!read (2000, rec=1)  (((ndzy1(ii,k,m), k=1,Ny), ii=1,Nx), m=1,Nz1) 

!open (3000, file = 'ndzz1.dat', form='unformatted', access='direct', recl=4*Nx*Ny*Nz1)
!read (3000, rec=1)  (((ndzz1(ii,k,m), k=1,Ny), ii=1,Nx), m=1,Nz1) 


!open (4000, file = 'netx1.dat', form='unformatted', access='direct', recl=4*Nx*Ny1*Nz)
!read (4000, rec=1)  (((netx1(ii,k,m), k=1,Ny1), ii=1,Nx), m=1,Nz)

!open (5000, file = 'nety1.dat', form='unformatted', access='direct', recl=4*Nx*Ny1*Nz)
!read (5000, rec=1)  (((nety1(ii,k,m), k=1,Ny1), ii=1,Nx), m=1,Nz)

!open (6000, file = 'netz1.dat', form='unformatted', access='direct', recl=4*Nx*Ny1*Nz)
!read (6000, rec=1)  (((netz1(ii,k,m), k=1,Ny1), ii=1,Nx), m=1,Nz)


!open (7000, file = 'nksx1.dat', form='unformatted', access='direct', recl=4*Nx1*Ny*Nz)
!read (7000, rec=1)  (((nksx1(ii,k,m), k=1,Ny), ii=1,Nx1), m=1,Nz) 

!open (8000, file = 'nksy1.dat', form='unformatted', access='direct', recl=4*Nx1*Ny*Nz)
!read (8000, rec=1)  (((nksy1(ii,k,m), k=1,Ny), ii=1,Nx1), m=1,Nz) 

!open (9000, file = 'nksz1.dat', form='unformatted', access='direct', recl=4*Nx1*Ny*Nz)
!read (9000, rec=1)  (((nksz1(ii,k,m), k=1,Ny), ii=1,Nx1), m=1,Nz) 



!open (1000, file = 'tau.dat')
!read (1000, *)  (((ndzz1(ii,k,m), k=1,Ny), ii=1,Nx), m=1,Nz+1) 

!close(11000)
!close(33000)
!close(44000)

!close(1000)
!close(2000)
!close(3000)
!close(4000)

!open (1001,file='ndzzRead.dat')
!write (1001,*) jac1
!do i=1,NNx
! do j=1,NNy
!	write(1001,*) i,j, Jac(i,j)
! enddo
!enddo
!ii=0
!do i=1,NNx
! do j=1,NNy
!i=51
!	if (abs(ndzz1(1,j,i)-yeta(i,j))>1.e-4) then
!		write(*,*) i,j, ndzz1(1,j,i), yeta(i,j)/(xksi(i,j)*yeta(i,j)-xeta(i,j)*yksi(i,j))
!	endif
!	if (abs(ndzy1(1,j,i)+xeta(i,j))>1.e-4) then
!		write(*,*) i,j, ndzy1(1,j,i),-xeta(i,j)/(xksi(i,j)*yeta(i,j)-xeta(i,j)*yksi(i,j))
!	endif
!	if (abs(nety1(1,j,i)-xksi(i,j))>1.e-4) then
!		write(*,*) i,j, nety1(1,j,i),xksi(i,j)
!	endif
!	if (abs(netz1(1,j,i)+yksi(i,j))>1.e-4) then
!		write(*,*) i,j, netz1(1,j,i),-yksi(i,j)
!	endif

!	if (abs(jac1(1,j,i)-Jac(i,j))>1.e-4) then
!	ii=ii+1
!		write(*,*) i,j, jac1(1,j,i), Jac(i,j) !, abs(jac1(1,j,i)-Jac(i,j))
!	endif

 !enddo
!enddo
!write(*,*)ii

end subroutine geom
