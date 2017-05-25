! выбор шага по времени

subroutine timeStep
use commonArrays
use commonVariables

implicit none

integer :: i,j
real*8 :: cc,taux,tauy,cfl_check,uk,vk,t,jjx,jjy,metrx,metry
real*4 :: tautau
t=0.

open (6, file = 'check_time.dat')

!do i=1,NNx
! do j=1,NNy
  

!  call crash(0,q(i,j,1),q(i,j,4),i,j)

!  cc=dsqrt(gam*q(i,j,4)/q(i,j,1))
!  taux=cfl*h_x(i,j)/max(abs(q(i,j,2)-cc),abs(q(i,j,2)),abs(q(i,j,2)+cc))
!  tauy=cfl*deta(i,j)/max(abs(q(i,j,3)-cc),abs(q(i,j,3)),abs(q(i,j,3)+cc))

!  tau(i,j)=min(taux,tauy)
!!  write (6,*) i,j, tau(i,j)
! enddo
!enddo

!tau_min=tau(1,1)

!do i=1,NNx
! do j=1,NNy
!  if (tau(i,j)<=tau_min) then
!   tau_min=tau(i,j)
!  endif
! enddo
!enddo

open(211000,file='tau.dat')

tau_min=1.e+10

do i=1,NNx
 do j=1,NNy
 
  uk=(yeta(i,j)*q(i,j,2)-xeta(i,j)*q(i,j,3))
  vk=(-yksi(i,j)*q(i,j,2)+xksi(i,j)*q(i,j,3))
 
!  write (6,'(a3,d13.7)') 'uk=',uk
!  write (6,'(i2,a1,i2,a3,d13.7)') i,' ',j,'vk=',vk
  
  cc=dsqrt(gam*q(i,j,4)/q(i,j,1))
  
   if (i==1.or.i==NNx) then
		jjx=Jac(i,j)
   else
		jjx=(Jac(i-1,j)*dksi(i,j)/2.+Jac(i,j)*dksi(i-1,j)/2.)/(dksi(i-1,j)/2.+dksi(i,j)/2.)
   endif


   if (j==1.or.j==NNy) then
		jjy=Jac(i,j)
   else
		jjy=(Jac(i,j-1)*deta(i,j)/2.+Jac(i,j)*deta(i,j-1)/2.)/(deta(i,j-1)/2.+deta(i,j)/2.)
   endif

!  jjy=xksi(i,j+1)*metry-metrx*yksi(i,j+1)
!  jjy=xksic(i,j)*metry-metrx*yksic(i,j)

!    write(6,*) 'yksic(i,j)',yksic(i,j)
!  write(6,*) 'xksic(i,j)',xksic(i,j)

!   write(6,*) i,j,1./jjy

 ! write (6,'(i2,a1,i2,a5,d13.7)') i,' ',j,' cc= ',cc
!  write(6,*)

  taux=cfl*dksi(i,j)/2./(abs(uk)+cc)*jjx
  tauy=cfl*deta(i,j)/2./(abs(vk)+cc)*jjy

!  read (211000,*) tautau

!  write (6,'(i2,a1,i2,a1,d13.7,a1,d13.7)') i,' ',j,' ', tautau,' ', taux
!  write(6,*)
!  if (abs(tautau-taux)>1.e-7) then !write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!	write (6,'(i2,a1,i2,a1,d13.7,a1,d13.7)') i,' ',j,' ', abs(tautau-taux),' ', taux/tautau
!	write (6,'(i2,a1,i2,a1,d13.7,a1,d13.7)') i,' ',j,' ', tautau,' ', taux
!	write(6,*)
!  endif


 ! write (6,'(a5,d13.7)') 'taux=',taux

!  write(6,'(a4,d13.7)') 'xksi',xksi(i,j) 
!  write(6,'(a4,d13.7)') 'yksi',yksi(i,j) 

!    write (6,'(a5,d13.7)') 'tauy=',tauy
!  write(6,*)

!  		  write(6,*) 'nku',cfl 
!		  write(6,'(a4,d13.7)') 'hdz ',dksi(i,j)/2. 
!		  write(6,'(a4,d13.7)') 'het ',deta(i,j)/2. 
	
!			write(6,*)
!		  write(6,'(a10,d13.7)') 'dz0(ijkp) ',1./jjx	!dsqrt(yeta(i,j)**2+xeta(i,j)**2) 
!		  write(6,'(a10,d13.7)') 'et0(ijkp) ',1./jjy	!dsqrt(xksi(i,j)**2+yksi(i,j)**2) 
!  write(6,*)

  if (taux<tau_min) tau_min=taux
  if (tauy<tau_min) tau_min=tauy

  tau(i,j)=min(taux,tauy)
!  tau(i,j)=tau_min

 ! write (6,*) 'tau',i,j,tau(i,j)

!   read(211000,*) tautau

!  if (abs(tautau-tauy)>1.e-7) write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!	write (6,'(i2,a1,i2,a1,d13.7,a1,d13.7)') i,' ',j,' ', abs(tautau-tauy),' ', tauy/tautau
!	write (6,'(i2,a1,i2,a1,d13.7,a1,d13.7)') i,' ',j,' ', tautau,' ', tauy
!	write(6,*)
!  endif
	

 enddo
enddo
!write(*,*) 't=',t

if (ExplicitSchemeType.ne.'God'.and.ImplicitSchemeType=='NONE') then 
	tau_min=0.5*tau_min
	tau=0.5*tau
endif

do i=1,NNx
 do j=1,NNy
!     read(211000,*) tautau

!  if (abs(tautau-tau(i,j))>1.e-7) write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!	write (6,'(i2,a1,i2,a1,d13.7,a1,d13.7)') i,' ',j,' ', abs(tautau-tau(i,j)),' ', tau(i,j)/tautau
!	write (6,'(i2,a1,i2,a1,d13.7,a1,d13.7)') i,' ',j,' ', tautau,' ', tau(i,j)
!	write(6,*)
!  endif
enddo
enddo
do i=1,NNx
 do j=1,NNy+1
!     read(211000,*) tautau
! write (6,*) xksi(i,j),tautau
 enddo
enddo

!  write (6,*) 'tau_min', tau_min
!do i=1,NNx
! do j=1,NNy

!  cc=dsqrt(gam*q(i,j,4)/q(i,j,1))
!  cfl_check=cc*tau_min/h_x(i,j)

!  write (6,*) i,j, cfl_check
! enddo
!enddo

! close (6)
end subroutine timeStep