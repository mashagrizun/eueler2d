! граничные услови€ на входе 

subroutine conditionsInlet
use commonArrays 
use commonVariables

implicit none

integer :: j

real*8 :: i0,di,s0,ds,du,dv,dp,dr,dalpha,alpha0,c0,p_dim,r_dim
real*8 :: un1,ut1,un2,ut2
real*8 :: q11,q12,q41,q42
real*8 :: left(4)=0.,right(4)=0.


if (flowType == 0) then ! дозвуковой поток

	do j=1,NNy
 
	 if (abs(q(1,j,3))<=del.and.abs(q(1,j,2))<=del)then    
	  alpha0=0.
	 else
	  alpha0=datan2(q(1,j,3),q(1,j,2))        ! угол в ближайшей к границе €чееке
	 endif
 
	 dalpha=alpha_in-alpha0

	 i0=gam/gamm*q(1,j,4)/q(1,j,1)+(q(1,j,2)**2+q(1,j,3)**2)/2. ! энтальпи€ в ближайщей к границе €чейке
	 di=i_in-i0


	 s0=q(1,j,4)/q(1,j,1)**gam              ! энтропийна€ функци€ в ближайшей к границе €чейке
	 ds=s_in-s0

	 c0=dsqrt(gam*q(1,j,4)/q(1,j,1))        ! скорость звука в ближайшей к границе €чейке


	 du=(di-q(1,j,2)*q(1,j,3)/dcos(alpha_in)**2*dalpha-gammo*q(1,j,1)**gamm*ds)/(c0+q(1,j,2)+q(1,j,3)*dtan(alpha_in))

	 dv=dtan(alpha_in)*du+q(1,j,2)/dcos(alpha_in)**2*dalpha

	! при вычислении du, dv вместо alpha0 (как это записанно в формулах) вз€то alpha_in (не на всех сетках устанавливаетс€ решение)

	 dp=q(1,j,1)*c0*du

	 dr=(q(1,j,1)*c0*du-q(1,j,1)**gam*ds)/c0**2

 	 q11= q(1,j,1)+dr  
     un1=(q(1,j,2)+du)*yeta(1,j)-(q(1,j,3)+dv)*xeta(1,j) 
     ut1=(q(1,j,2)+du)*xeta(1,j)+(q(1,j,3)+dv)*yeta(1,j) 
     q41= q(1,j,4)+dp  

	 
		 if (Ntime==NtimeMax) then

!		  open (1000,file='p1D.dat')
!		  write (1000,*) q41

!		  open (2000,file='r1D.dat')
!		  write (2000,*) q11

!		  open (3000,file='u11D.dat')
!		  write (3000,*) un1

		endif


	 q12= q(1,j,1)
     un2= q(1,j,2)*yeta(1,j)-q(1,j,3)*xeta(1,j) 
     ut2= q(1,j,2)*xeta(1,j)+q(1,j,3)*yeta(1,j) 
     q42= q(1,j,4)  

	 	 if (Ntime==NtimeMax) then

!		  open (7000,file='p2D.dat')
!		  write (7000,*) q42

!		  open (8000,file='r2D.dat')
!		  write (8000,*) q12

!		  open (9000,file='u22D.dat')
!		  write (9000,*) un2

		endif

	 call RR1d(q11,un1,ut1,q41,q12,un2,ut2,q42)

		 if (Ntime==NtimeMax) then

!		  open (14000,file='u11DRR.dat')
!		  write (14000,*) un1

!		  open (15000,file='r1DRR.dat')
!		  write (15000,*) q11

!		  open (13000,file='p1DRR.dat')
!		  write (13000,*) q41

		endif
	 
     qx(1,j,1)= q11
     qx(1,j,2)= un1*yeta(1,j)+ut1*xeta(1,j)   ! обратное преобразование скоростей
     qx(1,j,3)=-un1*xeta(1,j)+ut1*yeta(1,j)
     qx(1,j,4)= q41

	enddo

else if (flowType == 1) then  ! сверхзвуковой поток
	
	do j=1,NNy

		Mach=dsqrt(2./gamm*((p_in/p_out)**(gamm/gam)-1.))
		
		qx(1,j,1)=p_out/rg/(t_in/(1.+gamm/2.*Mach**2))    ! плотность

		qx(1,j,4)=p_out									  ! давление

		qx(1,j,2)=Mach*dsqrt(gam*p_out/qx(1,j,1))*dcos(alpha_in) ! компонента u

		qx(1,j,3)=Mach*dsqrt(gam*p_out/qx(1,j,1))*dsin(alpha_in) ! компонента v


	enddo

endif

do j=1,NNy
	 unx(1,j)=qx(1,j,2)*yeta(1,j)-qx(1,j,3)*xeta(1,j)  ! нормальна€ компонента скорости
enddo

! энтропийна€ функци€

do j=1,NNy
  s_dim(1,j)=qx(1,j,4)/qx(1,j,1)**gam
enddo

!if (Ntime==1.or.mod(Ntime,500)==0) then
! open (19, file = 'check_entropy.dat',status='old')
! write (19,'(a15,i5)') 'Time iteration=',Ntime
! write(19,*) 's_dim'
! do j=1,NNy
!  write(19,*) s_dim(1,j)
! enddo
! write(19,*)
!endif

end subroutine conditionsInlet
