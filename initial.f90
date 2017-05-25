! заполнение расчетного поля начaльными данными

subroutine initial 
use commonArrays
use commonVariables

implicit none

integer :: i,j,k,ii,NNN,m

NNN=Nx*Ny*Nz

write (*,*)
write (*,*) 'Start calculating initial field.'
write (*,*)

! везде задается статическое давление, выбрана нулевая скорость потока, следовательно,
! полное давление равно статическому, плотность наsйдена из условия изоэнтропичности.

write (*,*) 'Press 1 if you want to continue, press 0 to start new calculation'
read(*,*) mark

if (mark==1) then

	open (220, file='promNcontrol.dat')
	read (220,*) Ntime,(((q(i,j,k),i=1,NNx),j=1,NNy),k=1,4)

	open (320, file='promNcontrol2.dat')
	read (320,*) (((qx(i,j,k),i=1,NNx+1),j=1,NNy),k=1,4)

	open (420, file='promNcontrol3.dat')
	read (420,*) (((qy(i,j,k),i=1,NNx),j=1,NNy+1),k=1,4)

	open (520, file='promNcontrol4.dat')
	read (520,*) (((dqn(i,j,k),i=1,NNx),j=1,NNy),k=1,4)

	!Ntime=60000

	NtimeMax=Ntime+NtimeMax

	!open (21, file = 'r1.dat', form='unformatted', access='direct', recl=4*NNN)
	!read (21, rec=1)  (((r_f(ii,k,m), k=1,Ny), ii=1,Nx), m=1,Nz) 


	!open (22, file = 'u1.dat', form='unformatted', access='direct', recl=4*NNN)
	!read (22, rec=1)  (((u_f(ii,k,m), k=1,Ny), ii=1,Nx), m=1,Nz)


	!open (23, file = 'v1.dat', form='unformatted', access='direct', recl=4*NNN)
	!read (23, rec=1)  (((v_f(ii,k,m), k=1,Ny), ii=1,Nx), m=1,Nz)


	!open (24, file = 'w1.dat', form='unformatted', access='direct', recl=4*NNN)
	!read (24, rec=1)  (((w_f(ii,k,m), k=1,Ny), ii=1,Nx), m=1,Nz)


	!open (25, file = 'p1.dat', form='unformatted', access='direct', recl=4*NNN)
	!read (25, rec=1)  (((p_f(ii,k,m), k=1,Ny), ii=1,Nx), m=1,Nz)



	!do i=1,NNx
	! do j=1,NNy

	!	 q(i,j,1)=r_f(1,j,i)
	!	 q(i,j,2)=w_f(1,j,i)
	!	 q(i,j,3)=v_f(1,j,i)
	!	 q(i,j,4)=p_f(1,j,i)
 
	 ! enddo
!	enddo

endif

write (*,*) 'Ntime=', Ntime, ' NtimeMax=', NtimeMax

if (mark==0) then

	Ntime=0

	do i=1,NNx
		do j=1,NNy
		if (flowType == 1) then
			Mach=dsqrt(2./gamm*((p_in/p_out)**(gamm/gam)-1.))
			q(i,j,1)=p_out/rg/(t_in/(1.+gamm/2.*Mach**2))
		else
			q(i,j,1)=(p_out/s_in)**gamo
		endif
			q(i,j,2)=0.
!			q(i,j,2)=0.5
			q(i,j,3)=0.
			q(i,j,4)=p_out
			s_dim(i,j)=q(i,j,4)/q(i,j,1)**gam
		enddo
	enddo

endif

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


close (220);close (320);close (420);close (520);

end subroutine initial