! ������ ����������

subroutine readFiles
use commonVariables

implicit none

character*10 :: a
integer :: i,j,ii,k,m,r2

del=1.e-5
eps=0.05
pi=4.*atan(1.)

label=1

write (*,*)
write (*,*) 'Start reading files.'
write (*,*)

open (4, file = '3d_geom.dat', form='formatted', access='direct', recl=80)
read (4,*)
read (4,*) rmax
close(4)

!r2=r*2-1

 open (4, file = '3d_geom.dat', form='formatted', access='direct', recl=80)
 read (4,*)
! do ir=1,r2
!  open (4, file = '3d_geom.dat', form='formatted', access='direct', recl=80)
!  read (4,*) 
  read (4,*) 
! enddo
 read (4,*)  Nx, Ny, Nz, N1, N2
 close (4)

write(*,*)
write(*,'(a18,a5,x,i4,x,a5,x,i4,x,a5,x,i4,x,a4,x,i4,x,a4,x,i4)') 'Grid Parameters: ', 'Nx =', Nx, 'Ny =', Ny, 'Nz =', Nz, 'N1 =', N1, 'N2 =',  N2
write(*,*)

open (221,file='crash.dat')

Nx1=Nx+1; Ny1=Ny+1; Nz1=Nz+1
NN=Nx1*Ny1*Nz1
NNx=Nz; NNy=Ny

format1='(2x,a6,i1)'
format2='(a2,4x,200i13,2x)'
!format3='(a2,x,i3,x,200f10.5,2x)'
format3='(a2,x,i3,x,200d13.5,2x)'

open (11, file='initial.dat', form='formatted', status='old')

read (11,*) a, gam					! ���������� ��������
read (11,*) a, rg					! ������� ���������� 
read (11,*) a, flowType				! ���������� ��� ������������� �����
read (11,*) a, t_in					! ������ ����������� (�� �����)
read (11,*) a, p_in					! ������ �������� (�� �����)
read (11,*) a, alpha_in				! ���� ��������� ������ (�� �����)
read (11,*) a, p_out				! ����������� �������� (�� ������)
read (11,*) a, cfl					! ����� �������
read (11,*) a, charmin				! ����������� ����������� � ������������������ ����������
read (11,*) a, implicitInOutlet		! ������� ��������� ������� �� �����-�����
read (11,*) a, implicitPeriod		! ������� ��������� ������� �� �������� �������������
read (11,*) a, implicitWalls		! ������� ��������� ������� �� �������
read (11,*) a, coeff                ! ����������� ����� ������������ ���������
read (11,*) a, s_prec
read (11,*) a, epsilon				
read (11,*) a, NtimeMax				! ����� ����� �� �������
read (11,*) a, Ncontrol
read (11,*) a, Nderivatives
read (11,*) a, ResPeriod
read (11,*) a, SubitMax
read (11,*) a, epsilonSub
read (11,*) a, timeStepChoice		! ����� ���� ���� �� �������
read (11,*) a, GridChoice			! ����� �����
read (11,*) a, RRTypeChoice			! ����� ������� ���������� ������� (� �����, �� ��������)
read (11,*) a, ExplicitSchemeType	! ��� ����� �����
read (11,*) a, ImplicitSchemeType	! ��� ������� �����
read (11,*) a, timeDerivType		! 2-, 3-������� ����������� �� �������

close (11)



!if (flowType==1) then
!	write(*,*) 'Enter Mach number:'
!	read (*,*) Mach
!endif
!epsilon=1.e-7

!NtimeMax=3

! ��������� � ����������� ��������

gamm=gam-1.
gamp=gam+1.
gamo=1./gam
gammo=1./gamm
gampo=1./gamp
gampg=gamp/gam
gammg=gamm/gam
gammgp=gamm/gamp
gamgm=gam/gamm

if (alpha_in==45.) then
 alpha_in=atan(1.)
else 
 alpha_in=alpha_in*pi/180.
endif


! dimensionless

r_in=p_in/rg/t_in							! ��������� �� ����� �� ����������������

c_cr=dsqrt(2.*gam/gamp*p_in/r_in)           ! ����������� ��������

r_cr=(2.*gampo)**gammo*r_in					! ����������� ���������

r_in=r_in/r_cr								! ��������� �� ����� ���������������

!r_in=1.

p_in=p_in/r_cr/c_cr**2						! �������� �� ����� ���������������

!p_in=2.

t_in=p_in/r_in/rg							! ����������� �� ����� ���������������
											
i_in=gam/gamm*rg*t_in						! ��������� �� ����� ���������������

s_in=p_in/r_in**gam							! ����������� ������� �� ����� ���������������

p_out=p_out/r_cr/c_cr**2					! �������� �� ������ ���������������

!p_out=1.


end subroutine readFiles