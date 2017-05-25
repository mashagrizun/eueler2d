module commonVariables

! ���������� ����� �����
integer :: Nx,Ny,Nz,Nx1,Ny1,Nz1,NN,NNx,NNy,N1,N2,TW,rmax,label,charmin,limiter
integer :: implicitWalls,implicitInOutlet,implicitPeriod,mark
! ����� ���� ���� �� �������, ���������� ����� �� ������� (�������� �� �������� ����� initial.dat)
integer :: timeStepChoice, NtimeMax,Ntime,total=0,Ncontrol,SubitMax,Nderivatives,flowType
! ������� ����������, ���������� ��������,(�������� �� initial.dat), ��������� � ���
real*8 :: rg,gam,gamm,gamp,gamo,gammo,gampo,gampg,gammg,gammgp,gamgm,del,eps,epsilon,epsilonSub
! ��������� ������� (�������� �� initial.dat)
real*8 :: t_in,p_in,alpha_in,s_in,i_in,r_in,p_out
! ����������� ��������, ����������� ���������
real*8 :: c_cr,r_cr
! ����� ������� (�������� �� initial.dat), ����������� ��� �� ������� �� ���� �������, ����� pi
real*8 :: cfl,tau_min,pi,romax,s_prec,maxdro,coeff,Mach
! ��������� ���������� ��� ���������� ������
character*40 :: format1, format2, format3

character*40 :: RRTypeChoice,ExplicitSchemeType,ImplicitSchemeType,GridChoice,timeDerivType

real*8 :: maxdro1=0.,res=0., RoMaxMax, RoMaximum, RoMaximum10
integer :: ResPeriod


end module commonVariables