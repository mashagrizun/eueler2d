module commonArrays

! ���������� ����� x,y,z, ����������� �� �����
real*8,allocatable :: xold(:,:,:),yold(:,:,:),zold(:,:,:)
! ���������� ������ ������, ���������� �����, ���������� ������� "��������������" ������
real*8,allocatable :: xc(:,:),yc(:,:),xn(:,:),yn(:,:),xh(:,:),yh(:,:)
! ���������� ������� "������������" ������
real*8,allocatable :: xv(:,:),yv(:,:)
! ���������� �� x � y  �� "��������������" � "������������" ������
real*8,allocatable :: dxh(:,:),dyh(:,:),dxv(:,:),dyv(:,:)
! ����� ������
real*8,allocatable :: dksi(:,:),deta(:,:)
! ����������� ������������
real*8,allocatable :: xksi(:,:),yksi(:,:),xeta(:,:),yeta(:,:)
real*8,allocatable :: xksic(:,:),yksic(:,:),xetac(:,:),yetac(:,:)
! ���������� �� ������ �� �������� ����� ������ � ����� ���������� �������������� ������
real*8,allocatable :: dh_half_lt(:,:),dv_half_dn(:,:),dh_half_rt(:,:),dv_half_up(:,:),dh(:,:),dvert(:,:)
! ���������� �� ������� ������
real*8,allocatable :: dxh_mid(:,:),dxv_mid(:,:),dyh_mid(:,:),dyv_mid(:,:)

! ���������� ����� �������� �������� �����, ������� ������
real*8,allocatable :: dnoh(:,:),dnov(:,:),S(:,:),Jac(:,:)

! ������ ����������� ����������, ������ ��������� �����, ������ ���������� �� ������� (������ ������ - ����� ����������)
real*8,allocatable :: q(:,:,:),qn(:,:,:),dq(:,:,:),dqn(:,:,:)
real*8,allocatable :: dfi(:,:,:,:),delq(:,:,:)
! ������� ���������� �� �������� (����������� - �� 1 ������, ��� �����������), ���������t ���������� ��������
real*8,allocatable :: qx(:,:,:),qy(:,:,:),unx(:,:),uny(:,:),lamksi(:,:,:),lameta(:,:,:)
! ��������� ���, ���������� ��� ���. ���� �� �������
real*8,allocatable :: tau(:,:),h_x(:,:)
! ���������� �������� �������
real*8,allocatable :: flx(:,:,:),fly(:,:,:)
! ���������� ������� �� ��������
real*8,allocatable :: delta_god(:,:,:),s_dim(:,:)

real*4,allocatable :: uForF(:,:)

real*8,allocatable :: dqx(:,:,:),dqy(:,:,:),dqx_minus(:,:,:),dqy_minus(:,:,:),dqx_plus(:,:,:),dqy_plus(:,:,:),dqt(:,:,:)

real*8,allocatable :: ro1(:,:),ro2(:,:)

real*8,allocatable :: teta(:,:),delta_godn(:,:,:),par(:,:,:),ddq0(:,:,:),ddq1(:,:,:)
real*8,allocatable :: rq(:,:,:),rrq(:,:,:),dqtx(:,:,:),dqty(:,:,:)
real*4,allocatable :: r_f(:,:,:),u_f(:,:,:),v_f(:,:,:),w_f(:,:,:),p_f(:,:,:)

end module commonArrays