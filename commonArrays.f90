module commonArrays

! координаты сетки x,y,z, прочитанные из файла
real*8,allocatable :: xold(:,:,:),yold(:,:,:),zold(:,:,:)
! координаты центра ячейки, координаты узлов, координаты центров "горизонтальных" граней
real*8,allocatable :: xc(:,:),yc(:,:),xn(:,:),yn(:,:),xh(:,:),yh(:,:)
! координаты центров "вертикальных" граней
real*8,allocatable :: xv(:,:),yv(:,:)
! приращения по x и y  на "горизонтальных" и "вертикальных" гранях
real*8,allocatable :: dxh(:,:),dyh(:,:),dxv(:,:),dyv(:,:)
! длины граней
real*8,allocatable :: dksi(:,:),deta(:,:)
! метрические коэффициенты
real*8,allocatable :: xksi(:,:),yksi(:,:),xeta(:,:),yeta(:,:)
real*8,allocatable :: xksic(:,:),yksic(:,:),xetac(:,:),yetac(:,:)
! расстояния от центра до середины грани ячейки и между серединами противолежащих сторон
real*8,allocatable :: dh_half_lt(:,:),dv_half_dn(:,:),dh_half_rt(:,:),dv_half_up(:,:),dh(:,:),dvert(:,:)
! приращения по средним линиям
real*8,allocatable :: dxh_mid(:,:),dxv_mid(:,:),dyh_mid(:,:),dyv_mid(:,:)

! расстояния между центрами соседних ячеек, площадь ячейки
real*8,allocatable :: dnoh(:,:),dnov(:,:),S(:,:),Jac(:,:)

! вектор примитивных переменных, вектор фиктивных ячеек, вектор приращения по времени (третий индекс - номер компоненты)
real*8,allocatable :: q(:,:,:),qn(:,:,:),dq(:,:,:),dqn(:,:,:)
real*8,allocatable :: dfi(:,:,:,:),delq(:,:,:)
! векторы параметров на границах (размерность - на 1 больше, чем примитивных), нормальныt компоненты скорости
real*8,allocatable :: qx(:,:,:),qy(:,:,:),unx(:,:),uny(:,:),lamksi(:,:,:),lameta(:,:,:)
! временной шаг, приращения для выч. шага по времени
real*8,allocatable :: tau(:,:),h_x(:,:)
! компоненты векторов потоков
real*8,allocatable :: flx(:,:,:),fly(:,:,:)
! приращения потоков по Годунову
real*8,allocatable :: delta_god(:,:,:),s_dim(:,:)

real*4,allocatable :: uForF(:,:)

real*8,allocatable :: dqx(:,:,:),dqy(:,:,:),dqx_minus(:,:,:),dqy_minus(:,:,:),dqx_plus(:,:,:),dqy_plus(:,:,:),dqt(:,:,:)

real*8,allocatable :: ro1(:,:),ro2(:,:)

real*8,allocatable :: teta(:,:),delta_godn(:,:,:),par(:,:,:),ddq0(:,:,:),ddq1(:,:,:)
real*8,allocatable :: rq(:,:,:),rrq(:,:,:),dqtx(:,:,:),dqty(:,:,:)
real*4,allocatable :: r_f(:,:,:),u_f(:,:,:),v_f(:,:,:),w_f(:,:,:),p_f(:,:,:)

end module commonArrays