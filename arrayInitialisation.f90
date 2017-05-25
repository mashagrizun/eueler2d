subroutine arrayInitialisation
use commonVariables
use commonArrays
implicit none


allocate (xold(Nx1,Ny1,Nz1),yold(Nx1,Ny1,Nz1),zold(Nx1,Ny1,Nz1),&
&         xc(0:Nz+1,0:Ny+1),yc(0:Nz+1,0:Ny+1),xn(0:Nz+2,0:Ny+2),yn(0:Nz+2,0:Ny+2),&
&         xh(Nz,Ny1),yh(Nz,Ny1),xv(Nz1,Ny),yv(Nz1,Ny),&
&         dxh(Nz,Ny1),dyh(Nz,Ny1),dxv(Nz1,Ny),dyv(Nz1,Ny),&
&         dksi(Nz,Ny1),deta(Nz1,Ny),&
&         xksi(Nz,Ny1),yksi(Nz,Ny1),&
&         xeta(Nz1,Ny),yeta(Nz1,Ny),&
&         xksic(Nz,0:Ny1),yksic(Nz,0:Ny1),&
&         xetac(0:Nz1,Ny),yetac(0:Nz1,Ny),&
&         dh_half_rt(Nz,Ny),dv_half_dn(Nz,Ny),dh_half_lt(Nz,Ny),dv_half_up(Nz,Ny),dh(Nz,Ny),dvert(Nz,Ny),&
&         dxh_mid(Nz,Ny),dxv_mid(Nz,Ny),dyh_mid(Nz,Ny),dyv_mid(Nz,Ny),&

&         dnoh(1:Nz+1,Ny),dnov(Nz,1:Ny+1),S(Nz,Ny),Jac(Nz,Ny),&

&         q(Nz,Ny,4),qn(Nz,Ny,4),dq(Nz,Ny,4),dqn(Nz,Ny,4),&
&		  dfi(0:Nz+1,0:Ny+1,4,4),delq(Nz,Ny,4),&

&		  qx(Nz+1,Ny,4),qy(Nz,Ny+1,4),unx(Nz+1,Ny),uny(Nz,Ny+1),tau(Nz,Ny),&
&         flx(Nz+1,Ny,4),fly(Nz,Ny+1,4),lamksi(0:Nz+1,Ny,4),lameta(Nz,0:Ny+1,4),&

&		  delta_god(Nz,Ny,4),h_x(Nz,Ny),s_dim(Nz+1,Ny),uForF(Nz,Ny),&
&         dqx(Nz,Ny,4),dqy(Nz,Ny,4),dqx_minus(Nz,Ny,4),dqy_minus(Nz,Ny,4),dqx_plus(Nz,Ny,4),dqy_plus(Nz,Ny,4),dqt(Nz,Ny,4),&

&         ro1(Nz,Ny),ro2(Nz,Ny),&

&         teta(Nz,Ny),delta_godn(Nz,Ny,4),par(Nz+1,Ny+1,4),ddq0(Nz,Ny,4),ddq1(Nz,Ny,4),&
&		  rq(Nz+1,Ny+1,4),rrq(0:Nz+1,0:Ny+1,4),dqtx(Nz,Ny,4),dqty(Nz,Ny,4),&
&		  r_f(Nx,Ny,Nz),u_f(Nx,Ny,Nz),v_f(Nx,Ny,Nz),w_f(Nx,Ny,Nz),p_f(Nx,Ny,Nz))

delta_god=0.
dqn=0.
delq=0.
dfi=0.
dq=0.

delta_godn=0.
ddq0=0.
ddq1=0.

end subroutine arrayInitialisation

