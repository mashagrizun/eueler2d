subroutine arrayDeInitialisation
use commonArrays
implicit none

deallocate (xold,yold,zold,xc,yc,xn,yn,&
&         xh,yh,xv,yv,dxh,dyh,dxv,dyv,&
&         dksi,deta,xksi,yksi,xeta,yeta,&
&         xksic,yksic,xetac,yetac,&
&         dh_half_rt,dv_half_dn,dh_half_lt,dv_half_up,dh,dvert,&
&         dxh_mid,dxv_mid,dyh_mid,dyv_mid,&

&         dnoh,dnov,S,Jac,q,qn,dq,dqn,dfi,delq,qx,qy,unx,uny,lamksi,lameta,tau,flx,fly,&
&		  delta_god,h_x,s_dim,uForF,dqx,dqy,dqx_minus,dqy_minus,dqx_plus,dqy_plus,dqt,&
&         ro1,ro2,&
&         teta,delta_godn,par,ddq0,ddq1,&
&		  rq,rrq,dqtx,dqty,&
&		  r_f,u_f,v_f,w_f,p_f)
	
end subroutine arrayDeInitialisation