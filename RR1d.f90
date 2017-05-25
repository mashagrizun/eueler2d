! распад разрыва

subroutine RR1D(r11,u11,v11,p11,r22,u22,v22,p22)
use commonVariables

real*8, intent(inout) :: r11,u11,v11,p11,r22,u22,v22,p22
real*8 :: du,c11,c1c,c22,c2c,pkr,uu,uud
real*8 :: ura,uva,pa,a11,d11,d1c,a22,d22,d2c,dp,zap

integer :: ii,l

! О Р И Е Н Т А Ц И Я

l=0
if (p11.gt.p22) then
  l=1
  zap=p11
  p11=p22
  p22=zap
  zap=r11
  r11=r22
  r22=zap
  zap=u11
  u11=-u22
  u22=-zap
  zap=v11
  v11=v22
  v22=zap
endif

du=u11-u22
c11=sqrt(gam*p11/r11)
c22=sqrt(gam*p22/r22)


! Т О Ч Н Ы Й   Р А С П А Д

ii=0

pkr=(p11*r22*c22+p22*r11*c11+du*r11*c11*r22*c22)/(r11*c11+r22*c22)

uud=(p22-p11)/sqrt(.5*r11*(gamp*p22+gamm*p11))
ura=2.*c22*gammo*((p11/p22)**(.5*gammg)-1.)
uva=-2.*gammo*(c11+c22)

! КОНФИГУРАЦИЯ РАСПАДА

if (du.gt.uud) then
! Две ударных волны
 10	    continue
 ii=ii+1
 p=pkr
 pkr=p-((pkr-p11)/(r11*c11*sqrt(.5*(gampg*pkr/p11+gammg)))+(pkr-p22)/(r22*c22*sqrt(.5*(gampg*pkr/p22+gammg)))-du)/(.25*(gamp*pkr/p11+3.*gam-1.)/(gam*r11*c11*sqrt(.125*(gampg*pkr/p11+gammg)**3))+.25*(gamp*pkr/p22+3.*gam-1.)/(gam*r22*c22*sqrt(.125*(gampg*pkr/p22+gammg)**3)))
 if (abs(pkr-p)/p.gt.0.000001.and.ii.lt.10) goto 10
  a11=sqrt(.5*r11*(gamp*pkr+gamm*p11))
  a22=sqrt(.5*r22*(gamp*pkr+gamm*p22))
  d11=u11-a11/r11
  d22=u22+a22/r22
  if (d11*d22.gt.0.) then
   if (d11.lt.0.) goto 180
    goto 190
   else
	u11=(a11*u11+a22*u22+p11-p22)/(a11+a22)
	if (u11.gt.0.) goto 110
	 goto 120
	endif
   elseif (du.gt.ura) then

! Левая ударная волна и правая волна разрежения

    20	    continue
	ii=ii+1
	p=pkr
	pkr=p-((pkr-p11)/(r11*c11*sqrt(.5*(gampg*pkr/p11+gammg)))+2.*gammo*c22*((pkr/p22)**(.5*gammg)-1.)-du)/(.25*(gamp*pkr/p11+3.*gam-1.)/(gam*r11*c11*sqrt(.125*(gampg*pkr/p11+gammg)**3))+gamo/pkr*c22*(pkr/p22)**(.5*gammg))
	if (abs(pkr-p)/p.gt.0.000001.and.ii.lt.10) goto 20
     a11=sqrt(.5*r11*(gamp*pkr+gamm*p11))
     if (abs(1.-pkr/p22).le.1.0e-5) then
      a22=r22*c22
     else
      a22=.5*gammg*r22*c22*(1.-pkr/p22)/(1.-(pkr/p22)**(.5*gammg))
     endif
     uu=(a11*u11+a22*u22+p11-p22)/(a11+a22)
     d11=u11-a11/r11
     d22=u22+c22
     c2c=c22-.5*gamm*(u22-uu)
     d2c=uu+c2c
     if (d11*d22.gt.0.) then
      if (d11.lt.0.) goto 180
       goto 190
      else
       if (uu.gt.0.) then
        u11=uu
        goto 110
       elseif (d2c.gt.0.) then
        goto 130
       else
        goto 140
       endif
      endif
     elseif (du.gt.uva) then

! Две волны разрежения

      pkr=p11*((du-uva)/(ura-uva))**(2.*gamgm)
      if (abs(1.-pkr/p11).le.1.0e-05) then
       a11=r11*c11
      else
       a11=.5*gammg*r11*c11*(1.-pkr/p11)/(1.-(pkr/p11)**(.5*gammg))
      endif
      if (abs(1.-pkr/p22).le.1.0e-05) then
       a22=r22*c22
      else
       a22=.5*gammg*r22*c22*(1.-pkr/p22)/(1.-(pkr/p22)**(.5*gammg))
      endif
      uu=(a11*u11+a22*u22+p11-p22)/(a11+a22)
      d11=u11-c11
      c1c=c11+.5*gamm*(u11-uu)
      d1c=uu-c1c
      d22=u22+c22
      c2c=c22-.5*gamm*(u22-uu)
      d2c=uu+c2c
      if (d11*d22.gt.0.) then
       if (d11.lt.0.) goto 180
        goto 190
       else
        if (d1c.gt.0.) goto 150
        if (uu.gt.0.) goto 160
        if (d2c.gt.0.) goto 130
         goto 140
	   endif
	  else
! Вакуум
       goto 170
	  endif

! ПАРАМЕТРЫ НА ГРАНИ

!	Между левой ударной волной и контактным разрывом
      110	continue
      r11=r11*(gamp*pkr+gamm*p11)/(gamm*pkr+gamp*p11)
	  p11=pkr
      goto 200
!	Между контактным разрывом и правой ударной волной
      120	continue

	  r11=r22*(gamp*pkr+gamm*p22)/(gamm*pkr+gamp*p22)
	  p11=pkr
	  v11=v22
	  goto 200
!	Между контактным разрывом и правой волной разрежения
      130	continue

	  p11=pkr
	  r11=gam*pkr/c2c**2
	  u11=uu
	  v11=v22
	  goto 200
!	Внутри правой волны разрежения
      140	continue

	  c2c=2.*gampo*c22-gammgp*u22
	  p11=p22*(c2c/c22)**(2.*gamgm)
	  r11=gam*p11/c2c**2
	  u11=-c2c
	  v11=v22
	  goto 200
!	Внутри левой волны разрежения
      150	continue

	  c1c=2.*gampo*c11+gammgp*u11
	  p11=p11*(c1c/c11)**(2.*gamgm)
	  r11=gam*p11/c1c**2
	  u11=c1c
	  goto 200
!	Между левой волной разрежения и контактным разрывом
      160	continue

	  p11=pkr
	  r11=gam*pkr/c1c**2
	  u11=uu
	  goto 200

!	Вакуум
      170	continue
	  r11=(r11+r22)*(0.001/(p11+p22))**gamo
	  p11=0.001
	  u11=0.
	  v11=0.
	  goto 200

!	Волны ушли влево
      180	continue

	  p11=p22
	  r11=r22
	  u11=u22
	  v11=v22
      goto 200

!	Волны ушли вправо	
      190	continue

!	(сохраняются параметры P11,R11,U11,V11)

      200	continue
	  if (l.eq.1) u11=-u11
	  return
end
