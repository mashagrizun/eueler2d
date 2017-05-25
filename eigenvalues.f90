subroutine eigenvalues
use commonArrays
use commonVariables

implicit none

real*8 :: uksi,ueta,c,uu,vv,un,ut
integer :: i,j

 do i=0,NNx+1
   do j=1,NNy
    if (i==0) then

		c=dsqrt(gam*qx(1,j,4)/qx(1,j,1))
		uksi=qx(1,j,2)*yetac(0,j)-qx(1,j,3)*xetac(0,j)

		lamksi(0,j,1)=uksi-c
		lamksi(0,j,2)=uksi
		lamksi(0,j,3)=uksi
		lamksi(0,j,4)=uksi+c

    elseif (i==NNx+1) then

		c=dsqrt(gam*qx(NNx+1,j,4)/qx(NNx+1,j,1))
		uksi=qx(NNx+1,j,2)*yetac(NNx+1,j)-qx(NNx+1,j,3)*xetac(NNx+1,j)

		lamksi(NNx+1,j,1)=uksi-c
		lamksi(NNx+1,j,2)=uksi
		lamksi(NNx+1,j,3)=uksi
		lamksi(NNx+1,j,4)=uksi+c

	else

		c=dsqrt(gam*q(i,j,4)/q(i,j,1))
		uksi=q(i,j,2)*yetac(i,j)-q(i,j,3)*xetac(i,j)

		lamksi(i,j,1)=uksi-c
		lamksi(i,j,2)=uksi
		lamksi(i,j,3)=uksi
		lamksi(i,j,4)=uksi+c

	endif
   enddo
 enddo

 do i=1,NNx 

	 do j=1,NNy

		c=dsqrt(gam*q(i,j,4)/q(i,j,1))
		ueta=q(i,j,3)*xksic(i,j)-q(i,j,2)*yksic(i,j)

		lameta(i,j,1)=ueta-c
		lameta(i,j,2)=ueta
		lameta(i,j,3)=ueta
		lameta(i,j,4)=ueta+c
	 
	 enddo

     ueta=qy(i,1,3)*xksic(i,0)-qy(i,1,2)*yksic(i,0)
     c=dsqrt(gam*qy(i,1,4)/qy(i,1,1))

	 lameta(i,0,1)=ueta-c
	 lameta(i,0,2)=ueta
	 lameta(i,0,3)=ueta
	 lameta(i,0,4)=ueta+c

     ueta=qy(i,NNy+1,3)*xksic(i,NNy+1)-qy(i,NNy+1,2)*yksic(i,NNy+1)
     c=dsqrt(gam*qy(i,NNy+1,4)/qy(i,NNy+1,1))

	 lameta(i,NNy+1,1)=ueta-c
     lameta(i,NNy+1,2)=ueta
	 lameta(i,NNy+1,3)=ueta
	 lameta(i,NNy+1,4)=ueta+c
   
 enddo

end subroutine eigenvalues