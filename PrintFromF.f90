subroutine PrintFromF
use commonVariables
use commonArrays

implicit none

integer :: k,ii,m,j,i,mid,NNN

NNN=Nx*Ny*Nz
!NNN=1*Ny*Nz

!open (31, file = 'r_f.dat', status='old')
!open (21, file = 'r1f.dat', form='unformatted', access='direct', recl=4*NNN)
!read (21, rec=1)  (((r_f(ii,k,m), k=1,Ny), ii=1,Nx), m=1,Nz) 

!open (32, file = 'u_f.dat', status='old')
!open (22, file = 'u1f.dat', form='unformatted', access='direct', recl=4*NNN)
!read (22, rec=1)  (((u_f(ii,k,m), k=1,Ny), ii=1,Nx), m=1,Nz)

!open (33, file = 'v_f.dat', status='old')
!open (23, file = 'v1f.dat', form='unformatted', access='direct', recl=4*NNN)
!read (23, rec=1)  (((v_f(ii,k,m), k=1,Ny), ii=1,Nx), m=1,Nz)


!open (24, file = 'w1f.dat', form='unformatted', access='direct', recl=4*NNN)
!read (24, rec=1)  (((w_f(ii,k,m), k=1,Ny), ii=1,Nx), m=1,Nz)

!open (34, file = 'p_f.dat', status='old')
!open (25, file = 'p1f.dat', form='unformatted', access='direct', recl=4*NNN)
!read (25, rec=1)  (((p_f(ii,k,m), k=1,Ny), ii=1,Nx), m=1,Nz)

!mid=Nx1/2

!ii=mid

!write(31,*) 'r'
!write(31,format2) 'i=',(m,m=1,Nz)
! do k=1,Ny
!  write(31,format3) 'j=',k,(r_f(ii,k,m),m=1,Nz)
! enddo
!write(31,*)

!write(32,*) 'u'
!write(32,format2) 'i=',(m,m=1,Nz)
! do k=1,Ny
!  write(32,format3) 'j=',k,(u_f(ii,k,m),m=1,Nz)
! enddo
!write(32,*)

!write(33,*) 'v'
!write(33,format2) 'i=',(m,m=1,Nz)
! do k=1,Ny
!  write(33,format3) 'j=',k,(v_f(ii,k,m),m=1,Nz)
! enddo
!write(33,*)

!write(34,*) 'p'
!write(34,format2) 'i=',(m,m=1,Nz)
! do k=1,Ny
!  write(34,format3) 'j=',k,(p_f(ii,k,m),m=1,Nz)
! enddo
!write(34,*)
             
close(21)
close(22)
close(23)
close(24)
close(25)

!close(31)
!close(32)
!close(33)
!close(34)


! запись параметров в файлы формата решателя F для проверки

uForF=0.                  !! компонента скорости w из решателя F - здесь является компонентой u
Nx=1
!open (40,file='rer1.dat',form='unformatted',access='direct',recl=4*Nx*Ny*Nz)
open (40,file='r1.dat',form='unformatted',access='direct',recl=4*Nx*Ny*Nz)
write(40,rec=1) ((sngl(q(i,j,1)),j=1,Ny),i=1,Nz)

!open (41,file='reu1.dat',form='unformatted',access='direct',recl=4*Nx*Ny*Nz)
open (41,file='u1.dat',form='unformatted',access='direct',recl=4*Nx*Ny*Nz)
write(41,rec=1) ((uForF(i,j),j=1,Ny),i=1,Nz)

!open (42,file='rev1.dat',form='unformatted',access='direct',recl=4*Nx*Ny*Nz)
open (42,file='v1.dat',form='unformatted',access='direct',recl=4*Nx*Ny*Nz)
write(42,rec=1) ((sngl(q(i,j,3)),j=1,Ny),i=1,Nz)

!open (43,file='rew1.dat',form='unformatted',access='direct',recl=4*Nx*Ny*Nz)
open (43,file='w1.dat',form='unformatted',access='direct',recl=4*Nx*Ny*Nz)
write(43,rec=1) ((sngl(q(i,j,2)),j=1,Ny),i=1,Nz)

!open (44,file='rep1.dat',form='unformatted',access='direct',recl=4*Nx*Ny*Nz)
open (44,file='p1.dat',form='unformatted',access='direct',recl=4*Nx*Ny*Nz)
write(44,rec=1) ((sngl(q(i,j,4)),j=1,Ny),i=1,Nz)

close(40); close(41); close(42); close(43); close(44)

open (219, file='promNcontrol.dat')
write(219,*) Ntime,(((q(i,j,k),i=1,NNx),j=1,NNy),k=1,4)


open (319, file='promNcontrol2.dat')
write(319,*) (((qx(i,j,k),i=1,NNx+1),j=1,NNy),k=1,4)


open (419, file='promNcontrol3.dat')
write(419,*) (((qy(i,j,k),i=1,NNx),j=1,NNy+1),k=1,4)


open (519, file='promNcontrol4.dat')
write(519,*) (((dqn(i,j,k),i=1,NNx),j=1,NNy),k=1,4)

close (219);close (319);close (419);close (519)

!open (11010101,file='resultsr.dat')
!open (11010102,file='resultsu.dat')
!open (11010103,file='resultsv.dat')
!open (11010104,file='resultsw.dat')
!open (11010105,file='resultsp.dat')


! do i=1,NNx
!  do j=1,NNy

!	if (abs(r_f(1,j,i)-q(i,j,1))>1.e-10) then
!		write(11010101,*) 'ro',i,j, r_f(1,j,i), sngl(q(i,j,1)), abs(r_f(1,j,i)-sngl(q(i,j,1)))
!	endif

!	if (abs(u_f(1,j,i)-0.0)>1.e-10) then
!		write(11010102,*) 'u',i,j, u_f(1,j,i), 0.0, abs(u_f(1,j,i)-0.0)
!	endif

!	if (abs(v_f(1,j,i)-q(i,j,3))>1.e-10) then
!		write(11010103,*) 'v',i,j, v_f(1,j,i), sngl(q(i,j,3)), abs(v_f(1,j,i)-sngl(q(i,j,3)))
!	endif
!
!	if (abs(w_f(1,j,i)-q(i,j,2))>1.e-10) then
!		write(11010104,*) 'w',i,j, w_f(1,j,i), sngl(q(i,j,2)), abs(w_f(1,j,i)-sngl(q(i,j,2)))
!	endif

!	if (abs(p_f(1,j,i)-q(i,j,4))>1.e-10) then
!		write(11010105,*) 'p',i,j, p_f(1,j,i), sngl(q(i,j,4)), abs(p_f(1,j,i)-sngl(q(i,j,4)))
!	endif

!  enddo
! enddo   

!close(11010101)
!close(11010102)
!close(11010103)
!close(11010104)
!close(11010105)


end subroutine PrintFromF