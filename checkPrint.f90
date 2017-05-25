! процедура записи параметров в файлы(для проверки)

subroutine checkPrint(printCheck)
use commonArrays
use commonVariables

implicit none

character*6,intent(in) :: printCheck
character*10 :: name
integer :: i,j,k,fileNumb,numbx,numby

select case (printCheck)

case ('q0    ')                 ! примитивные переменные на начальном слое
 open (13, file = 'check_parameters1.dat',status='old')
 fileNumb=13
 name='q'
 numbx=NNx
 numby=NNy
 do i=1,NNx
  do j=1,NNy
   par(i,j,:)=q(i,j,:)
  enddo
 enddo
case ('q     ')                 ! примитивные переменные после расчета 
 open (14, file = 'check_parameters2.dat',status='old')
 fileNumb=14
 name='q'
 numbx=NNx
 numby=NNy
 do i=1,NNx
  do j=1,NNy
   par(i,j,:)=q(i,j,:)
  enddo
 enddo
case ('qx    ')                 ! значения на сторонах ячеек по направлению х
 open (15, file = 'check_qx_qy.dat',status='old')
 fileNumb=15
 name='qx'
 numbx=NNx+1
 numby=NNy
 do i=1,NNx+1
  do j=1,NNy
   par(i,j,:)=qx(i,j,:)
  enddo
 enddo
case ('qy    ')                 ! значения на сторонах ячеек по направлению у
 open (15, file = 'check_qx_qy.dat',status='old')
 fileNumb=15
 name='qy'
 numbx=NNx
 numby=NNy+1
 do i=1,NNx
  do j=1,NNy+1
   par(i,j,:)=qy(i,j,:)
  enddo
 enddo
case('flx   ')                  ! вектор потока по направлению х
 open (16, file = 'check_flx.dat',status='old')
 fileNumb=16
 name='flx'
 numbx=NNx+1
 numby=NNy
 do i=1,NNx+1
  do j=1,NNy
   par(i,j,:)=flx(i,j,:)
  enddo
 enddo
case ('fly   ')                 ! вектор потока по направлению у
 open (18, file = 'check_fly.dat',status='old')
 fileNumb=18
 name='fly'
 numbx=NNx
 numby=NNy+1
 do i=1,NNx
  do j=1,NNy+1
   par(i,j,:)=fly(i,j,:)
  enddo
 enddo
case('delta ')                  ! приращения по Годунову
 open (17, file = 'check_delta.dat',status='old')
 fileNumb=17
 name='delta'
 numbx=NNx
 numby=NNy
 do i=1,NNx
  do j=1,NNy
   par(i,j,:)=delta_god(i,j,:)
  enddo
 enddo
end select

write (fileNumb,'(a15,i5)') 'Time iteration=',Ntime
!do k=1,4
! write(fileNumb,format1) name,k
! write(fileNumb,format2) 'i=',(i,i=1,numbx)
! do j=1,numby
!  write(fileNumb,format3) 'j=',j,(par(i,j,k),i=1,numbx)
! enddo
! write(fileNumb,*)
!enddo

write(fileNumb,format2) 'i=',(i,i=1,numbx)
do k=1,4
 write(fileNumb,format1) name,k
  do j=1,numby
  write(fileNumb,format3) 'j=',j,(par(i,j,k),i=1,numbx)
 enddo
 write(fileNumb,*)
enddo

write(fileNumb,*) 'gam=', gam
write(fileNumb,*) 'rg=', rg
write(fileNumb,*) 't_in=', t_in
write(fileNumb,*) 'p_in=', p_in
write(fileNumb,*) 'alpha_in=', alpha_in
write(fileNumb,*) 'p_out=', p_out
write(fileNumb,*) 'cfl=', cfl
write(fileNumb,*) 'charmin=', charmin
write(fileNumb,*) 's_prec=', s_prec
write(fileNumb,*) 'epsilon=', epsilon
write(fileNumb,*) 'NtimeMax=', NtimeMax
write(fileNumb,*) 'Ncontrol=', Ncontrol
write(fileNumb,*) 'epsilonSub=', epsilonSub
write(fileNumb,*) 'timeStepChoice=', timeStepChoice
write(fileNumb,*) 'GridChoice=', GridChoice
write(fileNumb,*) 'RRTypeChoice=', RRTypeChoice
write(fileNumb,*) 'ExplicitSchemeType=', ExplicitSchemeType
write(fileNumb,*) 'ImplicitSchemeType=', ImplicitSchemeType

end subroutine checkPrint