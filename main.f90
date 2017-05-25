program main2D
use commonVariables
use commonArrays
USE DFPORT


implicit none
integer :: i,j,k,NtimeStop=0,timearr1(3),timearr2(3),hour,min,sec,m
character*8 ch_time1,ch_time2
real*8 :: maxdelro

Ntime=0

call readFiles  ! чтение параметров

call arrayInitialisation  ! инициализация массивов

call geom       ! расчет геометрии сетки, метрических коэффициентов

call initial    ! задание начальных условий

call PrintFromF

!call checkPrint('q0    ')  ! печать начального поля

call ITIME(timearr1)  ! время начала расчета

write(*,*) timearr1(1),':',timearr1(2),':',timearr1(3)

call time(ch_time1)  ! время начала расчета - вывод на экран

write(*,*)
write(*,*) 'Start: ', ch_time1 

! пока максимальная невязка на слое больше эпсилон, или шагов по времени сделано меньше заданного
! количества, проводить расчет

do while (romax>epsilon.or.Ntime<NtimeMax)

 qn=q
 dq=0.
 delq=0.

! выбор типа шага по времени

 select case (timeStepChoice) 
 case(1)                        ! локальный в каждой ячейке
  call timeStep 
 case(2)                        ! локальный на слое 
  call timeStep
  tau=tau_min
!  tau=2.00045e-3
! write(*,*) tau_min
 case(3)                        ! глобальный (cfl задается малым)
  if (Ntime==0) then
   call timeStep
   tau=tau_min
  endif
 end select

 Ntime=Ntime+1
 call calculation2D         ! расчет

 write(*,'(a19,i9,a1,i9)') 'Time layer number=',Ntime,'/',NtimeMax
 
enddo

write(*,*) 'total', total

call time(ch_time2)
call ITIME(timearr2)

if (timearr1(3)>timearr2(3)) then
 sec=timearr2(3)+60-timearr1(3)
 timearr2(2)=timearr2(2)-1
else
 sec=timearr2(3)-timearr1(3)
endif

if (timearr1(2)>timearr2(2)) then
 min=timearr2(2)+60-timearr1(2)
 timearr2(1)=timearr2(1)-1
else
 min=timearr2(2)-timearr1(2)
endif

hour=timearr2(1)-timearr1(1)

call PrintFromF

!call checkPrint('q     ')

write(*,*)
write(*,*) 'Start: ', ch_time1
write(14,*) 'Start: ', ch_time1
write(*,*) 'Stop: ', ch_time2
write(14,*) 'Stop: ', ch_time2
write(14,*) 'Ntime=',Ntime
write(14,*) 'Total=',total
write(14,*) hour,':',min,':',sec
write(*,*) hour,':',min,':',sec

call arrayDeinitialisation

close(13); close(14); close(15); close(16)
close(17); close(18); close(19); close(20)
close(50); close(60); close(100); close(101)
close(102);close(103);close(104);close(105)
close(106);close(107);close(221)

end program main2D
