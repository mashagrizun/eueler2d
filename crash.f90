subroutine crash (mark1,qq1,qq4,ii,jj)

use commonVariables

integer, intent(in) :: mark1
real*8, intent(in) :: qq1,qq4
integer,intent(in) :: ii,jj


		if (qq1<0.) then
			write(221,*) 'Ntime=', Ntime
			write(221,*) 'mark1=',mark1
			
			select case(mark1)

			case (0)
				write(221,*) 'before TimeStep'
			case (1)
				write(221,*) 'before Flux X'
			case (2)
				write(221,*) 'before Flux Y'
			case (3)
				write(221,*) 'before BoundsPeriod'
			case (4)
				write(221,*) 'before BoundsWallDown'
			case (5)
				write(221,*) 'before BoundsWallUp'

			end select

			write(221,*) ii,jj, ' ro=', qq1
		endif

		if (qq4<0.) then
			write(221,*) 'Ntime=', Ntime
			write(221,*) 'mark1=',mark1

			select case(mark1)

			case (0)
				write(221,*) 'before TimeStep'
			case (1)
				write(221,*) 'before Flux X'
			case (2)
				write(221,*) 'before Flux Y'
			case (3)
				write(221,*) 'before BoundsPeriod'
			case (4)
				write(221,*) 'before BoundsWallDown'
			case (5)
				write(221,*) 'before BoundsWallUp'

			end select
			write(221,*) ii,jj, ' p=', qq4
		endif

end subroutine crash