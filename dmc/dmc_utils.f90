module dmc_utils
	implicit none
	
	integer, private, parameter :: DP = kind(0.0D0)
	
	type :: walker(d)
		integer, len :: d
		real(kind=DP) :: pos(d)
		type(walker(d)), pointer :: next => null()
	end type walker

contains

	! append walker into list
	subroutine append(x, last, d)
		integer, intent(in) :: d
		type(walker(d)), pointer, intent(inout) :: last
		real(kind=DP), intent(in) :: x(d)
		
		allocate(last%next)
		last => last%next
		last%pos = x

	end subroutine append

	! remove walker from list
	subroutine remove(previous, current, d)
		integer, intent(in) :: d
		type(walker(d)), pointer, intent(inout) :: previous, current

		previous%next => current%next
		deallocate(current)
		current => previous%next

	end subroutine remove

end module
