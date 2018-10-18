  !Here I want to compare the speed difference in calling a subroutine if it qualifies for an if statement vs calling the subroutine with an if-statement withing

module test
  implicit none
contains
  subroutine testing( array, this )
    implicit none
    real, dimension(:,:), intent(in) :: array
    real, intent(in) :: this

    if (.true.) return
  end subroutine testing

end module test

program run
  use test
  implicit none
  real, dimension(1000,1000) :: array
  real :: T
  integer :: i

  real :: start, finish
  array = 0.0
  T = 0.0
  CALL CPU_TIME(start)
  do i = 1,1000000000
     if (.true. ) CALL testing(array, T)
  end do
  CALL CPU_TIME(finish)

  print*, "time (s): ", finish - start

end program run
