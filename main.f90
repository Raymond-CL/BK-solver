program main
  use BK
  implicit none
  real :: tstart,tfinish
  call cpu_time(tstart)


  call setBK()

  !call BKinfo()

  !call runBK()

  !call printBK()


  call cpu_time(tfinish)
  write(*,*) "time elapsed:",tfinish-tstart,"seconds."
end program main
