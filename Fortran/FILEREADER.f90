program FILEREADER

   real, dimension(:,:), allocatable :: Mrb
   integer :: n,m
   character (len=3) :: varname
   
   n = 2
   m = 2

   open (unit=99, file='array.txt', status='old', action='read')
!   read(99, *), n
!   read(99, *), m
   read(99, *), varname
   allocate(x(n,m))

   do I=1,n,1
      read(99,*) x(I,:)
      write(*,*) x(I,:)
   enddo

end