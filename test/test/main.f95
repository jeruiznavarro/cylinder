program main
    implicit none
    integer(8)::i=0
    integer(8)::j=0
    integer(8)::k=0
    integer(8)::ou=0
    integer(8)::array(2,2,2)=0.d0
    integer(8)::matrix(2,2000000)=0.d0
    integer(8)::vector(2000000)=0.d0
    real(8)::start=0.d0
    real(8)::endd=0.d0
    call cpu_time(start)
    forall(i=1:2,j=1,2000000)
        where(matrix(i,j)==3.d0) matrix(i,j)=1.d0
    end forall
    call cpu_time(endd)
    write(*,*) 'Time without loop:',endd-start
    call cpu_time(start)
    do i=1,2000000
        matrix(1,i)=1.d0
        matrix(2,i)=1.d0
    enddo
    call cpu_time(endd)
    write(*,*) 'Time with loop:',endd-start
end program main
