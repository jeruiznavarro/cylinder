program main
    implicit none
    integer(8)::i=0
    integer(8)::j=0
    integer(8)::k=0
    integer(8)::ou=0
    integer(8)::array(2,2,2)=0.d0
    integer(8)::matrix(2,2)=0.d0
    integer(8)::vector(2)=0.d0
    do i=1,2
        vector(i)=i
        do j=1,2
            matrix(i,j)=i+j
            do k=1,2
                array(i,j,k)=i+j+k
            enddo
        enddo
    enddo
    (((k,j,array(i,j,k),i=1,2),j=1,2),k=1,2)
    ou=array(i,j,k)
    write(*,*)
end
