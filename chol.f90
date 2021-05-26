program chol
    implicit none
    double precision, allocatable :: A(:,:), F(:,:), L(:,:), G(:,:); 
    integer :: i,j,k,n;
    real*8 :: start, finish;
    double precision :: temp, error;
    do n = 500,5000,500
    allocate(A(n,n));
    allocate(L(n,n));
    allocate(G(n,n));
    allocate(F(n,n));
    !< scalar cholesky decomposite 
    call generate_data;
    F = A;
    call cpu_time(start);
    call chol_decomp_col(A);
    call cpu_time(finish);
    do j=1,n
        L(j:n,j)=A(j:n,j);
    end do
    ! show matrix L, this matrix A = L*tL
    error = maxval(abs(abs(G)-abs(matmul(L, transpose(L)))));
    print *,n,error,finish-start;
    call cpu_time(start);
    call chol_decomp_row(F);
    call cpu_time(finish);
    do j=1,n
        L(j:n,j)=F(j:n,j);
    end do
    ! show matrix L, this matrix A = L*tL
    error = maxval(abs(abs(G)-abs(matmul(L, transpose(L)))));
    print *,n,error,finish-start;
    
    ! < row version 
    ! call cpu_time(start);
    ! call generate_data;
    ! call chol_decomp_row;
    ! call cpu_time(finish);
    
    deallocate(A);
    deallocate(L);
    deallocate(G);
    deallocate(F);
    end do
    contains
    subroutine chol_decomp_col(A)
    integer :: i,j,k;
    double precision, allocatable :: A(:,:)
    do  i = 1,n
        A(i,i) = SQRT (A(i, i));
        do  j = i+1, n
            A(j,i) = A(j,i)/A(i,i);
        end do
        do  k=i+1,n
            do j = k, n
                A(j,k) = A(j,k) - A(j,i)*A(k,i);
            end do
        end do
    end do
    end subroutine chol_decomp_col
    subroutine chol_decomp_col_with_temp()
        do k=1,n
            temp=sqrt(A(k,k));
            do i=k,n
                A(i,k)=A(i,k)/temp;
            enddo
            do j=k+1,n
                do i=j,n
                    A(i,j)=A(i,j)-A(i,k)*A(j,k);
                enddo
            enddo
        enddo
    end subroutine chol_decomp_col_with_temp
    subroutine generate_data()
        call random_number(A)
        A = 0.5*(A + transpose(A));
            do i=1,n
                A(i,i) = float(N) + A(i,i);
            enddo
            G=0;
            G=A;
    end subroutine Generate_Data
    subroutine chol_decomp_row(A)
    double precision, allocatable :: A(:,:)
    do k=1,n
        A(k:n,k)=A(k:n,k)/sqrt(A(k,k));
        do i=k+1,n
            A(i,k+1:i)=A(i,k+1:i)-A(i,k)*A(k+1:i,k);
        end do
    end do
    end subroutine chol_decomp_row
    subroutine chol_decomp_mpi

    end subroutine chol_decomp_mpi
end program chol
