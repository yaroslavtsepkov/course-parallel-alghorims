program chol_decomposite_basic
    implicit none
    real*8, allocatable :: A(:,:),C(:,:),G(:,:);
    real*8 mem;
    integer N;
    integer i,j,k;

    N = 3;
    !< allocate memory
    allocate(A(1:N,1:N)); 
    allocate(C(1:N,1:N));   
    allocate(G(1:N,1:N));

    call random_number(A);
    G=0;
    A=0.5*(A+transpose(A))
    do i=1,N
        A(i,i)=float(N)+A(i,i);
    enddo
    call disp("matrix A", A);
    C=A;
    ! !> vector version for chol decomposite
    ! do k=1,n
    !     A(k:n,k)=A(k:n,k)/sqrt(A(k,k));
    !     do j=k+1,n
    !         A(j:n,j)=A(j:n,j)-A(j:n,k)*A(j,k);
    !     enddo
    ! enddo
    ! > scalar version for chol decomposition
    do k=1,n
        mem=sqrt(A(k,k));
        do i=k,n
            A(i:n,k)=A(i,k)/sqrt(A(k,k));
        enddo
        do j=k+1,n
            do i=j,n
                A(i,j)=A(i,j)-A(i,k)*A(j,k);
            enddo
        enddo
    enddo

    ! do k=1,n
    !     A(k:n,k)=A(k:n,k)/sqrt(A(k,k));
    !     A(k+1:n,k+1)=(A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n));
    ! enddo
    do j=1,n
        G(j:n,j)=A(j:n,j)
    enddo
    print *, maxval(abs(C-matmul(G,transpose(G))))
    deallocate(A);
    deallocate(C);
    deallocate(G);
endprogram chol_decomposite_basic;