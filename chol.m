n=10;
A=rand(n); A=0.5*(A+A'); 
A=A+n*eye(n); 
B=A;

function A = scalar_chol(A, n)
  for j=1:n
    for k=1:j-1
      A(j:n,j)=A(j:n,j)-A(j:n,k)*A(j,k);
    endfor
    A(j:n,j)= A(j:n,j)/sqrt(A(j,j));
  endfor 
endfunction

function A = vector_chol(A,n)
  for j=1:n
   A(j:n,j)=A(j:n,j)-A(j:n,1:j-1)*A'(1:j-1,j);
   A(j:n,j)= A(j:n,j)/sqrt(A(j,j));
  endfor
endfunction

for i:n
  A = vector_chol(A, n)
  G=tril(A,0); % R=chol(B,"lower");
  ans[i] = norm(B-G*G')
endfor