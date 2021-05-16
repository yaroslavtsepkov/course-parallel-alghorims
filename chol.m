clear all;


# Функция генерирует и проводит матрицу к нужному нам виду
function A = preprocessing(n)
  A = rand(n);
  A = 0.5*(A+A');
  A = A + n*eye(n);
end
# Матричное разложение холецкого
function A = matrix_chol(A,n)
  for k=1:n
    A(k:n,k)=A(k:n,k)/sqrt(A(k,k));
    A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*(A(k+1:n,k))';
  end
end

# Строчный алгоритм разложения Холецкого
function A = row_chol(A,n)
  for k=1:n
    A(k:n,k)=A(k:n,k)/sqrt(A(k,k));
    for i=k+1:n
      A(i,k+1:i)=A(i,k+1:i)-A(i,k)*(A(k+1:i,k))';
    end
  end
end 

# Векторное разложение холецкого
function A = vector_chol(A,n)
  for k=1:n
    A(k:n,k)=A(k:n,k)/sqrt(A(k,k));
    for j=k+1:n
      A(j:n,j)=A(j:n,j)-A(j:n,k)*A(j,k);
    end
  end
end

# Функция скалярного разложения 
function A = scalar_chol(A, n)
  for k=1:n
    temp = sqrt(A(k,k));
    for i=k:n
      A(i,k)=A(i,k)/temp;
    end
    for j=k+1:n
      for i=j:n
        A(i,j)=A(i,j)-A(i,k)*A(j*k);
      end
    end
  end
end

# Функция вычисляет ошибку разложения
function error = calcError(A, B)
  G=tril(A,0);
  error = norm(B-G*G');
end

#minmat = 100;
#step=50;
#maxmat = 1000;
#for n=minmat:step:maxmat
n=30;
A = preprocessing(n);
B = A;
start = time();
A = scalar_chol(A,n);
finish = time() - start;
ans = calcError(A, B);

#end
