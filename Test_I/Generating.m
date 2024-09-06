function output = Generating(n,k,B)  %pw0,pw1
A=zeros(n,n); nk=n/k; 
for i=1:k
    %w0=binornd(1,pw0,nk*nk,1);   
    w0=binornd(1,B(i,i),nk*(nk+1)/2,1); 
    a1=length(w0);
    a2=(-1+sqrt(1+8*a1))/2 ; % 计算向量对应的上三角矩阵的维数
    w0a=zeros(a2,a2); % 以下把向量转化为上三角
for k1=1:a2
    for k2=k1:a2
        index=sum(a2:-1:a2-k1+2)+k2-k1+1;  
        w0a(k1,k2)=w0(index);      
    end
end
     w0a=w0a+w0a';
     w0a=w0a-diag(diag(w0a));
     A((i-1)*nk+1:i*nk,(i-1)*nk+1:i*nk)=w0a;
%A((i-1)*nk+1:i*nk,(i-1)*nk+1:i*nk)=reshape(w0,nk,nk);
end
for i=2:k
    for j=1:i-1
    %w1=binornd(1,pw1,nk*nk,1);
    w1=binornd(1,B(i,j),nk*(nk+1)/2,1);
    a1=length(w1);
    a2=(-1+sqrt(1+8*a1))/2 ; % 计算向量对应的上三角矩阵的维数
    w1a=zeros(a2,a2); % 以下把向量转化为上三角
for k1=1:a2
    for k2=k1:a2
        index=sum(a2:-1:a2-k1+2)+k2-k1+1;  
        w1a(k1,k2)=w1(index);      
    end
end
    w1a=w1a+w1a';
    w1a=w1a-diag(diag(w1a));
    A((i-1)*nk+1:i*nk,(j-1)*nk+1:j*nk)=w1a; 
    %A((i-1)*nk+1:i*nk,(j-1)*nk+1:j*nk)=reshape(w1,nk,nk); 
    end
end   
%%对称
for i=1:n
    for j=1:i-1
        A(j,i)=A(i,j);
    end
end
%A=sparse(A);
%dist=1e-6; A=A+dist*ones(n,n); 
A=A+diag(-diag(A));
output=A;

