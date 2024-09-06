function [output1,output2,output3,output4] = Corrupted_Zhang_BootAug(A0,n,k,z,m,p_value)
nk=round(n/k);B=zeros(k,k);
zk=n*z;
A=A0;
A(nk-zk+1:nk,:)=A0(2*nk-zk+1:2*nk,:);
for i=1:k
B(i,i)=sum(sum(A((i-1)*nk+1:i*nk,(i-1)*nk+1:i*nk)))/(nk*(nk-1));
end
for i=2:k
    for j=1:i-1
B(i,j)=sum(sum(A((i-1)*nk+1:i*nk,(j-1)*nk+1:j*nk)))/(nk*nk);
    end
end
diag_B=diag(diag(B));B=B+B';B=B-diag_B;dist=1e-6; B=B+dist*ones(k,k);
%% Aug
n_new=round(n/k/2);
sizeB=length(B);
In_p=max(diag(B));
diag_B=diag(diag(B)); diag_B=reshape(diag_B,sizeB^2,1);
vec_B=reshape(B,sizeB^2,1);
Out_p=sort((vec_B-diag_B)/2); %% most of between probability are 0
Out_p=Out_p(find(Out_p>min(Out_p),1));
%% added community
OA_new=zeros(n,n_new);
for i=1:k
    new_w1=binornd(1,Out_p,nk*n_new,1);
    OA_new((i-1)*nk+1:i*nk,1:n_new)=reshape(new_w1,nk,n_new);
end
new_w0=binornd(1,In_p,n_new*(n_new+1)/2,1); 
a1=length(new_w0);
a2=(-1+sqrt(1+8*a1))/2 ; % 计算向量对应的上三角矩阵的维数
new_w0a=zeros(a2,a2); % 以下把向量转化为上三角
for k1=1:a2
    for k2=k1:a2
        index=sum(a2:-1:a2-k1+2)+k2-k1+1;  
        new_w0a(k1,k2)=new_w0(index);      
    end
end
new_w0a=new_w0a+new_w0a'-diag(diag(new_w0a));
IA_new=new_w0a;
A_new=[A,OA_new;OA_new',IA_new];
A_new=A_new+diag(-diag(A_new));
%% new hat_B
B_new=zeros(k+1,k+1);
for i=1:k
B_new(i,i)=sum(sum(A_new((i-1)*nk+1:i*nk,(i-1)*nk+1:i*nk)))/(nk*(nk-1));
end
B_new(k+1,k+1)=sum(sum(A_new(k*nk+1:k*nk+n_new,k*nk+1:k*nk+n_new)))/(n_new*(n_new-1));
for i=2:k
    for j=1:i-1
B_new(i,j)=sum(sum(A_new((i-1)*nk+1:i*nk,(j-1)*nk+1:j*nk)))/(nk*nk);
    end
end
for j=1:k
B_new(k+1,j)=sum(sum(A_new(k*nk+1:k*nk+n_new,(j-1)*nk+1:j*nk)))/(nk*n_new);
end
diag_B_new=diag(diag(B_new));
B_new=B_new+B_new';
B_new=B_new-diag_B_new;
dist=1e-6; B_new=B_new+dist*ones(k+1,k+1);
%% added rho
rho_new=zeros(n+n_new,k+1);
for i=1:n+n_new
    i1=ceil(i/nk);
    for v=1:k
        if i1==v   
        rho_new(i,v)=sum((A_new(i,setdiff((v-1)*nk+1:v*nk,i))-B_new(i1,v))./sqrt(B_new(i1,v)*(1-B_new(i1,v))))/sqrt(nk-1);
        elseif i1<=k
        rho_new(i,v)=sum((A_new(i,(v-1)*nk+1:v*nk)-B_new(i1,v))./sqrt(B_new(i1,v)*(1-B_new(i1,v))))/sqrt(nk);
        else
        rho_new(i,k+1)=sum((A_new(i,setdiff(k*nk+1:k*nk+n_new,i))-B_new(i1,k+1))./sqrt(B_new(i1,k+1)*(1-B_new(i1,k+1))))/sqrt(n_new-1);
        end
    end  
end
%% Zhang AUG Test Statistics 
S_Ln_new=max(max(abs(rho_new)));
S_Test_aug=S_Ln_new^2-2*log(2*(n+n_new)*(k+1))+log(log(2*(n+n_new)*(k+1)));
%% Aug Boot
S_Test=zeros(m,1);
for loop=1:m
Ak=Generating(n,k,B);
%% added community
OA_new=zeros(n,n_new);
for i=1:k
    new_w1=binornd(1,Out_p,nk*n_new,1);
    OA_new((i-1)*nk+1:i*nk,1:n_new)=reshape(new_w1,nk,n_new);
end
new_w0=binornd(1,In_p,n_new*(n_new+1)/2,1); 
a1=length(new_w0);
a2=(-1+sqrt(1+8*a1))/2 ; % 计算向量对应的上三角矩阵的维数
new_w0a=zeros(a2,a2); % 以下把向量转化为上三角
for k1=1:a2
    for k2=k1:a2
        index=sum(a2:-1:a2-k1+2)+k2-k1+1;  
        new_w0a(k1,k2)=new_w0(index);      
    end
end
new_w0a=new_w0a+new_w0a'-diag(diag(new_w0a));
IA_new=new_w0a;
A_new=[Ak,OA_new;OA_new',IA_new];
A_new=A_new+diag(-diag(A_new));
%% new hat_B
B_new=zeros(k+1,k+1);
for i=1:k
B_new(i,i)=sum(sum(A_new((i-1)*nk+1:i*nk,(i-1)*nk+1:i*nk)))/(nk*(nk-1));
end
B_new(k+1,k+1)=sum(sum(A_new(k*nk+1:k*nk+n_new,k*nk+1:k*nk+n_new)))/(n_new*(n_new-1));
for i=2:k
    for j=1:i-1
B_new(i,j)=sum(sum(A_new((i-1)*nk+1:i*nk,(j-1)*nk+1:j*nk)))/(nk*nk);
    end
end
for j=1:k
B_new(k+1,j)=sum(sum(A_new(k*nk+1:k*nk+n_new,(j-1)*nk+1:j*nk)))/(nk*n_new);
end
diag_B_new=diag(diag(B_new));
B_new=B_new+B_new';
B_new=B_new-diag_B_new;
dist=1e-6; B_new=B_new+dist*ones(k+1,k+1);
%% added rho
S_rho_new=zeros(n+n_new,k+1);
for i=1:n+n_new
    i1=ceil(i/nk);
    for v=1:k
        if i1==v   
        S_rho_new(i,v)=sum((A_new(i,setdiff((v-1)*nk+1:v*nk,i))-B_new(i1,v))./sqrt(B_new(i1,v)*(1-B_new(i1,v))))/sqrt(nk-1);
        elseif i1<=k
        S_rho_new(i,v)=sum((A_new(i,(v-1)*nk+1:v*nk)-B_new(i1,v))./sqrt(B_new(i1,v)*(1-B_new(i1,v))))/sqrt(nk);
        else
        S_rho_new(i,k+1)=sum((A_new(i,setdiff(k*nk+1:k*nk+n_new,i))-B_new(i1,k+1))./sqrt(B_new(i1,k+1)*(1-B_new(i1,k+1))))/sqrt(n_new-1);
        end
    end  
end
%% AUG Test Statistics 
S_Ln=max(max(abs(S_rho_new)));
S_Test(loop,:)=S_Ln^2-2*log(2*(n+n_new)*(k+1))+log(log(2*(n+n_new)*(k+1)));
end
param0=[-2*log(sqrt(pi)),1];
options = optimset('GradObj','off'); %%检察参数估计
fun2= {@(param)(m*log(param(2))-sum(-(S_Test-param(1))./param(2))+sum(exp(-(S_Test-param(1))./param(2))))};
[param_mleS,~,~,~]=fminunc(fun2,param0,options);
hat_muS=param_mleS(1);hat_betaS=param_mleS(2);
Test_SBoot=-2*log(2*sqrt(pi))+2*(S_Test_aug-hat_muS)/hat_betaS;
%%
ev=-2*log(-2*sqrt(pi)*log(1-p_value));
if Test_SBoot>ev
    Size_S=1;
    else
    Size_S=0;
end
output1=Size_S;
output2=Test_SBoot;
output3=rho_new;
output4=B;

