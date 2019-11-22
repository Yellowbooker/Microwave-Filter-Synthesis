function [CTCM,flag]=to_CTCM(a,b,c)
N=a;%------order of filter
Tz=b./1i;%---------location of Tz
nfz=length(b);
M=c;%----------N+2orderWheelCM
%--------------------------------CT,afterWheel
%以前的N为矩阵的阶数，现在为滤波器的阶数

if (N/(nfz*2))==1
    NumTz=nfz-1;
    flag=1;
elseif (N/(nfz*2))>1
    NumTz=nfz;
    flag=1;
else
    disp('零点过多！');
    NumTz=0;
    flag=0;
end

fprintf('存在 %d 个完全三元组（CT）\n',NumTz);    
for i=1:NumTz
            theta=atan(M(N,N+1)/(M(N+1,N+1)+Tz(1,i)));
            R=eye(N+2);
            R(N,N)=cos(theta);
            R(N+1,N+1)=R(N,N);
            R(N,N+1)=-sin(theta);
            R(N+1,N)=-R(N,N+1);
            M=R*M*R';
    for j=1:N-1-2*i
            k=N-j-1;
            m=k+1;
            l=m+1;
            n=l;
            theta=atan(M(k+1,l+1)/(M(m+1,n+1)));
            R=eye(N+2);
            R(k+1,k+1)=cos(theta);
            R(m+1,m+1)=R(k+1,k+1);
            R(k+1,m+1)=-sin(theta);
            R(m+1,k+1)=-R(k+1,m+1);
            M=R*M*R';
    end
end
CTCM=M;
