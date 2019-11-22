function foldedCM=to_foldedCM(a,b)
%-----------------------------------------folder
N=a+2;%耦合矩阵的阶数
M=b;%横向耦合矩阵
for i=1:N-3 %有几个数相加
    if mod(i,2)==1
        c=-1;
        for j=1:N-3+1-i
            l=N-(i+1)/2-j+1;
            nn=l-1;
            k=(i+1)/2;
            m=k;
            theta=atan(c*M(k,l)/M(m,nn));
            R=eye(N);
            R(nn,nn)=cos(theta);
            R(l,l)=R(nn,nn);
            R(nn,l)=-sin(theta);
            R(l,nn)=-R(nn,l);
            M=R*M*R';
        end
    else
        c=1;
        for j=1:N-3+1-i
            k=3+i/2+j-2;
            m=k+1;
            l=N-i/2+1;
            nn=l;
            theta=atan(c*M(k,l)/M(m,nn));
            R=eye(N);
            R(k,k)=cos(theta);
            R(m,m)=R(k,k);
            R(k,m)=-sin(theta);
            R(m,k)=-R(k,m);
            M=R*M*R';
        end
    end
end
foldedCM=M;
