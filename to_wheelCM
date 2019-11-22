function wheelCM=to_wheelCM(a,b)
%------------------------------wheel
N=a+2;%耦合矩阵的阶数
M=b;%folded耦合矩阵
for i=1:N-3 %有几个数相加
        for j=1:N-3-i+1
            l=i-1+j+2;
            k=i;
            n=k+1;
            m=k;
            theta=-atan(M(k,l)/M(m,n));
            R=eye(N);
            R(n,n)=cos(theta);
            R(l,l)=R(n,n);
            R(n,l)=-sin(theta);
            R(l,n)=-R(n,l);
            M=R*M*R';
        end
end
wheelCM=M;
