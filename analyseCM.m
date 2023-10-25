%----------------------------画图验证耦合矩阵的正确性
function analyseCM(a,b,c,d,e,f)
n=a;%-----------------order of filter
Fstr=b;
Fsto=c;
M=d;
CF=e;
FBW=f;
R=zeros(n+2,n+2);
R(1,1)=1;
R(n+2,n+2)=1;
W=eye(n+2,n+2);
W(1,1)=0;
W(n+2,n+2)=0;
for k=1:1000
    wp(k)=((Fstr+(Fsto-Fstr)*k/1000)/CF-CF/(Fstr+(Fsto-Fstr)*k/1000))/FBW;
    Bd(k)=Fstr+(Fsto-Fstr)*k/1000;
    A=R.*(-1i)+W.*wp(k)+real(M);
    AP=inv(A);
    S21(k)=-2*1i*AP(n+2,1);
    S11(k)=1+2*1i*AP(1,1);
end

%画出带通滤波器的图
figure('name','wawa');
subplot(1,2,1);
dBS21=-20*log10(abs(S21));
dBS11=-20*log10(abs(S11));
plot(Bd,dBS21,'b',Bd,dBS11,'r','LineWidth',1.5);
set(gca,'YDir','reverse');
set(gca,'FontSize',12);
set(gca,'linewidth',1.2);
xlabel('频率(MHz)');
ylabel('抑制/回波损耗(dB)');
ylim([-2, 50]);
grid on;
title('带通幅度响应','fontsize',14);
subplot(1,2,2);
PS21=angle(S21).*180/pi;
PS11=angle(S11).*180/pi;
plot(Bd,PS21,'b',Bd,PS11,'r','LineWidth',1.5);
set(gca,'FontSize',12);
set(gca,'linewidth',1.2);
ylim([-180,180]);
grid on;
xlabel('频率(MHz)');
ylabel('相位(°)');
title('带通相位响应','fontsize',14);
