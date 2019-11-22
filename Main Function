%广义切比雪夫滤波器N+2阶耦合矩阵综合计算程序
%2019/08/30
%2019/10/15
%2019/10/16修改三元组提取
%2019/10/18四元组提取初步
%参考文献--------------------------------------------------------------------------------：
%1.Microwave Filters for Communication Systems Fundamentals, Design, and Applications 2nd
%2.General coupling matrix synthesis methods for Chebyshev filtering functions
%3.Advanced Coupling Matrix Synthesis Techniques for Microwave Filters
%4.Adaptive Synthesis and Design of Resonator Filters With Source/Load-Multiresonator Coupling
%5.An Analytical Technique for the Synthesis of Cascaded N-Tuplets Cross-Coupled Resonators Microwave Filters Using Matrix Rotations
%含有耦合矩阵Floder规范结构、wheel规范结构、CT结构、CQ结构及横向规范结构
%!小于三阶的全规范滤波函数综合时会出错，阶数过高（大于20）也可能会出错
%!2019/06/06，群时延计算错误
%!CT、CQ拓扑复数零点计算有误2019/10/29
clc;
clear all;

% eg1:
% n=4;
% S=[-3.7431j,-1.8051j,1.5699j,6.1910j];
% RL=22;

% eg2:
% n=4;
% S=[1.3217j,1.8082j];
% RL=22;

% eg3:
% n=7;
% S=[1.3958j,-1.3958j,1.0749,-1.0749];
% RL=23;

% n=input('输入滤波器阶数\nn=') %滤波器的阶数
% S=input('输入有限位置传输零点\nTz=') %有限位置传输零点位置，以向量形式输入归一化低通原型频率对应的传输零点，eg：S=[2j，-2j]
% RL=input('输入带内回波损耗\nRL=') %带内回波损耗的最小值，与波纹系数有关
% Start=input('起始频率\nStart Freq=') %MHz
% Stop=input('终止频率\nStop Freq=') %MHz
% Q=input('谐振器无载Q值\nQu=')
% Fstr=input('输入显示起始频率\nFstr=') %MHz
% Fsto=input('输入显示截至频率\nFsto=') %MHz

n=6;
S=[-2j,2j];
nfz=length(S);
RL=20;
Start=3450;
Stop=3550;
Q=200;

Fstr=3000;
Fsto=4000;

BW=Stop-Start;
CF=sqrt(Start*Stop);
FBW=BW/CF;

w=S*-1i;%复频率转化至实频率

%第一步：广义切比雪夫函数综合
%根据有限位置传输零点确定特征多项式P(s)/P(w)
P=1;
for k=1:length(w)
    P=conv(P,[1,-w(1,k)]);
end
P1=P.*1i.^[0:1:length(w)];
if mod((n-length(w)),2)==0
    P1=P1.*1i;
end
%P(s)/P(w)综合结束

%根据广义切比雪夫滤波函数以及P确定特征多项式F(s)/F(w)
w=[w,inf(1,n-length(w))];
U1=[1,-1/w(1,1)];
V1=[sqrt(1-1/(w(1,1)*w(1,1)))];
U2=U1;
V2=V1;
for k=2:n
    U2=conv(U1,[1,-1/w(1,k)])+conv([1,0,-1],V1.*sqrt(1-1/(w(1,k)*w(1,k))));
    V2=conv(V1,[1,-1/w(1,k)])+U1.*sqrt(1-1/(w(1,k)*w(1,k)));
    U1=U2;
    V1=V2;
end
F=U1./(U1(1,1));
V=V1./(V1(1,1));
V1=V.*1i.^[0:1:length(V)-1];
F1=F.*1i.^[0:1:n];
%F(s)/F(w)综合结束

%根据能量守恒公式以及P、F确定特征多项式E(s)/E(w)
P2=conv(P,P);
F2=conv(F,F);
P2=[zeros(1,length(F2)-length(P2)),P2];
ep=1/sqrt(10^(RL/10)-1)*abs(polyval(P,1)/polyval(F,1));

if length(S)==n
    epr=ep/sqrt(ep*ep-1);
else
    epr=1;
end
E2=P2./(ep*ep)+F2./(epr*epr);
RootE=roots(E2);
leftRootE=find(imag(RootE)>0);
RootE=RootE(leftRootE);
E=poly(RootE);
E1=E.*1i.^[0:1:n];
%E(s)/E(w)综合结束

for k=1:1000
    d(k)=-8+16/1000*k;
    r(k)=polyval(F,d(k))/polyval(E,d(k))/epr;
    t(k)=polyval(P,d(k))/polyval(E,d(k))/ep;
end

figure('name','无耗综合结果');
subplot(2,3,1);
I=-20*log10(abs(t));
R=-20*log10(abs(r));
plot(d,I,'b',d,R,'r','LineWidth',1.5);
set(gca,'YDir','reverse');
set(gca,'FontSize',12);
set(gca,'linewidth',1.2);
xlabel('归一化频率(rad/s)');
ylabel('抑制/回波损耗(dB)');
ylim([-2, 50]);
xlim([-8, 8]);
grid on;
title('低通原型幅度响应','fontsize',14);

%F3：反射零点，E3：反射/传输极点，P3：有限位置传输零点，V3：带内反射极大值
F3=roots(F1);
E3=roots(E1);
P3=roots(P1);
V3=roots(V1);
realF=real(F3);
realE=real(E3);
realP=real(P3);
realV=real(V3);
imagF=imag(F3);
imagE=imag(E3);
imagP=imag(P3);
imagV=imag(V3);
subplot(2,3,2);
plot(realF,imagF,'o',realP,imagP,'^','LineWidth',1.5);
axis([-2,2,-8,8]);
xlabel('实部');
ylabel('虚部');
set(gca,'FontSize',12);
set(gca,'linewidth',1.2);
title('反射零点与传输零点','fontsize',14);
grid on;
subplot(2,3,3);
plot(realE,imagE,'s',realV,imagV,'v','LineWidth',1.5);
axis([-2,2,-1.5,1.5]);
xlabel('实部');
ylabel('虚部');
set(gca,'FontSize',12);
set(gca,'linewidth',1.2);
title('传输/反射极点与带内反射极大值点','fontsize',14);
grid on;

%第二步：N+2阶耦合矩阵综合
%由S参数计算导纳多项式
M1=zeros(1,n);
N1=zeros(1,n);
if mod(n,2)==0
    for k=n+1:-2:1
        M1(1,k)=real(E1(1,k)+F1(1,k));
        N1(1,k)=1i*imag(E1(1,k)+F1(1,k));
    end
    for k=n:-2:2
        M1(1,k)=1i*imag(E1(1,k)+F1(1,k));
        N1(1,k)=real(E1(1,k)+F1(1,k));
    end
    lambdaK=roots(M1);
else
    for k=n+1:-2:2
        M1(1,k)=real(E1(1,k)+F1(1,k));
        N1(1,k)=1i*imag(E1(1,k)+F1(1,k));
    end
    for k=n:-2:1
        M1(1,k)=1i*imag(E1(1,k)+F1(1,k));
        N1(1,k)=real(E1(1,k)+F1(1,k));
    end
    lambdaK=roots(N1);
    T1=N1;
    N1=M1;
    M1=T1;
end

%计算MSL
Msl=sqrt((epr-1)/(epr+1));
if Msl==0
    MP=P1/ep;
else
    MP=P1/ep-1i*Msl*M1;
end

%根据特征值计算留数
for k=1:length(lambdaK)
    r22(k)=polyval(N1,lambdaK(k))/polyval(polyder(M1),lambdaK(k));
    r21(k)=polyval(MP,lambdaK(k))/polyval(polyder(M1),lambdaK(k));
end

%根据特征值与留数计算N+2阶耦合矩阵
M=zeros(n+2,n+2);
for k=1:n
    M(k+1,1)=sqrt(r22(k));
    M(1,k+1)=M(k+1,1);
    M(k+1,n+2)=r21(k)/M(k+1,1);
    M(n+2,k+1)=M(k+1,n+2);
    M(k+1,k+1)=lambdaK(k)*1i;
end
M(1,n+2)=Msl;
M(n+2,1)=Msl;
disp('横向拓扑：');
display(real(M));
%N+2阶未化简耦合矩阵计算结束

%第三步：由耦合矩阵计算带通滤波器S参数(N+2阶耦合矩阵的直接分析)
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
subplot(2,3,4);
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
subplot(2,3,5);
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

%滤波器群时延
tao21=-1*diff(PS21)/((Fstr-Fsto)/k);%%%%%%
tao21=[tao21,0];
subplot(2,3,6);
plot(Bd,tao21,'.b','LineWidth',1.5);
set(gca,'FontSize',12);
set(gca,'linewidth',1.2);
ylim([-0.5,12]);
grid on;
xlabel('频率(MHz)');
ylabel('群时延(ns)');
title('带通滤波器群时延','fontsize',14);
cursorMode=datacursormode(gcf);
set(cursorMode,'enable','on');

%包含损耗的计算
figure('name','低损耗综合结果');
subplot(1,3,1);
for k=1:1000
    wpp(k)=((Fstr+(Fsto-Fstr)*k/1000)/CF-CF/(Fstr+(Fsto-Fstr)*k/1000))/FBW-1i*1/Q/FBW;
    Bd(k)=Fstr+(Fsto-Fstr)*k/1000;
    Ap=R.*(-1i)+W.*wpp(k)+real(M);
    APp=inv(Ap);
    S21p(k)=-2*1i*APp(n+2,1);
    S11p(k)=1+2*1i*APp(1,1);
end
dBS21p=-20*log10(abs(S21p));
dBS11p=-20*log10(abs(S11p));
plot(Bd,dBS21p,'b',Bd,dBS11p,'r','LineWidth',1.5);
set(gca,'YDir','reverse');
set(gca,'FontSize',12);
set(gca,'linewidth',1.2);
xlabel('频率(MHz)');
ylabel('抑制/回波损耗(dB)');
ylim([-2, 50]);
grid on;
title('带损耗带通幅度响应','fontsize',14);
subplot(1,3,2);
plot(Bd,dBS21,'b',Bd,dBS21p,'r','LineWidth',1.5);
set(gca,'YDir','reverse');
set(gca,'FontSize',12);
set(gca,'linewidth',1.2);
xlabel('频率(MHz)');
ylabel('抑制/回波损耗(dB)');
ylim([0, 1]);
grid on;
title('有/无损耗带通幅度响应对比(窄带)','fontsize',14);

subplot(1,3,3);
plot(Bd,dBS21,'b--',Bd,dBS21p,'r','LineWidth',1.5);
set(gca,'YDir','reverse');
set(gca,'FontSize',12);
set(gca,'linewidth',1.2);
xlabel('频率(MHz)');
ylabel('抑制/回波损耗(dB)');
axis([3000,4000,-2,50]);
grid on;
title('有/无损耗带通幅度响应对比(宽带)','fontsize',14);
cursorMode=datacursormode(gcf);
set(cursorMode,'enable','on');

%------------------------folded coupling matrix
foldedCM=to_foldedCM(n,M);
disp('规范折叠拓扑：');
display(real(foldedCM));

%------------------------wheel coupling matrix
wheelCM=to_wheelCM(n,foldedCM);
disp('规范轮形拓扑：');
display(real(wheelCM));

%-----------------------CT coupling matrix
[CTCM,flag]=to_CTCM(n,S,wheelCM);
disp('三元组（CT）拓扑：');
display(real(CTCM));

%-----------------------CQ coupling matrix
CQCM=to_CQCM(flag,n,CTCM,nfz);
disp('四元组（CQ）拓扑：');
display(real(CQCM));

%-----------------------draw CM
analyseCM(n,Fstr,Fsto,CQCM,CF,FBW);

%-----------------------反归一化耦合矩阵（输出耦合带宽以及谐振频率）
unnormalizedCM=un_normalizeCM(CQCM,CF,FBW);
disp('反归一化耦合矩阵(MHz)：');
display(real(unnormalizedCM));
