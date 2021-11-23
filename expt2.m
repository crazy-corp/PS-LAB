clc
clear all

ip=fopen('Yins_input.m','r++');
op=fopen('Yins_output.m','w++');

nbus=fscanf(ip,'%d',1);
nline=fscanf(ip,'%d',1);

a=fscanf(ip,'%g',[8,nline]);
A=a';

LNo=A(:,1);
LP=A(:,2);
LQ=A(:,3);
R=A(:,4);
X=A(:,5);
Ycp=complex(0,A(:,6));
Ycq=complex(0,A(:,7));
Tap=A(:,8);

yshunt=fscanf(ip,'%f',[nbus 1]);
yshunt=complex(0.0,yshunt);


for i = 1:nbus
    NLCONT(i) = 0;
end

%Updation of NLCONT(I) Vector

for k = 1:nline
    P = LP(k);
    Q = LQ(k);
    NLCONT(P) = NLCONT(P) + 1;
    NLCONT(Q) = NLCONT(Q) + 1;  
end

ITAGF(1) = 1;
ITAGTO(1) = NLCONT(1);

for i = 2:nbus
    ITAGF(i) = ITAGTO(i-1) + 1;
    ITAGTO(i) = ITAGTO(i-1) + NLCONT(i);
end

for i = 1:nbus
    IFILL(i) = 0;
end

for k = 1:nline
    P = LP(k);
    Q = LQ(k);
    LPQ = ITAGF(P) +IFILL(P);
    LQP = ITAGF(Q) +IFILL(Q);
    ADJQ(LPQ) = Q;
    ADJL(LPQ) = k;
    ADJQ(LQP) = P;
    ADJL(LQP) = k;
    IFILL(P) = IFILL(P) + 1;
    IFILL(Q) = IFILL(Q) + 1;
end

%calculate admittances of lines

for k=1:nline
    yline(k) = 1/complex(R(k),X(k));
    if Tap(k) ~= 1
        t1 = 1-(1/Tap(k));
        t2 = -t1/Tap(k);
        Ycp(k)=t2*yline(k);
        Ycq(k)=t1*yline(k);
        yline(k)=yline(k)/Tap(k);
    end
end

%Formation of Ybus Matrix

%Diagonal elements

for i=1:nbus
    YPP(i)=complex(0,0);
end

for k=1:nline
    P=LP(k);
    Q=LQ(k);
    YPP(P)=YPP(P)+yline(k)+Ycp(k);
    YPP(Q)=YPP(Q)+yline(k)+Ycp(k);
end

for i=1:nbus
    YPP(i)=YPP(i)+yshunt(i);
end

%Off diagonal elements

for j=1:2*nline
    k=ADJL(j);
    Ybus(j)= -yline(k);
end