clc
clear all

ip=fopen('Yins_input.m','r++');
op=fopen('Yins_output.m','w++');

n=fscanf(ip,'%f',1);
m=fscanf(ip,'%f',1);

linedata=fscanf(ip,'%f',[7 m]);
linedata=linedata';

yshunt=fscanf(ip,'%f',[n 1]);
yshunt=complex(0.0,yshunt);

lp=linedata(:,1);
lq=linedata(:,2);
r=linedata(:,3);
x=linedata(:,4);
ycp=complex(0.0,linedata(:,5));
ycq=complex(0.0,linedata(:,6));
% if 0==0 
%     fprintf("hell yeh");
% end
    
tap=linedata(:,7);

for k=1:m
    yline(k) = 1/complex(r(k),x(k));
    if tap(k) ~= 1
        t1 = 1-(1/tap(k));
        t2 = -t1/tap(k);
        ycp(k)=t2*yline(k);
        ycq(k)=t1*yline(k);
        yline(k)=yline(k)/tap(k);
    end
end

for i=1:n
    for j=1:n
        ybus(i,j)=complex(0.0,0.0);
    end
end

for k=1:m
    p=lp(k);
    q=lq(k);
    ybus(p,p)=ybus(p,p)+yline(k)+ycp(k);
    ybus(q,q)=ybus(q,q)+yline(k)+ycp(k);
    ybus(p,q)=ybus(p,q)-yline(k);
    ybus(q,p)=ybus(p,q);
end

for i=1:n
    ybus(i,1)=ybus(i,1)+yshunt(i);
end

fprintf(op,'\tybus:\n');