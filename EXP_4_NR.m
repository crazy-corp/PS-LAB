clc
clear all
ip=fopen('NRinput.m','r++');
op=fopen('NRoutput.m','w++');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Read BASIC DATA %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = fscanf(ip,'%d',1);
nline = fscanf(ip,'%d',1);
nslack = fscanf(ip,'%d',1);
itermax = fscanf(ip,'%d',1);
epsilon = fscanf(ip,'%f',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Read LINE DATA %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linedata = fscanf(ip,'%f',[8 nline]);
linedata = linedata';
lno = linedata(:,1);
lp = linedata(:,2);
lq = linedata(:,3);
r = linedata(:,4);
x = linedata(:,5);
ycp = complex(0,linedata(:,6));
ycq = complex(0,linedata(:,7));
tap = linedata(:,8);
shunt = fscanf(ip,'%f',[n 1]);
shunt = complex(0,shunt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Read BUS DATA - 1 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
busdata1 = fscanf(ip,'%f',[5 n]);
busdata1 = busdata1';
bno = busdata1(:,1);
pgen = busdata1(:,2);
qgen = busdata1(:,3);
pload = busdata1(:,4);
qload = busdata1(:,5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Read BUS DATA - 2 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
busdata2 = fscanf(ip,'%f',[5 n]);
busdata2 = busdata2';
bno1 = busdata2(:,1);
itype = busdata2(:,2);
vsp = busdata2(:,3);
qmin = busdata2(:,4);
qmax = busdata2(:,5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% SPARSITY VECTORS FORMATION %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n
 nlcount(i) = 0;
end
for k = 1:nline
 p = lp(k);
 q = lq(k);
 nlcount(p) = nlcount(p)+1;
 nlcount(q) = nlcount(q)+1;
end
itagf(1) = 1;
itagto(1) = nlcount(1);
for i = 2:n
 itagf(i) = itagto(i-1) + 1;
 itagto(i) = itagto(i-1) + nlcount(i);
end
for i=1:n
 ifill(i) = 0;
end
for k = 1:nline
 p = lp(k);
 q = lq(k);
 lpq = itagf(p) + ifill(p);
 lqp = itagf(q) + ifill(q);
 adjq(lpq) = q;
 adjl(lpq) = k;
 adjq(lqp) = p;
 adjl(lqp) = k;
 ifill(q) = ifill(q)+1;
 ifill(p) = ifill(p)+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% LINE ADMITTANCE CALCULATION %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:nline
 yline(k) = 1/complex(r(k),x(k));
end
for k = 1:nline
 if (tap(k)~=1)
 t1 = 1-(1/tap(k));
 t2 = -t1/tap(k);
 ycp(k) = t2*yline(k);
 ycq(k) = t1*yline(k);
 yline(k) = yline(k)/tap(k);
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% YBUS FORMATION Using Sparsity Tech %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% DIAGONAL ELEMENTS OF YBUS %%%%%%%%%%%%%%%%
for i = 1:n
 ypp(i) = complex(0.0,0.0);
end
for k = 1:nline
 p = lp(k);
 q = lq(k);
 ypp(p) = ypp(p) + yline(k) + ycp(k);
 ypp(q) = ypp(q) + yline(k) + ycq(k);
end
for i = 1:n
 ypp(i) = ypp(i) + shunt(i);
end
%%%%%%% OFF DIAGONAL ELEMENTS OF YBUS %%%%%%%%%%%%%%%%
for j = 1:2*nline
 k = adjl(j);
 ypq(j) = -yline(k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Calculate Pinj, Qinj, Qmininj and Qmaxinj %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n
 pinj(i) = pgen(i) - pload(i);
 qinj(i) = qgen(i) - qload(i);
 qmininj(i) = qmin(i) - qload(i);
 qmaxinj(i) = qmax(i)- qload(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Assume Flat Voltage Start %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E(nslack) = complex(vsp(1),0);
for i = 2:n
 E(i) = complex(1,0);
end
%%%% Find the Magnitude and Phase angle of Voltages %%%%%%%%%%%
for i = 1:n
 v(i) = abs(E(i));
 delta(i) = angle(E(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Iterative Process Start %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:itermax

 deltapmax = 0;
 deltaqmax = 0;

 %%%%% Calculation of Pcal and Qcal %%%%%%%

 for i = 1:n

 ical = ypp(i)*E(i);

 jstart = itagf(i);
 jstop = itagto(i);

 for j = jstart:jstop
 q = adjq(j);
 ical = ical + ypq(j)*E(q);
 end

 scal = E(i)*conj(ical);

 pcal(i) = real(scal);

 qcal(i) = imag(scal);

 end
 %%%%%% Calculation of DelP, DelQ and other madifications to
 %%%%%% Vector

 for i = 1:n

 delp(i) = pinj(i) - pcal(i);
 delq(i) = qinj(i) - qcal(i);


 delp(nslack) = 0;
 delq(nslack) = 0;

 if (itype(i) == 2)
 delq(i) = 0;
 end

 %%%% Find the maximum DelP and DelQ values %%%%%

 if (abs(delp(i)) > deltapmax)
 deltapmax = abs(delp(i));
 end

 if (abs(delq(i)) > deltaqmax)
 deltaqmax = abs(delq(i));
 end
 end

 ideltapmax(iter) = deltapmax;
 ideltaqmax(iter) = deltaqmax;

 %%%%% Checking for convergence Criteria %%%%%%%%

 if (deltapmax <= epsilon && deltaqmax <= epsilon)

 fprintf(op,'Problem converged in %d iterations\n\n',iter);
 break;

 end

 %%%%% Initilization of Jacobian Matrix %%%%%%%%%%%

 for i = 1:(2*n)
 for j = 1:(2*n)
 A(i,j) = 0.0;
 end
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%% Updation of Jacobian Matrix %%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%% Diagonal Elements of Jacobian %%%%%%%%

 for i = 1:n
 qp = qcal(i);
 bpp = imag(ypp(i));
 v2 = v(i)*v(i);
 pp = pcal(i);
 gpp = real(ypp(i));
 A(i,i) = -qp - bpp*v2; %% H - SUB MATRIX
 A(i,i+n) = pp + gpp*v2; %% N - SUB MATRIX
 A(i+n,i) = pp - gpp*v2; %% M - SUB MATRIX
 A(i+n,i+n) = qp - bpp*v2; %% L - SUBMATRIX
 end

 %%%%% Off - Diagonal Elements of Jacobian %%%%%%%%

 for p = 1:n

 jstart = itagf(p);
 jstop = itagto(p);

 for j = jstart:jstop

 q = adjq(j);

 ep = real(E(p));
 eq = real(E(q));
 fp = imag(E(p));
 fq = imag(E(q));

 term2 = fp*eq - ep*fq;
 term1 = ep*eq + fp*fq;

 gpq = real(ypq(j));

 bpq = imag(ypq(j));

 A(p,q) = gpq*term2 - bpq*term1; %% H - SUB MATRIX
 A(p,q+n) = gpq*term1 + bpq*term2; %% N - SUB MATRIX
 A(p+n,q) = -A(p,q+n); %% M - SUB MATRXI
 A(p+n,q+n) = A(p,q); %% L - SUB MATRIX

 end
 end

 %%%%% Setting the Slack bus elements to ZERO %%%%%%%%%

 A(nslack,nslack) = 10^20; %% H - SUB MATRIX
 A(n+nslack,n+nslack) = 10^20; %% L - SUB MATRIX

 %%%% Setting the PV bus elements to ZERO %%%%%%%%%

 for i = 1:n
 if (itype(i) == 2)
 A(i+n,i+n) = 10^20;
 end
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%% Setting to DelPQ Vector %%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 for i = 1:n
 delpq(i) = delp(i);
 delpq(i+n) = delq(i);
 end

 %%%%%%%%% Gauss Elimination Technique %%%%%%%%%%%%%%%
 %%%%%% [M] = [J] [X] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 for i = 1:(2*n)
 fact = (1/A(i,i));
 A(i,i) = 1;
 for j = (i+1):(2*n)
 A(i,j) = A(i,j)*fact;
 end
 delpq(i) = delpq(i)*fact;

 for k = (i+1):(2*n)
 fact1 = A(k,i);
 A(k,i) = 0;
 for j = (i+1):(2*n)
 A(k,j) = A(k,j) - A(i,j)*fact1;
 end
 delpq(k) = delpq(k) - fact1*delpq(i);
 end
 end

 delx(2*n) = delpq(2*n);
 
 for i = ((2*n)-1):-1:1
 sum = 0;
 for j = (i+1):(2*n)
 sum = sum+A(i,j)*delx(j);
 end
 delx(i) = delpq(i) - sum;
 end

 %%%%%%%%% Updation of Bus Voltages %%%%%%%%%%%

 for i = 1:n
 delta(i) = delta(i) + delx(i);
 v(i) = v(i) + delx(i+n)*v(i);
 enew(i) = v(i)*cos(delta(i));
 fnew(i) = v(i)*sin(delta(i));
 E(i) = complex(enew(i),fnew(i));
 end

 %%%%% If problem has not Converged then print the following %%%%%%%%%

 if (iter == itermax)
 fprintf(op,'Problem not converged');
 end

end
%%%%%%%%%%%%%%% Convergence Plot %%%%%%%%%%%%%%%%
figure(1)
plot(ideltapmax)
hold on
plot(ideltaqmax)
title('Convergence Plot')
xlabel('Number of Iterations')
ylabel('DeltaPmax & DeltaQmax')
legend('deltapmax','deltaqmax')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Calculation of Line Flows %%%%%%%%%%%%%%
%%%%%%%%%%%%% Calculation of Line LOSSES %%%%%%%%%%%%%
%%%%%%%%%%%%% Calculation of Slack Bus Power %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TSloss=0;
for k = 1:nline
 p = lp(k);
 q = lq(k);

 Ipq = E(p)*ycp(k) + (E(p) - E(q))*yline(k);
 Spq(k) = E(p)*conj(Ipq); %% Line Flow from P to Q

 Iqp = E(q)*ycq(k) + (E(q) - E(p))*yline(k);
 Sqp(k) = E(q)*conj(Iqp); %% Line Flow from Q to P

 Sloss(k) = Spq(k) + Sqp(k); %%% Line Loss
 TSloss = TSloss + Sloss(k);
end
%%%%%%%%%%%%% Calculation of Slack Bus Power %%%%%%%%%%
Jstart = itagf(nslack);
Jstop = itagto(nslack);
sum = complex(0,0);
for j = Jstart:Jstop
 k = adjl(j);
 sum = sum + Spq(k);
end
pinj(nslack) = real(sum);
qinj(nslack) = imag(sum);
pgen(nslack) = pinj(nslack) + pload(nslack);
qgen(nslack) = qinj(nslack) + qload(nslack);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%% PRINT THE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Print the Bus voltages
for i = 1:n
 fprintf(op,' %d %f %f\n',i,abs(E(i)),angle(E(i)));
end
fprintf(op,'\n\n');
%%% Print the Line flows
for i=1:nline
 fprintf(op,'Spq(%d)=\t',i);
 if imag(Spq(i))<0
 fprintf(op,'\t%f-i(%f)\t',real(Spq(i)),abs(imag(Spq(i))));
 else
 fprintf(op,'\t%f+i(%f)\t',real(Spq(i)),abs(imag(Spq(i))));
 end
 fprintf(op,'\n');
end
fprintf(op,'\n');
for i=1:nline
 fprintf(op,'Sqp(%d)=\t',i);
 if imag(Sqp(i))<0
 fprintf(op,'\t%f-i(%f)\t',real(Sqp(i)),abs(imag(Sqp(i))));
 else
 fprintf(op,'\t%f+i(%f)\t',real(Sqp(i)),abs(imag(Sqp(i))));
 end
 fprintf(op,'\n');
end
fprintf(op,'\n');
%%% Print the Line flows
for i = 1:nline
 fprintf(op,'Sloss(%d)=\t',i);
 if imag(Sloss(i))<0
 fprintf(op,'\t%f-i(%f)\t',real(Sloss(i)),abs(imag(Sloss(i))));
 else
 fprintf(op,'\t%f+i(%f)\t',real(Sloss(i)),abs(imag(Sloss(i))));
 end
 fprintf(op,'\n');
end
if imag(TSloss)<0
 fprintf(op,'TotalSloss=\t\t%f-i(%f)\t',real(TSloss),imag(TSloss));
else
 fprintf(op,'TotalSloss=\t\t%f+i(%f)\t',real(TSloss),imag(TSloss));
end
fprintf(op,'\n\n');
%%%%%%%% Slack Bus Power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(op,'Pslack = %f\t',pgen(nslack));
fprintf(op,'Qslack = %f\t',qgen(nslack));