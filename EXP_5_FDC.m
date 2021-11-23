clc
clear all
ip=fopen('NRinput.m','r++');
op=fopen('FDCoutput.m','w++');
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
%%%%%%% Compensation %%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Formation of B1 Matrix %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initilization of B1 matrix to Zeros
for i = 1:n
 for j = 1:n
 B1(i,j) = 0.0;
 end
end
%%%% Updation of B1 matrix
for k = 1:nline
 bline(k) = 1/x(k);
end
for k = 1:nline
 p = lp(k);
 q = lq(k);
 B1(p,p) = B1(p,p) + bline(k);
 B1(q,q) = B1(q,q) + bline(k);
 B1(p,q) = B1(p,q) - bline(k);
 B1(q,p) = B1(q,p) - bline(k);
end
B1(nslack,nslack) = 10^20; %%% Slack Bus treatment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Chelosky Method for B1 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B1(1,1) = sqrt(B1(1,1));
for j = 2:n
 B1(1,j) = B1(1,j)/B1(1,1);
 B1(j,1) = B1(1,j);
end
for i = 2:n
 sum = 0;
 for k = 1:i-1
 sum = sum + B1(i,k)*B1(i,k);
 end
 B1(i,i) = sqrt(B1(i,i)-sum);
 for j = i+1:n
 sum = 0;
 for k = 1:i-1
 sum = sum + B1(i,k)*B1(k,j);
 end
 B1(i,j) = (B1(i,j)-sum)/B1(i,i);
 B1(j,i) = B1(i,j);
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Formation of B2 Matrix %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initilization of B2 matrix to Zeros
for i = 1:n
 for j = 1:n
 B2(i,j) = 0.0;
 end
end
%%%%% Updation of B2 matrix
for p = 1:n
 B2(p,p) = -imag(ypp(p));
end
for p = 1:n

 jstart = itagf(p);
 jstop = itagto(p);

 for j = jstart:jstop
 q = adjq(j);
 B2(p,q) = B2(p,q) - imag(ypq(j));
 end
end
B2(nslack,nslack) = 10^20; %%% Slack Bus treatment
for p = 1:n
 if (itype(p) == 2)
 B2(p,p) = 10^20; %%%% PV bus treatment
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Chelosky Method for B2 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B2(1,1) = sqrt(B2(1,1));
for j = 2:n
 B2(1,j) = B2(1,j)/B2(1,1);
 B2(j,1) = B2(1,j);
end
for i = 2:n
 sum = 0;
 for k = 1:i-1
 sum = sum + B2(i,k)*B2(i,k);
 end
 B2(i,i) = sqrt(B2(i,i)-sum);
 for j = i+1:n
 sum = 0;
 for k = 1:i-1
 sum = sum + B2(i,k)*B2(k,j);
 end
 B2(i,j) = (B2(i,j)-sum)/B2(i,i);
 B2(j,i) = B2(i,j);
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Assume Flat Voltage Start %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E(nslack) = complex(vsp(1),0); %%%% For Slack bus
for i = 2:n
 E(i) = complex(1,0); % For all other buses
end
%%%% Find the Magnitude and Phase angle of Bus Voltages %%%%%%%%%%%
for i = 1:n
 v(i) = abs(E(i));
 delta(i) = angle(E(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Iterative Process Start %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 0:itermax


 %%%%%%%%% One delta operation %%%%%%%%%%%%%


 deltapmaxd = 0;
 deltaqmaxd = 0;

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

 %%%%%% Calculation of DelP, DelQ and other madifications to Mismatch-
 %%%%%% Vector

 for i = 1:n

 delp(i) = pinj(i) - pcal(i);
 delq(i) = qinj(i) - qcal(i);


 delp(nslack) = 0; %% For Slack bus
 delq(nslack) = 0; %% For Slack bus

 if (itype(i) == 2)
 delq(i) = 0; %% For PV bus
 end

 %%%% Find the maximum DelP and DelQ values %%%%%

 if (abs(delp(i)) > deltapmaxd)
 deltapmaxd = abs(delp(i));
 end

 if (abs(delq(i)) > deltaqmaxd)
 deltaqmaxd = abs(delq(i));
 end
 end

 iter
 deltapmaxd
 deltaqmaxd

 %%%%% Checking for convergence Criteria %%%%%%%%

 if (deltapmaxd <= epsilon && deltaqmaxd <= epsilon)

 fprintf(op,'\n\nProblem converged in %d iterations\n\n',iter);
 fprintf(op,'\n');
 break;

 end

 %%%%%%%%% Solve for Delta %%%%%%%%%%%%%

 for i = 1:n
 delp1(i) = delp(i)/v(i);
 end

 del1(1) = delp(1)/B1(1,1);

 for i = 2:n
 sum = 0;
 for j = 1:i-1
 sum = sum + B1(i,j)*del1(j);
 end
 del1(i) = (delp(i)-sum)/B1(i,i);
 end

 delx(n) = del1(n)/B1(n,n);

 for i = n-1:-1:1
 sum = 0;
 for j = i+1:n
 sum = sum + B1(i,j)*delx(j);
 end
 delx(i) = (del1(i)-sum)/B1(i,i);
 end

 %%%%%%% Delta & Voltage Updation %%%%%%%%%%%

 for i = 1:n
 delta(i) = delta(i) + delx(i);
 e1(i) = v(i)*cos(delta(i));
 f1(i) = v(i)*sin(delta(i));
 E(i) = complex(e1(i),f1(i));
 end



 %%%%%%% First Half Iteration is over %%%%%%%


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 iter = iter + 0.5;

 %%%%%%%%% One V operation %%%%%%%%%%%%%

 deltapmaxv = 0;
 deltaqmaxv = 0;

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

 %%%%%% Calculation of DelP, DelQ and other madifications tom Mismatch-
 %%%%%% Vector

 for i = 1:n

 delp(i) = pinj(i) - pcal(i);
 delq(i) = qinj(i) - qcal(i);


 delp(nslack) = 0; %% For Slack Bus
 delq(nslack) = 0; %% For Slack Bus

 if (itype(i) == 2)
 delq(i) = 0; %% For PV Bus
 end

 %%%% Find the maximum DelP and DelQ values %%%%%

 if (abs(delp(i)) > deltapmaxv)
 deltapmaxv = abs(delp(i));
 end

 if (abs(delq(i)) > deltaqmaxv)
 deltaqmaxv = abs(delq(i));
 end
 end

 iter
 deltapmaxv
 deltaqmaxv

 %%%%% Checking for convergence Criteria %%%%%%%%

 if (deltapmaxv <= epsilon && deltaqmaxv <= epsilon)

 fprintf(op,'\n\nProblem converged in %d iterations\n\n',iter);
 fprintf(op,'\n');
 break;

 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%% Solve the Equation for DEL V %%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 for i = 1:n
 delq1(i) = delq(i)/v(i);
 end

 del2(1) = delq(1)/B2(1,1);

 for i = 2:n
 sum = 0;
 for j = 1:i-1
 sum = sum + B2(i,j)*del2(j);
 end
 del2(i) = (delq(i)-sum)/B2(i,i);
 end

 dely(n) = del2(n)/B2(n,n);

 for i = n-1:-1:1
 sum = 0;
 for j = i+1:n
 sum = sum + B2(i,j)*dely(j);
 end
 dely(i) = (del2(i)-sum)/B2(i,i);
 end

 %%%%%%%%%% Voltage updation %%%%%%%%%%5

 for i = 1:n
 v(i) = v(i) + dely(i);
 e2(i) = v(i)*cos(delta(i));
 f2(i) = v(i)*sin(delta(i));
 E(i) = complex(e2(i),f2(i));
 end

 if (iter == itermax)
 fprintf(op,'Problem not converged in %d iterations',itermax);
 break;
 end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Calculation of Line Flows %%%%%%%%%%%%%%
%%%%%%%%%%%%% Calculation of Line LOSSES %%%%%%%%%%%%%
%%%%%%%%%%%%% Calculation of Slack Bus Power %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:nline
 p = lp(k);
 q = lq(k);

 Ipq = E(p)*ycp(k) + (E(p) - E(q))*yline(k);
 Spq(k) = E(p)*conj(Ipq); %% Line Flow from P to Q

 Iqp = E(q)*ycq(k) + (E(q) - E(p))*yline(k);
 Sqp(k) = E(q)*conj(Iqp); %% Line Flow from Q to P

 Sloss(k) = Spq(k) + Sqp(k); %%% Line Loss
end
%%%%%%%%%%%%%% Calculation of System Loss %%%%%%%%
System_Loss = 0;
for k = 1:nline
 System_Loss = System_Loss + Sloss(k);
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
fprintf(op,' Bus Voltages \n\n');
for i = 1:n
 fprintf(op,' %d %f %f\n',i,abs(E(i)),angle(E(i)));
end
fprintf(op,'\n\n');
%%% Print the Line flows
fprintf(op,' Line Flows \n\n');
for i = 1:nline
 fprintf(op,'Spq(%d) = \t',i);
 if imag(Spq(1)<0)
 fprintf(op,'\t%f-i(%f)\t',real(Spq(i)),abs(imag(Spq(i))));
 else
 fprintf(op,'\t%f+i(%f)\t',real(Spq(i)),abs(imag(Spq(i))));
 end
 fprintf(op,'\n');
end
fprintf(op,'\n');
for i = 1:nline
 fprintf(op,'Sqp(%d) = \t',i);
 if imag(Sqp(1)<0)
 fprintf(op,'\t%f-i(%f)\t',real(Sqp(i)),abs(imag(Sqp(i))));
 else
 fprintf(op,'\t%f+i(%f)\t',real(Sqp(i)),abs(imag(Sqp(i))));
 end
 fprintf(op,'\n');
end
fprintf(op,'\n');
%%% Print the Line Loss
fprintf(op,' Line Losses \n\n');
for i = 1:nline
 fprintf(op,'Sloss(%d) = \t',i);
 if imag(Sloss(1)<0)
 fprintf(op,'\t%f-i(%f)\t',real(Sloss(i)),abs(imag(Sloss(i))));
 else
 fprintf(op,'\t%f+i(%f)\t',real(Sloss(i)),abs(imag(Sloss(i))));
 end
 fprintf(op,'\n');
end
fprintf(op,'\n\n');
%%%%%%%%%%%%%% Calculation of System Loss %%%%%%%%
fprintf(op,'Total System Loss = %f + %f\t',real(System_Loss),imag(System_Loss));
fprintf(op,'\n\n');
%%%%%%%% Slack Bus Power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(op,'Slack Bus Power = %f + %f\t',pgen(nslack),qgen(nslack));
fprintf(op,'\n\n');
toc