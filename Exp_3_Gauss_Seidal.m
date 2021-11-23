clc
clear all

ip = fopen('gsinput.m','r++');
op = fopen('gsoutput.m','w++');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Read   BASIC DATA    %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = fscanf(ip,'%d',1);

nline = fscanf(ip,'%d',1);

nslack = fscanf(ip,'%d',1);

itermax = fscanf(ip,'%d',1);

alpha = fscanf(ip,'%f',1);

epsilon = fscanf(ip,'%f',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Read    LINE DATA    %%%%%%%%%%%%%%%
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
%%%%%%%%  Read   BUS DATA - 1    %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

busdata1 = fscanf(ip,'%f',[5 n])
busdata1 = busdata1';
bno = busdata1(:,1);
pgen = busdata1(:,2);
qgen = busdata1(:,3);
pload = busdata1(:,4);
qload = busdata1(:,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Read   BUS DATA - 2    %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

busdata2 = fscanf(ip,'%f',[5 n]);
busdata2 = busdata2';
bno1 = busdata2(:,1);
itype = busdata2(:,2);
vsp = busdata2(:,3);
qmin = busdata2(:,4);
qmax = busdata2(:,5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  SPARSITY VECTORS FORMATION  %%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%  LINE ADMITTANCE CALCULATION  %%%%%%%%
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
%%%%%%%  YBUS FORMATION Using Sparsity Tech %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% DIAGONAL ELEMENTS OF YBUS  %%%%%%%%%%%%%%%%

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

%%%%%%% OFF DIAGONAL ELEMENTS OF YBUS  %%%%%%%%%%%%%%%%

for j = 1:2*nline
    k =  adjl(j);
    ypq(j) = -yline(k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Calculate Pinj, Qinj, Qmininj and Qmaxinj  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:n
    pinj(i) = pgen(i) - pload(i);
    qinj(i) = qgen(i) - qload(i);
    qmininj(i) = qmin(i) - qload(i);
    qmaxinj(i) = qmax(i)- qload(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  Assume Flat Voltage Start  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(nslack) = complex(vsp(1),0);

for i = 2:n
    E(i) = complex(1,0);
end

for i = 1:n
    Eold(i) = E(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  Iterative Process Start  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iter = 1:itermax
    
    delta_Emax = 0;
    
    for i = 1:n
        if(i ~= nslack)
            if(itype(i) == 2)    %% PV - Bus Treatmet
                e1 = real(E(i));
                f1 = imag(E(i));
                delta1 = atan(f1/e1);
                enew = vsp(i)*cos(delta1);
                fnew = vsp(i)*sin(delta1);
                Enew = complex(enew,fnew);
                
                Ical = ypp(i)*Enew;
                
                Jstart = itagf(i);
                Jstop = itagto(i);
                
                for j = Jstart:Jstop
                    q = adjq(j);
                    Ical = Ical + ypq(j)*E(q);
                end
                
                Scal = Enew*conj(Ical);
                
                Qcal = imag(Scal);
                
                qinj(i) = Qcal;
                
                if (Qcal < qmininj(i))
                    qinj(i) = qmininj(i);
                end
                
                if (Qcal > qmaxinj(i))
                    qinj(i) = qmaxinj(i);
                end
                
                if (qinj(i) == Qcal)
                    E(i) = Enew;
                end
            end
            
           %% PQ - Bus Treatmet
            
            S(i) = complex(pinj(i),-qinj(i));
            sum = S(i)/conj(E(i));
            
            Jstart = itagf(i);
            Jstop = itagto(i);
            
            for j = Jstart:Jstop
                q = adjq(j);
                sum = sum - ypq(j)*E(q);
            end
            
            E(i) = sum/ypp(i);                      %% Cal. of updated Voltages
            
            delE(i) = E(i) - Eold(i);               %% Cal. del E   values
            
            if (abs(delE(i)) > abs(delta_Emax))     %% Cal. of del Emax
                delta_Emax = abs(delE(i));
            end
            
                        
            E(i) = Eold(i) + alpha*delE(i);       %% Applying the Acceleration factor
            Eold(i) = E(i);                       %% Setting it to old value
            
        end
    end
    
    
    %%%%% Checking for Convergence   %%%%%%%%
    figure(1)
    plot(iter,abs(delta_Emax),'o')
    hold on
    if (abs(delta_Emax) <= epsilon)
        fprintf(op,'\n\nProblem  Converged in %d Iterations\n\n',iter);
        break;
    end
    
    if (iter == itermax)
        fprintf(op,'Problem not converged');
        break;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  Calculation of Line Flows  %%%%%%%%%%%%%%
%%%%%%%%%%%%%  Calculation of Line LOSSES  %%%%%%%%%%%%%
%%%%%%%%%%%%%  Calculation of Slack Bus Power %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:nline
    p = lp(k);
    q = lq(k);
    
    Ipq = E(p)*ycp(k) + (E(p) - E(q))*yline(k);
    Spq(k) = E(p)*conj(Ipq);                       %% Line Flow from P to Q
    
    Iqp = E(q)*ycq(k) + (E(q) - E(p))*yline(k);
    Sqp(k) = E(q)*conj(Iqp);                       %% Line Flow from Q to P
    
    Sloss(k) = Spq(k) + Sqp(k);         %%% Line Loss
end


%%%%%%%%%%%%%  Calculation of Slack Bus Power %%%%%%%%%%


Jstart  = itagf(nslack);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% PRINT THE RESULTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Print the Bus voltages

for i = 1:n
    fprintf(op,' %d   %f   %f\n',i,abs(E(i)),angle(E(i)));
end

fprintf(op,'\n\n');

%%% Print the Line flows

for i=1:nline
    fprintf(op,'Spq(%d)=\t',i);
    if imag(Spq(1)<0)
        fprintf(op,'\t%f-i(%f)\t',real(Spq(i)),abs(imag(Spq(i))));
    else
        fprintf(op,'\t%f+i(%f)\t',real(Spq(i)),abs(imag(Spq(i))));
    end
    fprintf(op,'\n');
end

fprintf(op,'\n');

for i=1:nline
    fprintf(op,'Sqp(%d)=\t',i);
    if imag(Sqp(1)<0)
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
    if imag(Sloss(1)<0)
        fprintf(op,'\t%f-i(%f)\t',real(Sloss(i)),abs(imag(Sloss(i))));
    else
        fprintf(op,'\t%f+i(%f)\t',real(Sloss(i)),abs(imag(Sloss(i))));
    end
    fprintf(op,'\n');
end

fprintf(op,'\n\n');


%%%%%%%% Slack Bus Power  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(op,'Pslack = %f\t',pgen(nslack));

fprintf(op,'Qslack = %f\t',qgen(nslack));