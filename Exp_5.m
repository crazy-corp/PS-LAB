tic;
clc;
clear all;
ip = fopen('ip.m','r++');
op = fopen('op_exp5.m','w++');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Read   BASIC DATA    %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BMVA = fscanf(ip,'%g',1);
BKV = fscanf(ip,'%g',1);
nbus = fscanf(ip,'%d',1);
nline = fscanf(ip,'%d',1);
nslack = fscanf(ip,'%d',1);
epsilon = fscanf(ip,'%f',1);
itermax = fscanf(ip,'%d',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Read    LINE DATA    %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linedata = fscanf(ip,'%f',[6,nline]);
linedata = linedata';
LNo = linedata(:,1);
LP = linedata(:,2);
LQ = linedata(:,3);

RRR = 1;

r = linedata(:,4)*RRR;
x = linedata(:,5);

Ltype = linedata(:,6);

bz = (BKV*BKV)/BMVA;

R = r/bz;
X = x/bz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Read   BUS DATA    %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bus_data = fscanf(ip,'%f',[5,nbus]);
bus_data = bus_data';

bus_data = bus_data./(BMVA*10^3);

BusNo = bus_data(:,1);
PGen = bus_data(:,2);
QGen = bus_data(:,3);

LA = 1;

Pload = bus_data(:,4)*LA;
Qload = bus_data(:,5)*LA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  SPARSITY VECTORS FORMATION  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%       NLCONT(I) Vectror Formation   %%%%%%%%%%%
%    Initilization of NLCONT(I) Vector

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

%%%      FORMATION OF RESERVATION CHART           %%%%%

ITAGF(1) = 1;
ITAGTO(1) = NLCONT(1);

for i = 2:nbus
    ITAGF(i) = ITAGTO(i-1) + 1;
    ITAGTO(i) = ITAGTO(i-1) + NLCONT(i);
end

%%%% FORMATION OF ADJQ(J) & ADJL(J) Vectors   %%%%%
% Initilization of IFILL(I) Vector  %%%%

for i = 1:nbus
    IFILL(i) = 0;
end

for k = 1:nline
    P = LP(k);
    Q = LQ(k);
    LPQ = ITAGF(P) + IFILL(P);
    LQP = ITAGF(Q) + IFILL(Q);
    ADJQ(LPQ) = Q;
    ADJL(LPQ) = k;
    ADJQ(LQP) = P;
    ADJL(LQP) = k;
    IFILL(P) = IFILL(P) + 1;
    IFILL(Q) = IFILL(Q) + 1;
end

% Calculation of Primitive impedance admitance of all lines

for k = 1:nline
    zline(k) = complex(R(k),X(k));
    yline(k) = 1/zline(k);
end


% Assume Flat Voltage Start for all Buses

for i = 1:nbus
    E(i) = complex(1,0);
    Eold(i) = E(i);
end


% ############# For ZIP Model of loads ################

 a0=1;      a1=0;    a2=0;                 % CP Model
%  a0=0;     a1=1;    a2=0;                 % CI Model
% a0=0;     a1=0;    a2=1;                 % CZ Model
%a0=0.4;   a1=0.3;  a2=0.3;               % Composite load ZIP


% Calculation of Pinj and Qinj at each BUS

for i = 2:nbus
    Pinj(i) = Pload(i) - PGen(i);
    Qinj(i) = Qload(i) - QGen(i);
end


% Iterative Process starts

for iter = 0:itermax
    
    delEmax = 0;
    
    %  Calculation of Current Injections
    
    for i = 1:nbus
        SoInj(i) = complex(Pinj(i),Qinj(i));         % Complex power injections
        
        Sinj(i) = SoInj(i)*(a0+a1*E(i)+a2*E(i)^2);
        
        Iinj(i) = conj(Sinj(i)/E(i));
    end
    
    
    %  Calculation of BRANCH CURRENTS in the respective lines
    
    for p = nbus:-1:1
        Ibr(p) = Iinj(p);
        for j = ITAGF(p):ITAGTO(p)
            q = ADJQ(j);
            if(q > p)
                Ibr(p) = Ibr(p) + Ibr(q);
            end
        end
    end
    
    %  Calculation of Voltages at Receving side
    
    for k = 1:nline
        P = LP(k);
        Q = LQ(k);
        E(Q) = E(P) - zline(k)*Ibr(Q);
    end
    
    
    %  Calculation of delE(i) avlues
    
    for i = 2:nbus
        delE(i) = E(i)-Eold(i);
        Eold(i) = E(i);
    end
    
    %  Finding of delEmax value
    
    for i = 2:nbus
        if (abs(delE(i)) > delEmax)
            delEmax = abs(delE(i));
        end
    end
    
    XX(iter+1) = iter;
    YY(iter+1) = delEmax;
    
    %      Check for convergence
    
    if (delEmax <= epsilon)
        fprintf(op, '\nThe Problem is converged in %d iterations\n',iter);
        break;
    end
    
end

% Calculation of Line flows and Line losses

Tot_loss = 0;
for k = 1:nline
    p = LP(k);
    q = LQ(k);
    Ipq = (E(p)-E(q))/zline(k);
    Spq = E(p)*conj(Ipq);
    Iqp = (E(q)-E(p))/zline(k);
    Sqp = E(q)*conj(Iqp);
    Sloss(k) = Spq + Sqp;
    Tot_loss = Tot_loss + Sloss(k);
end

Tot_Ploss = real(Tot_loss)*BMVA;
Tot_Qloss = imag(Tot_loss)*BMVA;


fprintf(op,'\nTotal Active Power Loss, P(MW) = %f\n',Tot_Ploss);
fprintf(op,'Total Reactive Power Loss, Q(MVAr) = %f\n',Tot_Qloss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Voltage Profile of the System %%%%%%%%%%%%

figure
plot(1:nbus,abs(E),'d-k','LineWidth',2)
xlabel('Bus No.');
ylabel('Voltage in P.U');

legend('Bus Voltages');

title('Voltage Profile of the System');


%%%%%%%%%%%%%%%% Convergence PLOT %%%%%%%%%%%%
figure
plot(XX,YY,'p-r','LineWidth',2)
xlabel('No. of Iterations');
ylabel('DelEmax');

legend('DelEmax');

title('Convergence Chars of DLF');



toc;