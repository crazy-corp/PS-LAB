clc
%% Point by Point Solution of Swing Equation
%% Problem - f = 50hz generator 50 MVA supplying 48 MW with inertia constant 'H' = 2.7 MJ/MVA at rated speed.
% E = 1.04 pu ,V = 1 pu, X1 = X2 = 0.4 pu. three phase fault on line 2.
% (a) plot swing curve for a susteained fault up to a time of 5 secs.
% (b) plot swing curve if fault is cleared by isolating line in 0.1 seconds.
% (c) Find the critical clearing angle
%% MVA base = 50 given
E = 1.04;
V = 1;
Xd = 0.2;
X1 = 0.4;
X2 = 0.4;
H = 2.7;
P = 48e3;
G = 50e3;
f = 50;
X12 = (X1/2) + Xd ; % Before fault - Reactance between generator and infinite bus
Xfault = 1.0; % During the fault - Reactance between generator and infine bus
X122 = X1 + Xd; % Postfault condition
Pmax = E*V/X12; % Pmax before fault
M = H/(180*f); % angular momentum = H/180*f
Po = P/G; % Power output in pu = 48 MW/50 MVA
delo = asind(Po/Pmax); % initial load angle in degrees //Pe = (E*V/X) sin(delo)
% Prefault condition Power angle Curve
del = 0:pi/10:pi;
Peo = (E*V/X12)*sin(del);
% During fault condition Power angle Curve
Pe1 = (E*V/Xfault)*sin(del);
% Post fault condition Power angle Curve
Pe2 = (E*V/X122)*sin(del); % Power curve after clearing fault
% Figure-1 Power angle Curve
plot(del,Peo);
set(gca,'XTick',0:pi/10:pi);
set(gca,'XTickLabel',{'0','','','','','90','','','','','180'});
title('Power Angle Curves');
xlabel('Load angle');
ylabel('Genpower');
text((2/3)*pi,(1.04/0.4)*sin((2/3)*pi),'\leftarrow Before fault (Pre-fault)','horizontalAlignment','left');
hold all
plot(del,Pe1);
text((2/3)*pi,(1.04/1.0)*sin((2/3)*pi),'\leftarrow During fault','HorizontalAlignment','left');
plot(del,Pe2);
text((2/3)*pi,(1.04/0.6)*sin((2/3)*pi),'\leftarrow Post-fault','HorizontalAlignment','left');
hold off
%% ------------
t = 0.05; % time step preferably 0.05 seconds
t1 = 0:t:0.5;
%% (a) sustained fault at t = 0
% for discontinuity at t = 0 , we take the average of accelerating power
% before and after the fault
% at t = 0-, Pa1 = 0
% at t = 0+, Pa2 = Pi - Pe2
% at t = 0 ,Pa = (Pa1 + Pa2)/2
Pao = (Po - ((E*V/Xfault)*sind(delo)))/2; % at the instant of fault del1 = delo
Pa(1) = Pao;
cdel(1) = 0;
dt = t^2/M;
for i = 1:11

 if i == 1

 d2(i) = dt*Pa(i);
 del(i) = delo;

 else
 cdel(i) = cdel(i-1) + d2(i-1);

 del(i) = del(i-1) + cdel(i);

 Pe(i) = (E*V/Xfault)*sind(del(i));

 Pa(i) = Po - Pe(i);

 d2(i) = dt*Pa(i);
 end

end
%% swing curve - 1 Plot
%figure (2);
%plot(t1,del);
%set(gca,'Xtick',0:0.05:0.5);
%set(gca,'XtickLabel',{'0','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50'});
%title('Swing Curve');
%xlabel('seconds');
%ylabel('degrees');
%text(0.30,110,' Sustained fault','HorizontalAlignment','right');
%text(0.001,130,' load angle increases with time -- Unstable state','HorizontalAlignment','left');
%% (b) Fault cleared in 0.10 seconds ,2nd step ---- 3rd element [1]0 [2]0.05,[3]0.10
Pafo = (Po - ((E*V/Xfault)*sind(delo)))/2; % at the instant of fault del1 = delo
Paf(1) = Pao;
cdelf(1) = 0;
d1f = t^2/M;
Pe2 = (E*V/X122);
for i = 1:2

 if i == 1

 d2f(i) = dt*Pa(i);
 delf(i) = delo;

 else
 cdelf(i) = cdelf(i-1)+d2f(i-1);

 delf(i) = delf(i-1)+cdelf(i);

 Pef(i) = (E*V/Xfault)*sind(delf(i));

 Paf(i) = Po - Pef(i);

 d2f(i) = dt*Paf(i);
 end

end
% after clearing fault, power curve shift to Pe3
for i = 3:11
 if i == 3

 cdelf(i) = cdelf(i-1)+d2f(i-1);
 delf(i) = delf(i-1)+cdelf(i) ;
 Pef(i) = (E*V/Xfault)*sind(delf(i));
 Paf(i) = Po - Pef(i);
 a1 = Paf(i);
 d2f(i) = dt*Paf(i);
 a2 = d2f(i);

 Pef(i) = Pe2*sind(delf(i));
 Paf(i) = Po - Pef(i);
 d2f(i) = dt*Paf(i);

 Paf(i) = (Paf(i)+ a1)/2;
 d2f(i) = (d2f(i) + a2)/2;

 else

 cdelf(i) = cdelf(i-1)+d2f(i-1);

 delf(i) = delf(i-1)+cdelf(i);

 Pef(i) = Pe2*sind(delf(i));

 Paf(i) = Po - Pef(i);

 d2f(i) = dt*Paf(i);
 end
end
%% ------
figure (3);
plot(t1,delf);
set(gca,'Xtick',0:0.05:0.5);
set(gca,'XtickLabel',{'0','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50'});
title('Swing Curve');
xlabel('seconds');
ylabel('degrees');
text(0.25,57,' Fault Cleared in 0.10 sec','HorizontalAlignment','right');
text(0.15,30,' load angle decreases with time -- Stable state','HorizontalAlignment','left');
%% (c) critical clearing angle
Pmax1 = E*V/X12;
Pmax2 = (E*V/Xfault);
Pmax3 = (E*V/X122);
r2 = Pmax3/Pmax1;
r1 = Pmax2/Pmax1;
R = r2-r1;
[x,y] = cart2pol(real(E),imag(E));
Emag = y;
Eangle = x; % angle in radians
delom = 180 - asind(Po/Pmax3);
X = (delom - delo)*pi/180;
Y = sind(delo);
Z = cosd(delom);
deltac = acos(((X*sind(delo))/R)-((r1*cosd(delo))/R)+((r2*cosd(delom)))/R);
fprintf('\nCritical clearing angle is')
deltac*180/pi