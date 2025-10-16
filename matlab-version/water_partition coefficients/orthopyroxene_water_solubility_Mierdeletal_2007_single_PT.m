T=  1353.03 + 273.15;    %% K
P= 65.085;      %% kbar

A= 0.01354; %% ppm/bar
H= -4563;   %% J/mol
V= 12.1;    %% cc/mol
R = 8.314; % J/(mol K)
Aal= 0.042; %% ppm/bar^0.5
Hal=-79685; %% J/mol
Val= 11.3;  %% cc/mol

fh2o =2243.505080784617*10000; % Gpa to bar      3.616834067114809e+06
c= A*fh2o*exp(-H/(R*T))*exp(-(V*P*100)/(R*T))
cal=Aal*power(fh2o,0.5)*exp(-Hal/(R*T))*exp(-(Val*P*100)/(R*T))
tc=c+cal
