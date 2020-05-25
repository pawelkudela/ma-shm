function em=environment_degradation(em20,Temp,Moist)


load krzywaET.txt; % Young's modulus temperature dependence experimental data
% two column: temperature, Young's modulus
XY=krzywaET;clear krzywaET;
X=XY(:,1);
Y=XY(:,2);% unscaled
n=6;
[p]=polyfit(X,Y,n);
% F=polyval(p,X);
% x=0:1:100;
% f=polyval(p,x);
% plot(x,f,'-',X,Y,'o');
% T - temperature (0..100 C)
% E - actual Young modulus in temperature T

E20=polyval(p,20);% Youngs modulus at temperature 20C
% degradation factor due to temperature
E=polyval(p,Temp);
dft=E/E20;

load krzywaEM.txt; % Young's modulus moisture dependence experimental data
% two column: moisture, Young's modulus
XY=krzywaEM;clear krzywaEM;
X=XY(:,1);
Y=XY(:,2);% unscaled
n=7;
[p]=polyfit(X,Y,n);
% F=polyval(p,X);
% x=0:0.01:10;
% f=polyval(p,x);
% plot(x,f,'-',X,Y,'o');
% degradation factor due to moisture
E0=polyval(p,0);% Youngs modulus at 0% moisture
E=polyval(p,Moist);
dfm=E/E0;
dfm=1; % curve fitting must be improved, now no effect of moisture
em=em20*dft*dfm;
