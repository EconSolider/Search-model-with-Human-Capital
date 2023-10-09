%% parameters
clc,clear; close all;
bet=0.9; delta=0.1;
chi=0.02; Am=0.4;
xi=0.5; delta_e=0.1;
alph=0.35; eta=0.5;
rhoa=0.9; rhos=0.9;

kappa=0.1;
options = optimset('Display','off','TolX',1e-10,'TolFun',1e-10);
%% steady state
s_list=linspace(0.001,0.067,40);
p_list=zeros(size(s_list));
E_list=zeros(size(s_list));
H_list=zeros(size(s_list));
y_list=zeros(size(s_list));
Welfare_list=zeros(size(s_list));
K_list=zeros(size(s_list));
z_list=zeros(size(s_list));
phi_list=zeros(size(s_list));

for ii=1:length(s_list)
r=1/bet-1+delta;
A=1;
s=s_list(ii);
n=0.6;
EplusS=1/3;
m=chi*n;
v=(m/Am/(s*(1-n))^xi)^(1/(1-xi));
theta=v/(1-n);
p=m/(1-n);
E=EplusS-s;
H=E/delta_e;
k=(r/alph/A)^(1/(alph-1))*H;
y=A*k^alph*H^(1-alph);
Y=n*y;
K=n*k;
I=delta*K;

% solve for z,lambda,kappa,phi
x0=[0.05,0.5,0.05];
fun=@(x) main_fun_a(x,eta,theta,bet,p,y,alph,E,s,chi,Y,I,n,v,delta_e,H,kappa);
[x,Fval,exitflag]=fsolve(fun,x0,options);
z=x(1);
lambda=x(2);
phi=x(3);

if exitflag <= 0
    error('fsolve fails');
end


c=1/lambda;
J=theta/bet/p*kappa;
w=(1-alph)*y/H+(bet-bet*chi-1)*theta/bet/p*kappa/H;

%solve for W and U
a0=[5,5];
fun2=@(x) main_fun_b(x,w,H,bet,chi,z,phi,lambda,E,s,p);
a=fsolve(fun2,a0,options);
W=a(1);
U=a(2);
lambdaE=phi*(1-n)/(E+s);

% contain results to lists
p_list(ii)=p;
E_list(ii)=E;
H_list(ii)=H;
y_list(ii)=y;
Welfare_list(ii)=W-U;
K_list(ii)=K;
z_list(ii)=z;
phi_list(ii)=phi;
end


figure;
plot(s_list, z_list, 'LineWidth', 2);
hold on;
yline(0, 'r-', 'LineWidth', 2);
xlabel('Searching time (s)', 'FontSize', 14);
ylabel('Reservation wage (z)', 'FontSize', 14);
title('Relation between reservation wage and search time', 'FontSize', 16);
ax = gca;
ax.FontSize = 14;


figure;
plot(s_list,phi_list);
figure;
plot(s_list,E_list);
figure;
plot(s_list,H_list);
figure;
plot(s_list,K_list);
figure;
plot(s_list,y_list);
figure;
plot(s_list,Welfare_list);

