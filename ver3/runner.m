%% parameters
clc,clear; close all;
bet=0.9; delta=0.1;
chi=0.02; Am=0.4;
xi=0.5; delta_e=0.1;
alph=0.35; eta=0.5;
rhoa=0.9; rhos=0.9;

kappa=0.1;
options = optimset('Display','Final','TolX',1e-10,'TolFun',1e-10);
%% steady state
r=1/bet-1+delta;
A=1;
s=0.02;
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
x=fsolve(fun,x0,options);
z=x(1);
lambda=x(2);
phi=x(3);

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
%% save results
ms=m; sss=s; ns=n; vs=v; thetas=theta;
ps=p; cs=c; rs=r; Es=E; ws=w;
lambdaEs=lambdaE; Ws=W; Us=U; ks=k; Hs=H;
Js=J; Ks=K; Ys=Y; yss=y; As=A;
Is=I; lambdas=lambda; EplusSs=EplusS;

save param.mat ms sss ns vs thetas ...
               ps cs rs Es ws ...
               lambdaEs Ws Us ks Hs ...
              Js Ks Ys yss As ...
               Is lambdas EplusSs;
save param.mat rhoa bet delta chi Am xi ...
              delta_e alph eta phi z kappa -append;
          
%% running dynare for stcochastic simulation
% i.e. IRF(impulse response function) for exogenous shocks 
dynare b_analysisA_stoc.mod nolog
% results are contained in the oo_.structure
Y=oo_.irfs.Y_epsa;
figure('Name','Y_epsa');
plot(Y(6:20));

%% running dynare for determinstic simulation
dynare c_analysisB_deter.mod nolog
% results are contained in the simulated_time_series. structure.
Y=Simulated_time_series.Y.data;
E=Simulated_time_series.E.data;
theta=Simulated_time_series.theta.data;

figure('Name','Y');
plot(Y(6:40));
figure('Name','E');
plot(E(6:40));
figure('Name','\theta');
plot(theta(6:40));

