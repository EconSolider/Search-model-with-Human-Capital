%% parameters
clc,clear;
beta=0.94; delta=0.1;
chi=0.02; Am=0.4;
xi=0.5; deltaE=0.1;
alpha=0.35; eta=0.5;
%% steady state
r=1/beta-1+delta;
A=1;
n=0.5;
theta=1.1;
EplusS=1/3;
v=theta*(1-n);
s=(chi*n/Am*(1-n)^(-xi)*v^(xi-1))^(1/xi);
m=Am*(s*(1-n))^xi*v^(1-xi);
p=m/(1-n);
E=EplusS-s;
H=E/deltaE;
k=(r/alpha/A)^(1/(alpha-1))*H;
y=A*k^alpha*H^(1-alpha);
Y=n*y;
K=n*k;
I=delta*K;
% solve for z,lambda,kappa,phi
x0=[0.05,0.5,0.05,0.4];
fun=@(x) main_fun_a(x,eta,theta,beta,p,y,alpha,E,s,chi,Y,I,n,v,deltaE,H);
x=fsolve(fun,x0);
z=x(1);
lambda=x(2);
kappa=x(3);
phi=x(4);

c=1/lambda;
J=theta/beta/p*kappa;
w=(1-alpha)*y/H+(beta-beta*chi-1)*theta/beta/p*kappa/H;

%solve for W and U
a0=[5,5];
fun2=@(x) main_fun_b(x,w,H,beta,chi,z,phi,lambda,E,s,p);
a=fsolve(fun2,a0);
W=a(1);
U=a(2);
fun_la=@(lambdaE) beta*(lambda*w*n+lambdaE*(1-deltaE))-lambdaE;
lambdaE=fsolve(fun_la,0.5);
%% save results
ms=m; ss=s; ns=n; vs=v; thetas=theta; 
ps=p; cs=c; rs=r; Es=E; ws=w;
lambdaEs=lambdaE; Ws=W; Us=U; ks=k; Hs=H;
Js=J; Ks=K; Ys=Y; ys=y; As=A;
Is=I; lambdas=lambda;

save param.mat ms ss ns vs thetas ... 
               ps cs rs Es ws ...
               lambdaEs Ws Us ks Hs ...
               Js Ks Ys ys As ...
               Is lambdas;
save param.mat beta delta chi Am xi ...
               deltaE alpha eta phi z kappa -append;
 %% running dynare
 dynare Search_Human.mod