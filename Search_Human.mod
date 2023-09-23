%variables
var m s n v theta
    p c lambda r E
    w lambdaE W U k
    H J K Y y
    A I
    ;
    
%exogenous shocks
 varexo epsa;
 
% parameters
parameters ms ss ns vs thetas 
               ps cs lambdas rs Es 
               ws lambdaEs Ws Us ks 
               Hs Js Ks Ys ys 
               As Is 
               beta delta chi Am xi
               deltaE alpha eta rhoa phi
               z kappa;
%Calibrations
@#include "Set_param_value.inc"
rhoa=0.9;

% model
model;
%1
m=Am*(s*(1-n))^xi*v^(1-xi);
%2
theta=v/(1-n);
%3
p=m/(1-n);
%4
1/c=lambda;
%5
beta*lambda(+1)*(1+r(+1)-delta)=lambda;
%6
phi*(1-n)/(E+s)=lambdaE;
%7
beta*(lambda(+1)*w(+1)*n(+1)+lambdaE(+1)*(1-deltaE))=lambdaE;
%8
phi/lambda/(E+s)=beta*lambda(+1)/lambda*p/s*(W(+1)-U(+1));
%9
r=A*alpha*(k/H)^(alpha-1);
%10
W=w*H+beta*lambda(+1)/lambda*((1-chi)*W(+1)+chi*U(+1));
%11
U=z-phi/lambda*log(E+s)
        +beta*lambda(+1)/lambda*(p*W(+1)+(1-p)*U(+1));
%12
J=A*k^alpha*H^(1-alpha)-w*H-r*k
        +beta*lambda(+1)/lambda*((1-chi)*J(+1));
%13
0=-kappa+beta*lambda(+1)/lambda*p/theta*J(+1);
%14
eta*J=(1-eta)*(W-U);
%15
K=n*k;
%16
Y=n*y;
%17
y=A*k^alpha*H^(1-alpha);
%18
Y=c+I+z*(1-n)+kappa*v;
%19
n=(1-chi)*n(-1)+m;
%20
K=(1-delta)*K(-1)+I;
%21
H=(1-deltaE)*H(-1)+E;
%22
log(A)=rhoa*log(A(-1))+epsa;
end;

steady_state_model;
@#include "Steady.inc"
end;
steady;
check;

shocks;
var epsa=0.01^2;
end;
stoch_simul(order=1,irf=20,periods=0);