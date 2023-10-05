%variables
var m s n v theta
    p c lambda r E
    w lambdaE W U k
    H J K Y y
    A I EplusS
    ;
    
%exogenous shocks
 varexo epsa;
 
% parameters
parameters ms ss ns vs thetas 
               ps cs lambdas rs Es 
               ws lambdaEs Ws Us ks 
               Hs Js Ks Ys yss 
               As Is EplusSs
               bet delta chi Am xi
               delta_e alph eta rhoa phi
               z kappa;
%Calibrations
@#include "Set_param_value.inc"

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
bet*lambda(+1)*(1+r-delta)=lambda;
%6
phi*(1-n)/(E+s)=lambdaE;
%7
bet*(lambda(+1)*w*n+lambdaE(+1)*(1-delta_e))=lambdaE;
%8
phi/lambda/(E+s)=bet*lambda(+1)/lambda*p/s*(W(+1)-U(+1));
%9
r=A*alph*(k/H)^(alph-1);
%10
W=w*H+bet*lambda(+1)/lambda*((1-chi)*W(+1)+chi*U(+1));
%11
U=z-phi/lambda*log(E+s)
        +bet*lambda(+1)/lambda*(p*W(+1)+(1-p)*U(+1));
%12
J=A*k^alph*H^(1-alph)-w*H-r*k
        +bet*lambda(+1)/lambda*((1-chi)*J(+1));
%13
0=-kappa+bet*lambda(+1)/lambda*p/theta*J(+1);
%14
eta*J=(1-eta)*(W-U);
%15
K=n*k;
%16
Y=n*y;
%17
y=A*k^alph*H^(1-alph);
%18
Y=c+I+z*(1-n)+kappa*v;
%19
n=(1-chi)*n(-1)+m;
%20
K=(1-delta)*K(-1)+I;
%21
H=(1-delta_e)*H(-1)+E;
%22
log(A)=rhoa*log(A(-1))+epsa;
%23
EplusS=E+s;
end;
steady;
check;








