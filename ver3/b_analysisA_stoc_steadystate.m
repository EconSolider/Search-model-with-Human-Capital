function [ys,params,check] = b_analysisA_stoc_steadystate(ys,exo,M_,options_)
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options   [structure] Dynare options structure
%
% Output: 
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] 0 if steady state computation worked and to
%                        1 of not (allows to impose restrictions on parameters)

%% Step 0: initialize indicator and set options for numerical solver
check = 0;
options = optimset('Display','Final','TolX',1e-10,'TolFun',1e-10);

%% Step 1: read out parameters to access them with their name
for ii = 1:M_.param_nbr
  eval([ M_.param_names{ii} ' = M_.params(' int2str(ii) ');']);
end

%% Step 2: Enter model equations here
r=1/bet-1+delta;
A=As;
s=sss;
n=ns;
EplusS=EplusSs;
m=chi*n;
v=(m/Am*(s*(1-n))^(-xi))^(1/(1-xi));
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

%% Step 3: Update parameters and variables
params=NaN(M_.param_nbr,1);
for iter = 1:M_.param_nbr %update parameters set in the file
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

for ii = 1:M_.orig_endo_nbr %auxiliary variables are set automatically
  eval(['ys(' int2str(ii) ') = ' M_.endo_names{ii} ';']);
end

end