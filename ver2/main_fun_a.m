
function F=main_fun_a(x,eta,theta,beta,p,y,alpha,E,s,chi,Y,I,n,v,deltaE,H)
z=x(1);
lambda=x(2);
kappa=x(3);
phi=x(4);

F(1)=eta/(1-eta)*theta/beta/p*kappa ...
    -(1-alpha)*y + z ...
    -phi/lambda*log(E+s) ...
    - (beta-beta*chi-1)*theta/beta/p*kappa ...
    - beta*(1-chi-p)*theta/beta/p*kappa;

F(2)=(1-n)*phi/(E+s) - (1/(1/beta-1+deltaE))*lambda*n ...
    *((1-alpha)*y+(beta-beta*chi-1)*theta/beta/p*kappa)/H;

F(3)=phi/lambda/(E+s) - 1/s*(eta/(1-eta))*theta*kappa;

F(4)= Y - 1/lambda - I - z*(1-n)- kappa*v;

end