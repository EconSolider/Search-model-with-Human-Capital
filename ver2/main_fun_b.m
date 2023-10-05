function F=main_fun_b(x,w,H,beta,chi,z,phi,lambda,E,s,p)
W=x(1);
U=x(2);

F(1)=w*H + beta*((1-chi)*W+chi*U) - W;
F(2)=z-phi/lambda*log(E+s) + beta*(p*W+(1-p)*U) - U;
end