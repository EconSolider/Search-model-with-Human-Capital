@#include "a_param_model.mod"
%*****************
% 分析A　IRFs
%*****************
shocks;
var epsa=0.01^2;
end;
stoch_simul(order=1,irf=20,periods=0,nograph);
