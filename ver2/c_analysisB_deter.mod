@#include "a_param_model.mod"
%*****************************************
% 分析B　parameterの変化による定常解の変化
%****************************************
% 学習とサーチ時間が永続的に1/3から1/2に上昇させて
% 求人倍率を1.1倍から1.5倍に上昇させる
EplusSs=1/2;
endval;
EplusS=1/2;
end;
steady;
% 他変数の変化を計算,結果はsimulated_time_seriesに保存される。
perfect_foresight_setup(periods=50);
perfect_foresight_solver(stack_solve_algo=0); 
