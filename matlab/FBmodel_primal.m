function [J,A0,W0,Q]=FBmodel_primal(X,Y,covar,params)
%Solves |J|^2 + lambda_1 |Y - XJ-CW0|^2 + lambda_2 |X^T -JY^T-A0C^T|^2
[n,d]=size(X);
% [~,k]=size(covar);
lam1=params.lam1;
lam2=params.lam2;
C=[ones(n,1) covar];
CCC=C/(C'*C);
CCCC=(C/(C'*C))*C';
%% 
J = (eye(d) + lam1*(X'*X) - lam1*X'*CCCC*X+lam2*(Y'*Y)*eye(d)-lam2*(Y'*CCCC*Y)*eye(d))...
    \((lam1+lam2)*(X'*Y-X'*CCCC*Y));
A0 = X'*CCC - J*Y'*CCC;
W0 = CCC'*Y - CCC'*X*J;
Q = (eye(d) + lam1*(X'*X) - lam1*X'*CCCC*X+lam2*(Y'*Y)*eye(d)-lam2*(Y'*CCCC*Y)*eye(d))...
    \((lam1+lam2)*(X'-X'*CCCC));
end