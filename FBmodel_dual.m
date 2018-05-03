function [J,A0,W0,Q]=FBmodel_dual(X,Y,covar,params)
%Solves |J|^2 + lambda_1 |Y - XJ-CW0|^2 + lambda_2 |X^T -JY^T-A0C^T|^2
[n,~]=size(X);
[~,k]=size(covar);
lam1=params.lam1;
lam2=params.lam2;
C=[ones(n,1) covar];
CCC=C/(C'*C);
CCCC=(C/(C'*C))*C';
% obj = @(J,W0,A0)(norm(J,2)^2 + lam1*norm(Y-X*J-C*W0,2)^2 + lam2*norm(X'-J*Y'-A0*C',2)^2);
%% 
K = X*X';
% M = [zeros(k+1,k+1) C' zeros(k+1,n);C -K/(1+lam2*Y'*Y-lam2*Y'*CCCC*Y) -eye(n); zeros(n,k+1) eye(n) -lam1*eye(n)];
M2 = [-K/(1+lam2*(Y'*Y)-lam2*Y'*CCCC*Y)-eye(n)/lam1 C; C' zeros(k+1,k+1)];
b2=[eye(n) + (lam2*K*(CCCC-eye(n)))/(1+lam2*(Y'*Y)-lam2*Y'*CCCC*Y); zeros(k+1,n)];
iM=inv(M2);
% b = [zeros(k+1,1);Y+(-lam2*K*Y+lam2*K*CCCC*Y)/(1+lam2*Y'*Y-lam2*Y'*CCCC*Y);zeros(n,1)];
v = M2\b2*Y;
L = v(1:n,:);
W0 = v(n+1:end);
J = (lam2*X'*Y - lam2*X'*CCCC*Y-X'*L)/(1+lam2*(Y'*Y)-lam2*Y'*CCCC*Y);
A0 = X'*CCC - J*Y'*CCC;

Q= (lam2*X'-lam2*X'*CCCC-X'*iM(1:n,1:n)*b2(1:n,1:n))/(1+lam2*(Y'*Y)-lam2*Y'*CCCC*Y);
% J2=Q*Y;
end