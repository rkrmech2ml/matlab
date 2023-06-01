function Xs = NewtonRoot1
imax=1000;
a_err=10^-3;
r_err=10^-10;
i=0;
Xest=0,X1= 0;
F0= funE(X1);
for i = 1:imax

F1= funE(Xest);
FD1= funD(Xest);
if FD1==0
    fprintf( ' Not possible to execute Newtons Methode.\n');
    break
end    

Xi = Xest - (funE(Xest)/funD(Xest));
    if (abs(Xi-Xest)) <((r_err*abs(F0))+a_err)
    Xs = Xi;
    fprintf('The root has found in %i iterations. \n',i);
    fprintf('The root of the equation is');
    disp (Xs);
    break
    end
Xest = Xi;
A(i)=Xi;
end
disp (A);
if i== imax

fprintf('Solution was not cbtained in %i iterations. \n' ,imax);
fprintf('Please change problem parameters. \n');
end
%function defenition
function y = funE(x)
y = exp(-x)-(10^-9);
% Derivative
function y = funD(x)
y = -exp(-x); 	