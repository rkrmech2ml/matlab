function Xs = simple_NewtonRoot1
imax=2000;
err=10^-3;
i=0;
Xest=0,X1= 0;
Xnext=1;
F0= funE(X1);
disp (F0);
record_point=[];
for i = 1:imax

F1= funE(Xest);
   

Xi =Xest-(((Xest-Xnext)/funE(Xest)-funE(Xnext))*funE(Xnext));
    if (abs(Xi-Xest)) <=err
    Xs = Xi;
    fprintf('The root has found in %i iterations. \n',i);
    fprintf('The root of the equation is');
    disp (Xs);
    break
    end
Xest=Xi;
frecrd = funE(Xest);
record_point=[record_point,Xest];
disp (Xest);
A(i)=Xi;
end
fplot(@ (Xest) exp(-Xest)-10^-9);
hold on
plot(record_point,frecrd,"b*")
%disp (A);
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