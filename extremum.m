close all; clear; clc;
% Initialisation of variables to use
x0 = [1.21;-1.15];
tol = 1e-4;
maxits = 50;
% function
f = @(x,y)((x-1)^4 + 2*(x-1)^2*(y-1)^2-2*(y-1)^2-(2*(y+1))^2+1+(y+1)^4);
% Gradient
g= @(x,y) [[4* (x-1) ^3+4*((y-1) ^2)* (x-1)]; [4*((x-1)^2)*(y-1)-4*(y-1)-4*(2*y+1)+4*(y+1)^3]];
% Hessian
h= @(x,y) [[ 12*((x-1)^2)+4*((y-1)^2), 8*(x-1)*(y-1) ];[8*(x-1)*(y-1),12*((y+1)^2)+4*((x-1)^2)-12];];
% Call to newton's function and displaying our results accordingly
[r, iters, flag] = newton_min(g,h,x0,tol,maxits);
fprintf ("<strong>Newton's method</strong>\n\n");
switch (flag)
    case 0
        fprintf ("There was a convergence on f\n\n");
        fprintf("The minima found is: \n");
        disp(r);
        fprintf("It took %d iterations.\n\n",iters);
    case 1
        fprintf ("There was a convergence on x\n\n");
        fprintf("The minima found is: \n");
        disp(r);
        fprintf("It took %d iterations.\n\n",iters);
    otherwise
        fprintf ("There was no convergence\n\n");
        
end
function [r, iters, flag] = newton_min(dg,ddg,x0,tol,maxits)
    x = x0(1); y = x0(2);
    r = NaN;
    flag = -1;
    
    for iters = 1 : maxits
    
        x_old = [x;y];
        
        x_new = x_old - (ddg(x,y)\dg(x,y));
        
        if norm(dg(x,y)) < tol
            
            flag = 0;
            r = x_new;
            return;
        end
        
        if norm(x_new - x_old) <= (tol + eps*norm(x_new))
            
            flag = 1;
            r = x_new;
            return;
        
        end
        
        x = x_new(1);
        y = x_new(2);
    
    end
end