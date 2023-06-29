%optimising energy in a string under external force applied(mathematical approach) 
%writen Rahul KR
%Assignment for Optimization 1

function step_size = armijo_goldstein
    % Armijo-Goldstein step size rule.
clc 
    % Set default values for optional arguments
        sigma = 0.8;
        beta = 0.8;
    x =[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];%10 variable
    %x =[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
    function fun = fun1(x,n)
         N1 = n;
         a1 = 2;
         b1 = -1;
         c1 = -1;
         B1 = diag(a1*ones(1,N1)) + diag(b1*ones(1,N1-1),1) + diag(c1*ones(1,N1-1),-1);
         A1=B1*(n+1);
         a1 = 4;
         b1 = 1;
         c1 = 1;
         C1 = diag(a1*ones(1,N1)) + diag(b1*ones(1,N1-1),1) + diag(c1*ones(1,N1-1),-1);
         M= C1*(1/(6*(n+1)));
        O=ones(n,1);
        fun=0.5*(x'*A1*x)-O'*M*x;

    end
    function  grad = grad1(x)
         N1 = n;
         a1 = 2;
         b1 = -1;
         c1 = -1;
         B1 = diag(a1*ones(1,N1)) + diag(b1*ones(1,N1-1),1) + diag(c1*ones(1,N1-1),-1);
         A1=B1*(n+1);
         a1 = 4;
         b1 = 1;
         c1 = 1;
         C1 = diag(a1*ones(1,N1)) + diag(b1*ones(1,N1-1),1) + diag(c1*ones(1,N1-1),-1);
         M= C1*(1/(6*(n+1)));
         O=ones(n,1);
         grad= A1*x-M*O;
    end
    function  hess = hess1(n)
        N1 = n;
         a1 = 2;
         b1 = -1;
         c1 = -1;
         B1 = diag(a1*ones(1,N1)) + diag(b1*ones(1,N1-1),1) + diag(c1*ones(1,N1-1),-1);
         A1=B1*(n+1);
         hess =inv(A1);
    end
    m=0;
    step_size = beta^m;
    n=20;
    
    step_sizes = [];
while(m<100 && step_size>0.001)
    direction= grad1(x);
    step_size = beta^m;
    if fun1(x + step_size * -direction,n) < fun1(x,n) + sigma * step_size * grad1(x)' * -direction
        final_sz=step_size;
        break;
    end
    step_sizes = [step_sizes, step_size];
    m=m+1;
end  
fprintf('Step size variation is: %.4f\n',step_sizes);
disp(final_sz)
    
%optimizing the cost function with steepest decent


epsi=10^-4;
v=x;
n=20;%look for x
grad_ini=grad1(v);
fprintf('initial gradient is:\n');
disp(grad_ini);
rhs= epsi*(1+norm(grad_ini));
xnew_arr=[];
funct=[];
m=0;
while(m<100)
%     dk=-grad_ini;
    dk=hess1(n)*grad1(x);
    xnew=x-final_sz*dk;
    if (norm(grad_ini)<=rhs)
        break;
        
    end
    grad_ini =grad1(xnew);
    xnew_arr =[xnew_arr,xnew];
    funct=[funct,fun1(xnew,n)];
    m=m+1;
    x=xnew;
end
fprintf('solution is:\n');
disp(xnew);
fprintf('Dimension of the system is:\n');
disp(size(xnew));
fprintf('Iterations for convergence:');
disp(m);
%disp(xnew_arr)
fprintf('minimum value using SD:')
value=fun1(xnew,n);
disp(value);
fprintf('norm of the gradient at the minimum point:')
disp(norm(grad1(xnew)));

fprintf('Verifying it with  pcg methode:\n')

N1 = 20;
a1 = 2;
b1 = -1;
c1 = -1;
B1 = diag(a1*ones(1,N1)) + diag(b1*ones(1,N1-1),1) + diag(c1*ones(1,N1-1),-1);
A1=B1*(N1+1);
a1 = 4;
b1 = 1;
c1 = 1;
C1 = diag(a1*ones(1,N1)) + diag(b1*ones(1,N1-1),1) + diag(c1*ones(1,N1-1),-1);
M1= C1*(1/(6*(N1+1)));
O=ones(N1,1);
M2=M1*O;

y=pcg(A1,M2);

display(y);

fprintf('minimum of cost function is(using pcg methode:')
fvalue_p=fun1(y,n);
disp(fvalue_p);
fprintf('norm of the gradient at the minimum point:')
disp(norm(grad1(y)));



end

