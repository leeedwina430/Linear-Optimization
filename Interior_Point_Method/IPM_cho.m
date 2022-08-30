function [x,y,s,iter,rtime,obj] = IPM_cho(A,b,c)
% Solver for standard form LP: min c'*x, Ax=b, x>=0
% Using Predictor-Corrector Interior Point Method
% Using inbuilt function "cho()" for Standard Cholesky Decomposition in MATLAB

t0 = cputime; 
[m,n] = size(A); 

% Initial point
x_ = A'*((A*A')\b);
y_ = ((A*A')\A)*c;
s_ = c-A'*y_;
xhat = x_ + max(-1.5*min(x_), 0) * ones(n,1);
shat = s_ + max(-1.5*min(s_), 0) * ones(n,1);
deltax = 0.5 * (xhat'*shat) / (ones(1,n) * shat);
deltas = 0.5 * (xhat'*shat) / (ones(1,n) * xhat);
x = xhat + deltax * ones(n,1);
y = y_; s = shat + deltas * ones(n,1);

bcsize = max(norm(b), norm(c));
flag = {0};
% Stop Criteria I: Maximum Iteration Number
for iter = 1:100
    Rc = A'*y + s - c;
    Rb = A * x - b; 
    Rxs = x.*s; 
    residual = norm([Rc;Rb;Rxs])/(1+bcsize);

    fprintf('iter %2i: residual = %0.2e', iter, residual);
    fprintf(' obj = %12.6e\n', c'*x);

    % Stop Critetia II: Enough Small Residual
    if residual < 5.e-8
        flag{1} = 1;
        break; 
    end

    % Normal Equation
    d = min(5.e+15, x./s);
    B_ = A * diag(d.^0.5);
    B = B_ * B_'; 
    t1 = x.*Rc - Rxs; 
    t2 = -(Rb + A*(t1./s));     % RHS

    dy = B\t2; 
    dx = (x.*(A'*dy)+t1)./s; 
    ds = -Rc - A'*dy;

    % Step Length
    ap = -1/min(min(dx./x), -1); 
    ad = -1/min(min(ds./s), -1);
    mu = mean(Rxs);
    muaff = [(x + ap*dx)'*(s + ad*ds)]/n;
    sigma = power(muaff/mu,3);
    Rxs = Rxs - sigma*mu + dx.*ds;

    t1 = x.*Rc - Rxs; 
    t2 = -(Rb + A*(t1./s));

    % Search Direction
    % Cholesky Decomposition
    try U = chol(B);
        v = (U')\t2;
        dy = U\v;
    catch
        dy = B\t2;
    end

    dx = (x.*(A'*dy)+t1)./s; 
    ds = -Rc - A'*dy;

    ap = min(0.999*min(-x(dx<0)./dx(dx<0)), 1);
    ad = min(0.999*min(-s(ds<0)./ds(ds<0)), 1);
    
    % Update
    x = x + ap * dx; 
    s = s + ad * ds; 
    y = y + ad * dy;
end

rtime = cputime-t0; obj = c'*x;
fprintf('----------------------------------------\n')
if flag{1} == 1
    fprintf('[IPM] Status: Solved. Obj = %0.2f\n[m, n] = [%g %g], Iterarion = %i, CPU = %g\n\n', obj,m,n,iter,rtime);
else 
    fprintf("Status: Infeasible/Unbounded.\n[m, n] = [%g %g], Iterarion = %i, CPU = %g\n\n", m,n,iter,rtime);
end
end
