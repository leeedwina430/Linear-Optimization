%% Examples from book
% One can run each part, or uncomment executing line to play on different
% version of IPM implementations.

load rsolvable.mat
IPM(A,b,c);
% IPM_cho(A,b,c);

%% Example from book: Example 3.5
% -------------------------------------
% Simplex Performance                 |
% -------------------------------------
% State:  Solved                      |
% Optimal Value = -136.000000         |
% Total iterations = 3                |
% Running Time = 0.013521             |
% -------------------------------------

A = [1 2 2 1 0 0; 2 1 2 0 1 0; 2 2 1 0 0 1];
b = [20; 20; 20];
c = [-10 -12 -12 0 0 0];

IPM(A,b,c');
% IPM_cho(A,b,c');

%% Example from book: Example 3.8 (Degeneracy Case)
% -------------------------------------
% Simplex Performance                 |
% -------------------------------------
% State:  Solved                      |
% Optimal Value = 1.750000            |
% Total iterations = 4                |
% Running Time = 0.001518             |
% -------------------------------------

A = [1 2 3 0; -1 2 6 0; 0 4 9 0; 0 0 3 1];
b = [3; 2; 5; 1];
c = [1 1 1 0];

IPM(A,b,c');
% IPM_cho(A,b,c');

%% Example from book: Exercise 3.17
% -------------------------------------
% Simplex Performance                 |
% -------------------------------------
% State:  Solved                      |
% Optimal Value = -3.000000           |
% Total iterations = 5                |
% Running Time = 0.002846             |
% -------------------------------------

A = [1 3 0 4 1; 1 2 0 -1 1; -1 -4 3 0 0];
b = [2; 2; 1];
c = [2 3 3 1 -2];

IPM(A,b,c');
% IPM_cho(A,b,c');

%% Given Examples on Elearning
% Try for given nondegenerate cases
% The final outcome table can be seen in the shell.
%
% "scrs8.mat" 
% CVX(SDPT3): 37 iterations; Obj: 904.297
% IPM: 20 iterations; Obj: 904.306
%
% "e226.mat" 
% CVX(SeDuMi): 22 iterations; Obj: -18.7519
% IPM: 19 iterations; Obj: -18.7519
%
% "nug08.mat" 
% CVX(SeDuMi): 8 iterations; Obj: 203.5
% IPM: 9 iterations; Obj: 203.5
% IPM(Modified cholesky): 9 iterations; Obj: 203.5

cvx_solver SeDuMi
filen = {'scrs8','e226','nug08'};

obj = zeros(3,2); time = zeros(3,2); it = zeros(3,2);
for i=1:3
    load(filen{i});
    [~,n] = size(A);
    cvx_begin quiet
    variable x(n)
        minimize(dot(c,x))
        subject to
            A* x == b;
            x >= 0;
    cvx_end
    obj(i,1) = cvx_optval;
    time(i,1) = cvx_cputime;
    it(i,1) = cvx_slvitr;
    [~,~,~,iter,rtime,optv] = IPM(full(A),b,c);
    obj(i,2) = optv;
    time(i,2) = rtime;
    it(i,2) = iter;
    
end

Solver = {'CVX';'IPM'};
scrs8_OPV = [obj(1,1); obj(1,2)];
scrs8_Time = [time(1,1); time(1,2)];
scrs8_Iter = [it(1,1); it(1,2)];

e226_OPV = [obj(2,1); obj(2,2)];
e226_Time = [time(2,1); time(2,2)];
e226_Iter = [it(2,1); it(2,2)];

nug08_OPV = [obj(3,1); obj(3,2)];
nug08_Time = [time(3,1); time(3,2)];
nug08_Iter = [it(3,1); it(3,2)];

table(Solver, scrs8_OPV, scrs8_Iter, scrs8_Time)
table(Solver, e226_OPV, e226_Iter, e226_Time)
table(Solver, nug08_OPV, nug08_Iter, nug08_Time)

%% Comparison with CVX(SDPT3/SeDuMi) on Random Solvable Examples
% The final outcome table can be seen in the shell.
% Compare the average iteration number and running time between my IPM and CVX
% on random generated feasible and bounded examples.

rep = 20;
m = 500; n = 1000;

time = zeros(rep,2); it = zeros(rep,2);
for i=1:rep

    A = rand(m,n)-0.5; 
    x = rand(n,1); s = rand(n,1); y = rand(m, 1);
    c = A'*y + s;
    b = A*x;

    cvx_begin quiet
    variable x(n)
        minimize(dot(c,x))
        subject to
            A* x == b;
            x >= 0;
    cvx_end
    time(i,1) = cvx_cputime;
    it(i,1) = cvx_slvitr;

    [~,~,~,iter,rtime,optv] = IPM(full(A),b,c);
    time(i,2) = rtime;
    it(i,2) = iter;
    
end


Solver = {'CVX';'IPM'};
ave_Time = [mean(time(:,1)); mean(time(:,2))];
ave_Iter = [mean(it(:,1)); mean(it(:,2))];
table(Solver, ave_Iter, ave_Time)

