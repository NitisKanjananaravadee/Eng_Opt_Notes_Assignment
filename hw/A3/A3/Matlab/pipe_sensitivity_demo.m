%% PIPE_SENSITIVITY_DEMO -- a complete sensitivity analysis, in one page
%
% Engineering Optimisation (M-EOP).  Reference code for Lecture 5 and for
% Assignment 3, Task 1.  It is deliberately MINIMAL and deliberately COMPLETE:
% every line you need is here, on the five-pipe network of the lecture, and
% nothing else is.  Read it, run it, then write your own for the assignment's
% larger network.
%
%   state equation   A(x) p = q ,   A(x) = sum_e x_e * A_e
%   objective        f = c' p                        (c selects and signs)
%   sensitivities    df/dx_e , by BOTH exact methods
%
% Requires: nothing beyond base MATLAB.

clear; clc

%% ---- 1. the network ---------------------------------------------------
edges   = [1 2; 1 3; 2 3; 2 4; 3 4];   % each row: the two junctions a pipe joins
x       = [2; 2; 1; 1; 1];             % conductances -- the design variables
nNode   = 4;
outNode = 4;                           % output junction, held at p = 0

n    = numel(x);
free = setdiff(1:nNode, outNode);      % junctions with unknown pressure
m    = numel(free);
q    = zeros(m,1);  q(1) = 1;          % unit flow injected at junction 1

%% ---- 2. assemble, and solve once --------------------------------------
% A is LINEAR in x, so dA/dx_e is exactly Ae{e}: nothing to differentiate.
Ae = cell(n,1);   A = zeros(m);
for e = 1:n
    T = zeros(nNode);
    i = edges(e,1);  j = edges(e,2);
    T([i j],[i j]) = T([i j],[i j]) + [1 -1; -1 1];
    Ae{e} = T(free,free);
    A     = A + x(e)*Ae{e};
end
p = A \ q;                             % THE ANALYSIS -- one linear solve

%% ---- 3. the objective -------------------------------------------------
% Minimise f = -p_2: junction 2 must be well supplied, and this course
% minimises everything, so a quantity to be made large enters with a minus.
c = zeros(m,1);  c(2) = -1;
f = c' * p;

%% ---- 4a. DIRECT differentiation -- one solve per DESIGN VARIABLE -------
%       A * dp/dx_e = -Ae{e} * p ,   then   df/dx_e = c' * dp/dx_e
g_direct = zeros(n,1);
for e = 1:n
    dpdx        = A \ (-Ae{e}*p);      % <-- a solve, inside the loop
    g_direct(e) = c' * dpdx;
end

%% ---- 4b. ADJOINT -- one solve per RESPONSE ----------------------------
% Lagrangian L = f + lam'(A p - q)  (Lecture 4's sign convention), so
%       A' lam = -c ,   then   df/dx_e = lam' * Ae{e} * p
lam = A' \ (-c);                       % <-- one solve, OUTSIDE any loop
g_adjoint = zeros(n,1);
for e = 1:n
    g_adjoint(e) = lam' * Ae{e} * p;   % no solve here: a product
end

%% ---- 5. report --------------------------------------------------------
fprintf('pressures p        = [%s]\n', sprintf('%.4f ', p));
fprintf('objective  f = -p_2 = %.4f\n\n', f);
fprintf('  pipe    direct      adjoint\n');
for e = 1:n
    fprintf('   %d   %+9.5f   %+9.5f\n', e, g_direct(e), g_adjoint(e));
end
fprintf('\n  the two agree to %.1e -- same algebra, different bookkeeping\n', ...
        max(abs(g_direct - g_adjoint)));
fprintf('  SOLVES:  direct %d (= n_design)   adjoint %d (= n_response)\n', n, 1);
fprintf('  We minimise, so the sign says WHICH WAY TO MOVE, both ways:\n');
fprintf('    negative -> WIDENING that pipe improves the design;\n');
fprintf('    positive -> NARROWING it does.  A positive entry is an instruction,\n');
fprintf('                not a prohibition -- and on this network it is the largest\n');
fprintf('                improvement on offer (see Lecture 5, step 4).\n');
