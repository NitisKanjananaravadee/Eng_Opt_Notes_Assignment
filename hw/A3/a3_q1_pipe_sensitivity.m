%% Clear everything from memory and start from scratch
clear all; close all; clc;

%% ======================================================================= %%
%  A3, TASK 1 -- SENSITIVITY ANALYSIS OF AN EIGHT-JUNCTION PIPE NETWORK
%
%  Engineering Optimisation (M-EOP), Assignment 3, Task 1.
%
%  This is pipe_sensitivity_demo.m (Lecture 5) on the twelve-pipe network of
%  the assignment, with ONE addition: the assignment asks for two settings,
%  and they are chosen by opposite methods, so both are implemented here and
%  each is cross-checked against the other.
%
%      state equation   A(x) p = q ,   A(x) = sum_e x_e * A_e        [m^3/s]
%      Setting A        f = p_1            1 response , 12 variables -> ADJOINT
%      Setting B        p_k , k = 1..7     7 responses,  2 variables -> DIRECT
%
%  Units: x_e conductance [m^3/(s.Pa)] (consistent, so p is in Pa and f in Pa;
%  the numbers are the assignment's non-dimensional ones, q = 1 at junction 1).
%
%  ASSUMPTION, stated as one: the assignment's table gives L_e to 3 decimals
%  (1.118 = sqrt(5)/2 etc.).  Those printed values are used verbatim for the
%  volume budget, NOT recomputed from a geometry, because the budget V = 22 is
%  quoted against them.
%
%  Requires: nothing beyond base MATLAB.
%
%  RUN AND VERIFIED in MATLAB R2026a on 2026-11-14.
%% ======================================================================= %%

%% ---- 1. the network ----------------------------------------------------
% Row e of edges: the two junctions pipe e joins.  Figure 1 of the assignment.
edges = [1 2      % pipe  1
         1 3      % pipe  2
         2 3      % pipe  3
         2 4      % pipe  4
         3 5      % pipe  5
         4 5      % pipe  6
         4 6      % pipe  7
         5 7      % pipe  8
         6 7      % pipe  9
         6 8      % pipe 10
         7 8      % pipe 11
         3 4];    % pipe 12

x = [3.0; 2.5; 1.0; 2.0; 2.0; 1.2; 1.8; 1.6; 0.9; 1.4; 1.1; 0.8];  % conductances
L = [1.118; 1.118; 1.000; 1.000; 1.000; 1.000; 1.000; 1.000; ...   % lengths [m]
     1.000; 1.118; 1.118; 1.414];

nNode   = 8;
outNode = 8;                           % the outlet, held at p = 0
V       = 22;                          % pipe volume budget
pmin    = 0.35;                        % Setting B's contractual minimum

n    = numel(x);                       % 12 pipes
free = setdiff(1:nNode, outNode);      % junctions with unknown pressure, 1..7
m    = numel(free);                    % 7 free junctions
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
nSolveAnalysis = 1;

%% ---- 3. SETTING A -- 12 variables, 1 response ---------------------------
% f = p_1, the pump head.  c selects junction 1; we MINIMISE p_1, so c is +e_1
% (Lecture 5's demo carried a minus only because it wanted p_2 made large).
c = zeros(m,1);  c(1) = 1;
f = c' * p;

% ---- 3a. the adjoint -- one solve per RESPONSE, and there is one ---------
% L = c'p + lam'(Ap - q)  =>  A' lam = -c , then df/dx_e = lam' Ae p.
lam = A' \ (-c);                       % <-- one solve, OUTSIDE any loop
nSolveAdjoint = 1;
gAdjoint = zeros(n,1);
for e = 1:n
    gAdjoint(e) = lam' * Ae{e} * p;    % no solve here: a product
end

% ---- 3b. direct differentiation -- one solve per DESIGN VARIABLE --------
% Computed ONLY to cross-check 3a: 12 solves against the adjoint's 1.
gDirect = zeros(n,1);
for e = 1:n
    dpdx       = A \ (-Ae{e}*p);       % <-- a solve, inside the loop
    gDirect(e) = c' * dpdx;
end
nSolveDirectA = n;

% ---- 3c. this objective is SELF-ADJOINT, and that is checked ------------
% f = p_1 = q'p, so c = q, and A is symmetric: A lam = -c = -q gives
% lam = -p exactly.  The adjoint solve is the analysis again, so it is free,
% and the gradient collapses to a pressure drop squared:
%     df/dx_e = -p' Ae p = -(p_i - p_j)^2  <= 0 .
selfAdjointResidual = norm(lam + p)/norm(p);
pFull = zeros(nNode,1);  pFull(free) = p;          % p_8 = 0
dropSquared = -(pFull(edges(:,1)) - pFull(edges(:,2))).^2;

%% ---- 4. SETTING B -- 2 variables, 7 responses ---------------------------
% Responses: every junction pressure p_k, k = 1..7, each needing the analysis.
designB = [1 2];                       % only pipes 1 and 2 may be re-laid
nB      = numel(designB);

% ---- 4a. direct differentiation -- the whole Jacobian in nB solves ------
jacobianDirect = zeros(m, nB);
for k = 1:nB
    e = designB(k);
    jacobianDirect(:,k) = A \ (-Ae{e}*p);   % column k = dp/dx_e, ALL SEVEN rows
end
nSolveDirectB = nB;

% ---- 4b. the adjoint -- one solve per response, so seven --------------
% Computed ONLY to cross-check 4a: 7 solves against the direct method's 2.
jacobianAdjoint = zeros(m, nB);
for k = 1:m                            % one adjoint solve PER RESPONSE p_k
    ck      = zeros(m,1);  ck(k) = 1;
    lamk    = A' \ (-ck);
    for kk = 1:nB
        jacobianAdjoint(k,kk) = lamk' * Ae{designB(kk)} * p;
    end
end
nSolveAdjointB = m;

%% ---- 5. the constraints that need no analysis --------------------------
% Setting A's budget sum_e L_e x_e <= V is an EXPLICIT function of x: the
% gradient is the length vector itself, known before the network is solved.
volumeUsed     = L' * x;
volumeGradient = L;

% Setting B's objective L_1 x_1 + L_2 x_2 is explicit in the same way.
objectiveGradientB = L(designB);

%% ---- 6. report --------------------------------------------------------
fprintf('================ THE ANALYSIS ================\n');
fprintf('%d junctions, %d pipes, %d free junctions\n', nNode, n, m);
fprintf('pressures p = [%s]\n', sprintf('%.5f ', p));
fprintf('pump head  f = p_1 = %.5f\n', f);
fprintf('volume used  sum L_e x_e = %.4f of the V = %g allowed  (%s)\n', ...
        volumeUsed, V, ternaryString(volumeUsed <= V, 'feasible', 'INFEASIBLE'));

fprintf('\n============ SETTING A -- 12 variables, 1 response ============\n');
fprintf('chosen method: THE ADJOINT (1 response beats 12 variables)\n\n');
fprintf('  pipe   from-to      L_e      x_e        df/dx_e      -(dp_e)^2    direct (check)\n');
for e = 1:n
    fprintf('   %2d     %d-%d      %6.3f   %5.2f   %+12.6e  %+12.6e  %+12.6e\n', ...
            e, edges(e,1), edges(e,2), L(e), x(e), ...
            gAdjoint(e), dropSquared(e), gDirect(e));
end
fprintf('\n  adjoint vs direct      : %.1e   (same algebra, different bookkeeping)\n', ...
        max(abs(gAdjoint - gDirect)));
fprintf('  adjoint vs -(dp_e)^2   : %.1e   (the self-adjoint closed form)\n', ...
        max(abs(gAdjoint - dropSquared)));
fprintf('  ||lam + p||/||p||      : %.1e   (lam = -p: f = q''p is SELF-ADJOINT)\n', ...
        selfAdjointResidual);
fprintf('\n  SOLVES ACTUALLY PERFORMED\n');
fprintf('    analysis                      %d\n', nSolveAnalysis);
fprintf('    adjoint equation              %d   (and even this one is free: lam = -p)\n', nSolveAdjoint);
fprintf('    ------------------------------ \n');
fprintf('    adjoint route, total          %d   (%d if self-adjointness is used)\n', ...
        nSolveAnalysis + nSolveAdjoint, nSolveAnalysis);
fprintf('    direct route would have cost  %d   (1 analysis + %d, one per pipe)\n', ...
        nSolveAnalysis + nSolveDirectA, nSolveDirectA);
fprintf('\n  THE BUDGET NEEDED NO ANALYSIS\n');
fprintf('    d/dx_e (sum L_e x_e) = L_e = [%s]\n', sprintf('%.3f ', volumeGradient));
fprintf('    It is an explicit linear function of x -- p never enters it, so\n');
fprintf('    there is nothing to differentiate through and no solve to perform.\n');

fprintf('\n============ SETTING B -- 2 variables, 7 responses ============\n');
fprintf('chosen method: DIRECT DIFFERENTIATION (2 variables beats 7 responses)\n\n');
fprintf('  the 7 x 2 Jacobian  dp_k/dx_e        and the adjoint cross-check\n');
fprintf('  junction k      dp_k/dx_1      dp_k/dx_2   |   p_k      p_k - pmin\n');
for k = 1:m
    fprintf('      %d        %+12.6e  %+12.6e   | %7.4f   %+8.4f\n', ...
            k, jacobianDirect(k,1), jacobianDirect(k,2), p(k), p(k) - pmin);
end
fprintf('\n  direct vs adjoint      : %.1e\n', ...
        max(max(abs(jacobianDirect - jacobianAdjoint))));
fprintf('\n  SOLVES ACTUALLY PERFORMED\n');
fprintf('    analysis                      %d\n', nSolveAnalysis);
fprintf('    direct differentiation        %d   (one per design variable)\n', nSolveDirectB);
fprintf('    ------------------------------ \n');
fprintf('    direct route, total           %d\n', nSolveAnalysis + nSolveDirectB);
fprintf('    adjoint route would have cost %d   (1 analysis + %d, one per response)\n', ...
        nSolveAnalysis + nSolveAdjointB, nSolveAdjointB);
fprintf('\n  THE OBJECTIVE NEEDED NO ANALYSIS\n');
fprintf('    d/dx (L_1 x_1 + L_2 x_2) = (%.3f, %.3f)\n', objectiveGradientB);
fprintf('    Material used is an explicit function of the design.  It is the\n');
fprintf('    seven PRESSURE constraints that are implicit, not the objective.\n');

fprintf('\n================ READING THE ANSWER ================\n');
[~, best] = min(gAdjoint);             % most negative = biggest reduction
fprintf('  widen PIPE %d (junctions %d-%d): df/dx = %+.6e, the most negative,\n', ...
        best, edges(best,1), edges(best,2), gAdjoint(best));
fprintf('  so one unit of conductance there buys the largest drop in pump head.\n');
fprintf('  Its pressure drop is the network''s largest: |p_%d - p_%d| = %.5f.\n', ...
        edges(best,1), edges(best,2), sqrt(-dropSquared(best)));
fprintf('  every df/dx_e <= 0: max = %+.3e, and df/dx_e = -(p_i - p_j)^2.\n', max(gAdjoint));
fprintf('  A_e is positive semidefinite, so the square can never be negative:\n');
fprintf('  widening ANY pipe lowers the inlet pressure, never raises it.\n');

%% ---- 7. figures -------------------------------------------------------
figure(1); clf;
barh(1:n, gAdjoint, 'FaceColor', [0.15 0.30 0.45], 'EdgeColor', 'w');
set(gca, 'YDir', 'reverse', 'YTick', 1:n);
yticklabels(compose('%d  (%d-%d)', (1:n)', edges(:,1), edges(:,2)));
xlabel('df/dx_e   [pump head per unit conductance]');  ylabel('pipe');
title(sprintf('Setting A: all twelve sensitivities (every one <= 0), f = p_1 = %.4f', f));
grid on;
lightFigure(gcf);
exportgraphics(gcf, 'a3_q1_sensitivities.png', 'Resolution', 300);

figure(2); clf;
imagesc(jacobianDirect);  colorbar;  colormap(flipud(parula));
set(gca, 'XTick', 1:nB, 'XTickLabel', compose('dx_%d', designB'), 'YTick', 1:m);
xlabel('design variable');  ylabel('junction k');
title('Setting B: the 7 x 2 Jacobian  dp_k/dx_e, from two solves');
for k = 1:m
    for kk = 1:nB
        text(kk, k, sprintf('%+.4f', jacobianDirect(k,kk)), ...
             'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end
lightFigure(gcf);
exportgraphics(gcf, 'a3_q1_jacobian.png', 'Resolution', 300);

%% ---- 8. hand-checkable assertions -- this script's test suite ----------
% Each is an identity that must hold, not a value copied from a previous run.
tol = 1e-10;
assert(max(abs(gAdjoint - gDirect)) < tol, ...
       'adjoint and direct disagree in Setting A');
assert(max(abs(gAdjoint - dropSquared)) < tol, ...
       'the self-adjoint closed form -(dp)^2 disagrees');
assert(selfAdjointResidual < tol, 'lam = -p failed: f = q''p should be self-adjoint');
assert(max(max(abs(jacobianDirect - jacobianAdjoint))) < tol, ...
       'direct and adjoint disagree in Setting B');
assert(all(gAdjoint <= tol), 'a sensitivity came out positive: A_e is PSD, so it cannot');
assert(abs(jacobianDirect(1,1) - gAdjoint(1)) < tol, ...
       'row 1 of Setting B''s Jacobian is dp_1/dx_1, which IS Setting A''s df/dx_1');
fprintf('\nall assertions passed\n');

%% ================================================================ %%
%        DO NOT MODIFY BELOW THIS LINE - UNLESS YOU ARE BRAVE        %
%  ================================================================  %

function s = ternaryString(cond, a, b)
    if cond, s = a; else, s = b; end
end

function lightFigure(fig)
% MATLAB R2025a and later follow the DESKTOP THEME, so on a dark desktop every
% exported figure comes out dark-backed with white text -- illegible once it is
% placed on white paper.  Force the light theme just before writing the file.
    try
        theme(fig, 'light');                       % R2025a and later
    catch
        try
            set(fig, 'Color', 'w');
            set(findall(fig, 'Type', 'axes'), 'Color', 'w');
        catch
        end
    end
    drawnow;
end
