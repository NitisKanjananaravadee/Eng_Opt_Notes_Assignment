%% A2 Q1 SKELETON -- three carts by MINIMUM POTENTIAL ENERGY
%
% Engineering Optimisation (M-EOP), Assignment 2, Task 1, part 3.
%
% Requires: Optimization Toolbox (fminunc).

clear; clc

%% =======================================================================
%%  the data
%% =======================================================================
% Spring rates, N/m (assignment sheet).
k1 = 5000;  k2 = 1500;  k3 = 2000;  k4 = 1000;
k5 = 2500;  k6 =  500;  k7 = 3000;  k8 = 3500;

P  = [1000; 2000; 3000];          % loads on carts 1, 2, 3   [N]
x0 = [0; 0; 0];                   % unloaded configuration   [m]

%% =======================================================================
%%  the potential energy
%% =======================================================================
% Grounded springs stretch by the displacement of the cart they hold;
% coupling springs stretch by the DIFFERENCE of the two they join.
% k3, k7 and k8 all ground cart 3 (k7 and k8 both run to the right wall).
V = @(x) 0.5*k1*x(1)^2 ...                 % k1: left wall -- cart 1
       + 0.5*k2*x(2)^2 ...                 % k2: left wall -- cart 2
       + 0.5*k3*x(3)^2 ...                 % k3: left wall -- cart 3
       + 0.5*k4*(x(2) - x(1))^2 ...        % k4: cart 1 -- cart 2
       + 0.5*k5*(x(3) - x(1))^2 ...        % k5: cart 1 -- cart 3
       + 0.5*k6*(x(3) - x(2))^2 ...        % k6: cart 2 -- cart 3
       + 0.5*k7*x(3)^2 ...                 % k7: cart 3 -- right wall
       + 0.5*k8*x(3)^2 ...                 % k8: cart 3 -- right wall
       - P'*x;                             % work done by the loads

%% ---- solver settings, shared so the comparison is fair -----------------
% Identical settings for both runs, so the ONLY difference is HessUpdate.
%
%   OptimalityTolerance  ON THE SLOPE.  Tests max|grad V|, the first-order
%                        optimality measure -- Lecture 1's grad V = 0 turned
%                        into a tolerance.  "The ground is flat enough."
%                        Stops with exit flag 1.
%
%   FunctionTolerance    ON THE OBJECTIVE.  Tests how much V changed over the
%                        step (for quasi-Newton, the decrease PREDICTED for
%                        the next one).  "V has stopped moving."
%                        Stops with exit flag 3 or 5.
%
%   StepTolerance        ON THE DESIGN.  Tests the length of the step in x.
%                        "x has stopped moving."   Stops with exit flag 2.
%
%   MaxFunctionEvaluations  NOT a tolerance -- a BUDGET.  How many times V may
%                        be evaluated before the solver gives up.  Stops with
%                        exit flag 0, which means NOT CONVERGED: the answer
%                        must not be trusted.  (MaxIterations is the same idea
%                        counted in iterations; left at its default of 400
%                        here, since steepest descent needs about 41.)
%
% The three tolerances can disagree about when to stop, so the settings below
% are chosen deliberately: FunctionTolerance and StepTolerance are set TIGHT
% so they do not end the run early, leaving the optimality test in charge.
% And MaxFunctionEvaluations is NOT decoration -- its default is
% 100*numberOfVariables = 300, while steepest descent needs about 800 here, so
% at the default the run would stop at exit flag 0 having not converged.

tol = 1e-10;                      % optimality tolerance, both methods

%% =======================================================================
%%  run both methods
%% =======================================================================
  [x_sd, V_sd, flag_sd, out_sd] = fminunc(V, x0, ...        % steepest descent
      optimoptions('fminunc', 'Display', 'none', ...
                   'Algorithm', 'quasi-newton', 'HessUpdate', 'steepdesc', ...
                   'OptimalityTolerance', tol, 'FunctionTolerance', 1e-14, ...
                   'StepTolerance', 1e-16, 'MaxFunctionEvaluations', 1e6));

  [x_qn, V_qn, flag_qn, out_qn] = fminunc(V, x0, ...        % quasi-Newton
      optimoptions('fminunc', 'Display', 'none', ...
                   'Algorithm', 'quasi-newton', 'HessUpdate', 'bfgs', ...
                   'OptimalityTolerance', tol, 'FunctionTolerance', 1e-14, ...
                   'StepTolerance', 1e-16, 'MaxFunctionEvaluations', 1e6));


%% ---- report ------------------------------------------------------------
fprintf('                     iterations   f evals   flag        x1         x2         x3\n');
fprintf('  steepest descent %10d %9d %6d  %9.6f %10.6f %10.6f\n', ...
        out_sd.iterations, out_sd.funcCount, flag_sd, x_sd);
fprintf('  quasi-Newton     %10d %9d %6d  %9.6f %10.6f %10.6f\n', ...
        out_qn.iterations, out_qn.funcCount, flag_qn, x_qn);
fprintf('\n  energies: %.6f and %.6f J\n', V_sd, V_qn);
fprintf('  the two answers agree to %.1e m\n', norm(x_sd - x_qn, inf));
fprintf('\n  max|grad V| at the stop: %.2e (SD) and %.2e (BFGS)\n', ...
        out_sd.firstorderopt, out_qn.firstorderopt);
fprintf('  exit flag 5 = the PREDICTED decrease fell below FunctionTolerance.\n');
fprintf('  With a finite-differenced gradient, |grad V| bottoms out near 1e-5,\n');
fprintf('  so OptimalityTolerance = %.0e is never reached and the function test\n', tol);
fprintf('  is what actually stops both runs -- the same test for both, which is\n');
fprintf('  what makes the comparison fair.\n');

%% =======================================================================
%%  The Hessian of V -- for part 4 only
%% =======================================================================
% Neither method above used this.  It is computed here purely so that part 4
% has something to compare the stiffness matrix against.
%
% Central second differences with h = 1e-4.  Not smaller: a second-difference
% formula divides by h^2, so a tiny h amplifies round-off badly.
h = 1e-4;  e = eye(3);  HessV = zeros(3);
for i = 1:3
    for j = 1:3
        HessV(i,j) = (V(x_qn+h*e(:,i)+h*e(:,j)) - V(x_qn+h*e(:,i)-h*e(:,j)) ...
                    - V(x_qn-h*e(:,i)+h*e(:,j)) + V(x_qn-h*e(:,i)-h*e(:,j)))/(4*h^2);
    end
end
fprintf('\n  Hessian of V at the solution:\n');
fprintf('     %9.1f %9.1f %9.1f\n', HessV');
fprintf('  eigenvalues %.1f  %.1f  %.1f\n', sort(eig(HessV)));

%% ---- now write your answer --------------------------------------------
% Report both rows of the table, and answer part 3(d).
% KEEP the Hessian above -- part 4 asks you to compare it with the stiffness
% matrix that the assembly script builds.
