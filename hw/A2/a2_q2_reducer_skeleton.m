%% A2 Q2 SKELETON -- speed reducer (gear train), Rao 7.22.4 p. 472
%
% Engineering Optimisation (M-EOP), Assignment 2, Task 2.
%
% Seven design variables (Rao Fig. 7.24):
%   x1 = b,  face width            x5 = l2, shaft-2 length between bearings
%   x2 = m,  tooth module          x6 = d1, shaft-1 diameter
%   x3 = z,  teeth on the pinion   x7 = d2, shaft-2 diameter
%   x4 = l1, shaft-1 length
%
% The objective and the bounds are given.  The CONSTRAINT FUNCTION is yours:
% write the eleven constraints from the problem statement in the form
%
%       g(x) = value/limit - 1 <= 0
%
% Requires: Optimization Toolbox (fmincon).

clear; clc

%% ---- objective: total weight (Rao p. 473) -- given --------------------
f = @(x) 0.7854*x(1)*x(2)^2*(3.3333*x(3)^2 + 14.9334*x(3) - 43.0934) ...
       - 1.508*x(1)*(x(6)^2 + x(7)^2) + 7.477*(x(6)^3 + x(7)^3) ...
       + 0.7854*(x(4)*x(6)^2 + x(5)*x(7)^2);

%% ---- bounds (the ranges in the problem statement) -- given --------------------------------
lb = [2.6 0.7 17 7.3 7.3 2.9 5.0]';
ub = [3.6 0.8 28 8.3 8.3 3.9 5.5]';

%% =======================================================================
%%  the constraints, EVERY ONE as value/limit - 1 <= 0
%% =======================================================================
% Write a function returning the column vector c of the eleven inequality
% constraint values.  Task 1 requires all of them in the non-dimensional form
% above -- not because it is tidier, but because a multiplier is a price PER
% UNIT OF ITS OWN CONSTRAINT.  Mix stresses, lengths and ratios and the
% multipliers come out in different currencies and cannot be ranked, so part 2
% would be unanswerable.  See L04's l04_demo_heatx.m for a worked example.
%
% SIX are already compared against 1, so you only subtract the 1.  For example
%
%     g1: 27 / (x1 x2^2 x3) <= 1     becomes     27/(x1*x2^2*x3) - 1 <= 0
%
% FIVE are compared against something else, so DIVIDE BY THE LIMIT first:
%
%     g5   stress in shaft 1  <= 1100   ->   ( ... )/1100 - 1 <= 0
%     g6   stress in shaft 2  <=  850   ->   ( ... )/850  - 1 <= 0
%     g7   x2*x3              <=   40   ->   x2*x3/40     - 1 <= 0
%     g9   x1/x2              <=   12   ->   (x1/x2)/12   - 1 <= 0
%     g8   x1/x2              >=    5   ->   5/(x1/x2)    - 1 <= 0
%
% Note g8 carefully: it is a LOWER bound, so the limit goes on TOP -- that is
% what keeps it dimensionless and of order one, and it is a different move from
% simply writing 5 - x1/x2 <= 0.
%
% The eleven, in order:
%   g1  bending stress of gear teeth        g2  surface (contact) stress
%   g3  transverse deflection of shaft 1    g4  transverse deflection of shaft 2
%   g5  stress in shaft 1                   g6  stress in shaft 2
%   g7  x2*x3 <= 40 (space envelope)        g8  x1/x2 >= 5  (cuttable teeth)
%   g9  x1/x2 <= 12 (cuttable teeth)        g10 (1.5x6 + 1.9)/x4 <= 1
%   g11 (1.1x7 + 1.9)/x5 <= 1
%
% Written as a local function at the bottom of this file.

%% ---- solve -- given ---------------------------------------------------
x0   = (lb + ub)/2;
opts = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
                    'MaxFunctionEvaluations',2e5,'MaxIterations',5e3, ...
                    'OptimalityTolerance',1e-10,'ConstraintTolerance',1e-10);
[x, fv, flag, out, lam] = fmincon(f, x0, [], [], [], [], lb, ub, @reducercon, opts);

names = {'g1  tooth bending','g2  surface stress','g3  deflection shaft 1', ...
         'g4  deflection shaft 2','g5  stress shaft 1','g6  stress shaft 2', ...
         'g7  x2*x3 <= 40','g8  x1/x2 >= 5','g9  x1/x2 <= 12', ...
         'g10 hub/bearings shaft 1','g11 hub/bearings shaft 2'};
vars  = {'x1 b  face width','x2 m  module','x3 z  teeth', ...
         'x4 l1 span shaft 1','x5 l2 span shaft 2','x6 d1 dia shaft 1','x7 d2 dia shaft 2'};

fprintf('PART 1 -- the optimum\n\n');
fprintf('  weight f* = %.4f\n', fv);
fprintf('  exit flag %d, %d iterations, %d function evaluations\n\n', ...
        flag, out.iterations, out.funcCount);
fprintf('  %-22s %10s %8s %8s   %s\n', 'design variable', 'x*', 'lb', 'ub', 'at a bound?');
for i = 1:7
    atlb = abs(x(i)-lb(i)) < 1e-6;  atub = abs(x(i)-ub(i)) < 1e-6;
    tag = '';  if atlb, tag = 'AT LOWER'; elseif atub, tag = 'AT UPPER'; end
    fprintf('  %-22s %10.6f %8.2f %8.2f   %s\n', vars{i}, x(i), lb(i), ub(i), tag);
end

c = reducercon(x);
fprintf('\n  %-26s %12s   %s\n', 'constraint', 'g(x*)', 'status');
for j = 1:11
    if c(j) > -1e-6, st = 'ACTIVE'; else, st = 'inactive'; end
    fprintf('  %-26s %12.6f   %s\n', names{j}, c(j), st);
end

% Every g is dimensionless, so every lambda is in the same currency (weight
% units per unit of constraint) and the list can be ranked.
fprintf('\nPART 2 -- multipliers of the ACTIVE constraints, largest first\n\n');
lamg = lam.ineqnonlin;
act  = find(c > -1e-6);
[~, ord] = sort(lamg(act), 'descend');
for j = act(ord)'
    fprintf('  %-26s lambda = %10.2f\n', names{j}, lamg(j));
end
fprintf('\n  active bounds (a bound is a constraint too):\n');
for i = 1:7
    if lam.lower(i) > 1e-6
        fprintf('  %-26s lambda = %10.2f  (at its lower bound)\n', vars{i}, lam.lower(i));
    elseif lam.upper(i) > 1e-6
        fprintf('  %-26s lambda = %10.2f  (at its upper bound)\n', vars{i}, lam.upper(i));
    end
end

xr = [3.5 0.7 17.0 7.3 7.3 3.35 5.29]';
cr = reducercon(xr);
fprintf('\nPART 3 -- Rao''s printed design\n\n');
fprintf('  Rao states f* = 2985.22 at xr = [3.5 0.7 17 7.3 7.3 3.35 5.29]\n');
fprintf('  (a) this objective evaluated at xr gives f(xr) = %.4f\n', f(xr));
fprintf('      difference from his printed value: %.4f\n\n', f(xr) - 2985.22);
fprintf('  (b) his design against the eleven constraints:\n');
for j = 1:11
    if cr(j) > 1e-9, st = '<-- VIOLATED'; else, st = ''; end
    fprintf('  %-26s %12.6f   %s\n', names{j}, cr(j), st);
end
fprintf('\n  our own answer: f* = %.4f, feasible to %.1e\n', fv, max(c));
fprintf('  his: f  = %.4f, but infeasible by %.4f on the worst constraint\n', f(xr), max(cr));

function [c, ceq] = reducercon(x)
ceq = [];
c = [ 27/(x(1)*x(2)^2*x(3)) - 1                                           % g1
      397.5/(x(1)*x(2)^2*x(3)^2) - 1                                      % g2
      1.93*x(4)^3/(x(2)*x(3)*x(6)^4) - 1                                  % g3
      1.93*x(5)^3/(x(2)*x(3)*x(7)^4) - 1                                  % g4
      sqrt((745*x(4)/(x(2)*x(3)))^2 + 16.9e6)/(0.1*x(6)^3*1100) - 1       % g5
      sqrt((745*x(5)/(x(2)*x(3)))^2 + 157.5e6)/(0.1*x(7)^3*850) - 1       % g6
      x(2)*x(3)/40 - 1                                                    % g7
      5/(x(1)/x(2)) - 1                                                   % g8  (a LOWER bound: limit on top)
      (x(1)/x(2))/12 - 1                                                  % g9
      (1.5*x(6) + 1.9)/x(4) - 1                                           % g10
      (1.1*x(7) + 1.9)/x(5) - 1 ];                                        % g11
end
