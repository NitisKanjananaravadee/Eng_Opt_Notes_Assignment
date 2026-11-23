%% Clear everything from memory and start from scratch
clear all; close all; clc;

%% ======================================================================= %%
%  A4, TASK 2 -- RELIABILITY-BASED PLANNING OF THE GROCERY STORE
%
%  Engineering Optimisation (M-EOP), Assignment 4, Task 2.
%
%  Assignment 1's grocery store (Lecture 2's linear programme) with the four
%  farm supplies and the five vegetable requirements made random.  The five
%  steps are Lecture 8's, and the code is l08_2_workshop.m's with new data.
%
%      minimise    sum_f sum_v  c_fv x_fv                   monthly bill [$]
%      subject to  gbar_i(x) + beta_i*sigma_gi(x) <= 0,  i = 1..9
%                  x_fv >= 0
%
%  x_fv = tons of vegetable v bought from farm f -- 20 design variables.
%  y = (S_1..S_4, D_1..D_5), the nine independent normal random variables.
%
%  THE ONE THING THAT MAKES THIS EASIER THAN THE LECTURE'S WORKSHOP: every
%  random variable enters its constraint with coefficient -1 or +1, never
%  multiplied by a design variable.  So sigma_gi does NOT depend on x, the
%  square root is a CONSTANT, and the problem stays LINEAR -- linprog, not
%  fmincon.  See item 1(c).
%
%  Requires: Optimization Toolbox, Statistics and Machine Learning Toolbox.
%
%  RUN AND VERIFIED in MATLAB R2026a on 2026-11-14.
%% ======================================================================= %%

%% ---- the data (Assignment 1's table, now with scatter) ----------------
%                 potato tomato  okra  eggpl spinach
price      = [      200    600   1600    800   1200     % farm 1
                    300    550   1400    850   1100     % farm 2
                    250    650   1500    700   1000     % farm 3
                    150    500   1700    900   1300 ];  % farm 4

supplyMean   = [180; 200; 100; 120];      % S_f mean                   [ton]
supplySd     = [ 18;  20;  10;  12];      % sigma_Sf  (10 % of the mean)
demandMean   = [100;  60;  20;  80;  40]; % D_v mean                   [ton]
demandSd     = [ 10;   6;   2;   8;   4]; % sigma_Dv  (10 % of the mean)

farmName = ["farm 1  "; "farm 2  "; "farm 3  "; "farm 4  "];
vegName  = ["potato  "; "tomato  "; "okra    "; "eggplant"; "spinach "];

nF = numel(supplyMean);                   % 4 farms
nV = numel(demandMean);                   % 5 vegetables
nX = nF*nV;                               % 20 design variables
nG = nF + nV;                             % 9 constraints

constraintName = [compose("Supply, %s", farmName); compose("Requirement, %s", vegName)];

%% ---- ITEM 1 -- the nine constraints, their means and their sigmas -----
% 1(a)  g_f      = sum_v x_fv - S_f          <= 0 ,  f = 1..4
%       g_{4+v}  = D_v - sum_f x_fv          <= 0 ,  v = 1..5
%
% 1(b)  For a FIXED design each g_i is a weighted sum of the nine random
%       variables; the mean takes the weights, the variance the weights
%       squared (Lecture 8).  Here exactly one random variable appears in
%       each g_i, with weight -1 (a supply) or +1 (a requirement):
%
%       gbar_f     = sum_v x_fv - Sbar_f ,     sigma_gf     = sigma_Sf
%       gbar_{4+v} = Dbar_v - sum_f x_fv ,     sigma_g{4+v} = sigma_Dv
%
% 1(c)  sigma_gi does NOT depend on x.  The design variables multiply only
%       the DETERMINISTIC prices and the +-1 incidence coefficients; no
%       random variable is ever multiplied by an x.  In the lecture's
%       workshop a_ij*x_j made sigma a function of x -- here it cannot be.

% The incidence matrix G: row i of G*x gives the x-part of g_i.
G = zeros(nG, nX);
for f = 1:nF
    for v = 1:nV
        G(f, (v-1)*nF + f) = 1;             % sum_v x_fv          (supply row)
        G(nF+v, (v-1)*nF + f) = -1;         % -sum_f x_fv    (requirement row)
    end
end
% x is ordered exactly as in Assignment 1: x(:) of a 4 x 5 matrix, column by
% column, so entry (v-1)*nF + f is "farm f, vegetable v".
gConst = [-supplyMean; demandMean];         % the part of gbar with no x
sigmaG = [supplySd;    demandSd];           % CONSTANT, by item 1(c)

gbar  = @(x) G*x(:) + gConst;               % 9 means
sigma = @(x) sigmaG;                        % 9 standard deviations, fixed

cost = @(x) price(:)'*x(:);

%% ---- ITEM 2 -- the reliability-based problem, in Lecture 8's form -----
% 2(a)  minimise  c'x   s.t.  gbar_i(x) + beta_i*sigma_gi <= 0 ,  x >= 0.
% 2(b)  sigma_gi is a constant, so every constraint is LINEAR in x.  Moved
%       into linprog's A*x <= b form the beta term is just a shifted
%       right-hand side -- the supplies shrink and the requirements grow:
%
%           sum_v x_fv <=  Sbar_f - beta*sigma_Sf        (supply, tightened)
%          -sum_f x_fv <= -(Dbar_v + beta*sigma_Dv)      (requirement, raised)
%
%       So linprog, not fmincon: it is Assignment 1's LP with new numbers.
Aineq = G;
bOf   = @(beta) -(gConst + beta.*sigmaG);   % b = -(gConst + beta*sigma)
lb    = zeros(nX,1);
opts  = optimoptions('linprog', 'Display', 'none');

solveFor = @(beta) linprog(price(:), Aineq, bOf(beta), [], [], lb, [], opts);

%% ---- the Monte Carlo sample: 200 000 months, drawn ONCE ---------------
rng(0);  N = 200000;
Ssample = supplyMean + supplySd.*randn(nF, N);      % 4 x N
Dsample = demandMean + demandSd.*randn(nV, N);      % 5 x N
% g_i <= 0 in a given month, for a given plan:
holdsOf = @(x) [ (G(1:nF,:)*x(:))      <= Ssample           % supply rows
                 (-G(nF+1:end,:)*x(:)) >= Dsample ];        % requirement rows

%% ---- ITEM 3 -- solve at p = 0.95 --------------------------------------
pRequired = 0.95;
betaRequired = norminv(pRequired);                   % 1.6449
xStar = solveFor(betaRequired*ones(nG,1));
plan  = reshape(xStar, nF, nV);

fprintf('================ ITEM 3 -- THE PLAN AT p = %.0f%% ================\n', 100*pRequired);
fprintf('  beta = norminv(%.2f) = %.4f\n\n', pRequired, betaRequired);
fprintf('  Bought (ton)   %8s %8s %8s %8s %8s %10s\n', vegName, "farm tot");
for f = 1:nF
    fprintf('  %-14s %8.2f %8.2f %8.2f %8.2f %8.2f %10.2f\n', ...
            farmName(f), plan(f,:), sum(plan(f,:)));
end
fprintf('  %-14s %8.2f %8.2f %8.2f %8.2f %8.2f %10.2f\n', "veg total", ...
        sum(plan,1), sum(plan(:)));
fprintf('\n  Total cost: $%.2f\n', cost(xStar));

%% ---- ITEM 4 -- check it by Monte Carlo --------------------------------
gAt    = gbar(xStar);
betaAt = -gAt./sigmaG;                               % beta_i(x*)
holds  = holdsOf(xStar);
mcEach = mean(holds, 2);
mcAll  = mean(all(holds, 1));

fprintf('\n================ ITEM 4 -- MONTE CARLO (%d months) ================\n', N);
fprintf('  %-22s %10s %10s %13s %13s\n', 'Constraint', 'gbar_i(x*)', ...
        'beta_i(x*)', 'Phi(beta), %', 'Monte Carlo %');
for i = 1:nG
    fprintf('  %-22s %10.3f %10.4f %12.3f%% %12.3f%%\n', constraintName(i), ...
            gAt(i), betaAt(i), 100*normcdf(betaAt(i)), 100*mcEach(i));
end
fprintf('  %-22s %10s %10s %12.3f%% %12.3f%%\n', 'All nine hold', '--', '--', ...
        100*prod(normcdf(betaAt)), 100*mcAll);
fprintf('\n  4(b): the formula row is the PRODUCT of the nine Phi values above it;\n');
fprintf('  the Monte Carlo row counts the months in which all nine held at once.\n');
fprintf('  The product is legitimate here because each g_i contains exactly ONE\n');
fprintf('  random variable and the nine are independent, so the nine events are\n');
fprintf('  independent too.  Nine constraints at 95%% give 0.95^9 = %.1f%%, not 95%%.\n', ...
        100*0.95^nG);

%% ---- ITEM 5 -- the price of reliability -------------------------------
fprintf('\n================ ITEM 5 -- THE PRICE OF RELIABILITY ================\n');
fprintf('  %10s %12s %14s %18s\n', 'p (%)', 'beta', 'total cost $', 'MC all nine (%)');
pList = [0.50 0.90 0.95 0.99];
costList = zeros(size(pList));  mcList = zeros(size(pList));
for k = 1:numel(pList)
    b  = norminv(pList(k));
    xk = solveFor(b*ones(nG,1));
    costList(k) = cost(xk);
    mcList(k)   = mean(all(holdsOf(xk), 1));
    fprintf('  %10.0f %12.4f %14.2f %17.2f%%\n', 100*pList(k), b, costList(k), 100*mcList(k));
end
fprintf('\n  p = 50%% (beta = 0) is Assignment 1''s deterministic answer, $%.0f.\n', costList(1));
fprintf('  Reliability costs money and the price ACCELERATES: %+.1f%% for 90%%,\n', ...
        100*(costList(2)/costList(1)-1));
fprintf('  %+.1f%% for 95%%, %+.1f%% for 99%% -- because beta grows faster than p near 1.\n', ...
        100*(costList(3)/costList(1)-1), 100*(costList(4)/costList(1)-1));
fprintf('  And the store as a whole always lags the per-constraint figure.\n');

%% ---- ITEM 6 [Bonus] -- plan the whole STORE at 95 %% ------------------
% THE FORMULATION.  Require P[every g_i <= 0] >= 0.95, not each separately.
% Each g_i contains exactly one of the nine independent random variables, so
% the nine events are INDEPENDENT and the joint probability is the product
% EXACTLY -- no bound, no approximation:
%
%     prod_i Phi(beta_i(x)) >= 0.95 ,      beta_i(x) = -gbar_i(x)/sigma_gi
%
% Taking logarithms gives one smooth constraint, written for fmincon as
%
%     c(x) = log(0.95) - sum_i log Phi(beta_i(x)) <= 0 .
%
% beta_i is LINEAR in x and log Phi is concave, so the feasible set is convex
% and the objective linear: the optimum fmincon returns is global.
% The optimiser now SHARES the reliability out instead of being told each
% constraint's share -- that is the whole difference from items 1-5.
pStore = 0.95;

    function [c, ceq] = storeReliability(x, G, gConst, sigmaG, pStore)
        betaI = -(G*x(:) + gConst)./sigmaG;
        logPhi = log(max(normcdf(betaI), realmin));   % realmin: never log(0)
        c   = log(pStore) - sum(logPhi);              % c <= 0
        ceq = [];
    end

nonlcon = @(x) storeReliability(x, G, gConst, sigmaG, pStore);
optsNL  = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...                % as in Lecture 4
    'OptimalityTolerance', 1e-10, ...      % stop when the KKT conditions hold
    'ConstraintTolerance', 1e-10, ...      % ... and the constraint is met
    'MaxFunctionEvaluations', 2e5, ...
    'Display', 'none');
% Start from the per-constraint p = 99 % plan: already close, and feasible
% for the joint requirement is not guaranteed, so fmincon must do real work.
x0Store = solveFor(norminv(0.99)*ones(nG,1));
[xStore, ~, flagStore] = fmincon(@(x) price(:)'*x(:), x0Store, ...
                                 [], [], [], [], lb, [], nonlcon, optsNL);
planStore = reshape(xStore, nF, nV);

gStore    = gbar(xStore);
betaStore = -gStore./sigmaG;
holdsS    = holdsOf(xStore);
mcEachS   = mean(holdsS, 2);
mcAllS    = mean(all(holdsS, 1));

fprintf('\n================ ITEM 6 [Bonus] -- THE WHOLE STORE AT 95%% ================\n');
fprintf('  exit flag %d; constraint used: prod_i Phi(beta_i) >= %.2f\n', flagStore, pStore);
fprintf('\n  Bought (ton)   %8s %8s %8s %8s %8s %10s\n', vegName, "farm tot");
for f = 1:nF
    fprintf('  %-14s %8.2f %8.2f %8.2f %8.2f %8.2f %10.2f\n', ...
            farmName(f), planStore(f,:), sum(planStore(f,:)));
end
fprintf('\n  %-22s %10s %10s %13s %13s\n', 'Constraint', 'gbar_i(x*)', ...
        'beta_i(x*)', 'Phi(beta), %', 'Monte Carlo %');
for i = 1:nG
    fprintf('  %-22s %10.3f %10.4f %12.3f%% %12.3f%%\n', constraintName(i), ...
            gStore(i), betaStore(i), 100*normcdf(betaStore(i)), 100*mcEachS(i));
end
fprintf('  %-22s %10s %10s %12.3f%% %12.3f%%\n', 'All nine hold', '--', '--', ...
        100*prod(normcdf(betaStore)), 100*mcAllS);
fprintf('\n  Total cost: $%.2f\n', cost(xStore));

fprintf('\n  6(c) THE ANALYSIS\n');
fprintf('  cost at item 3 (each constraint 95%%): $%.2f, store works %.2f%% of months\n', ...
        cost(xStar), 100*mcAll);
fprintf('  cost here        (the STORE at 95%%): $%.2f, store works %.2f%% of months\n', ...
        cost(xStore), 100*mcAllS);
fprintf('  the whole store at 95%% costs %+.2f%% more than item 3''s plan.\n', ...
        100*(cost(xStore)/cost(xStar) - 1));
fprintf('  Reliability is no longer shared equally: the optimiser buys it where\n');
fprintf('  it is cheap.  Constraints with slack sit at ~100%% for free; the\n');
fprintf('  binding ones carry the whole 5%% risk between them.\n');

%% ---- assertions: the claims above are predictions, so test them -------
tol = 3e-3;                      % Monte Carlo noise at N = 200 000 is ~0.1 %
assert(all(abs(mcEach - normcdf(betaAt)) < tol), ...
       'a per-constraint Monte Carlo reliability disagrees with Phi(beta)');
assert(abs(mcAll - prod(normcdf(betaAt))) < tol, ...
       'the product rule failed: the nine events are not independent after all');
assert(abs(mcAllS - pStore) < tol, ...
       'the bonus plan does not deliver 95%% at the store level');
assert(flagStore == 1, 'the bonus optimisation did not converge');
assert(abs(costList(1) - 173000) < 1, ...
       'p = 0.50 must reproduce Assignment 1''s $173,000 deterministic optimum');
assert(all(sigma(xStar) == sigma(xStore)), 'sigma moved with x; item 1(c) is wrong');
fprintf('\nall assertions passed\n');

%% ---- figure: what reliability costs -----------------------------------
figure(1); clf;
pFine = [0.50 0.60 0.70 0.80 0.90 0.95 0.975 0.99 0.995 0.999];
costFine = zeros(size(pFine));  mcFine = zeros(size(pFine));
for k = 1:numel(pFine)
    xk = solveFor(norminv(pFine(k))*ones(nG,1));
    costFine(k) = cost(xk);
    mcFine(k)   = mean(all(holdsOf(xk), 1));
end
yyaxis left;
plot(100*pFine, costFine/1000, '-o', 'LineWidth', 1.5);
ylabel('total cost  [$1000]');  xlabel('required reliability of EACH constraint  p  [%]');
yyaxis right;
plot(100*pFine, 100*mcFine, '-s', 'LineWidth', 1.2);
hold on; plot(100*pFine, 100*pFine, ':', 'LineWidth', 1);
ylabel('months in which ALL NINE hold  [%]');
grid on;
legend({'cost', 'store reliability (Monte Carlo)', 'the p you asked for'}, ...
       'Location', 'northwest');
title('The price of reliability -- and why the store always lags each constraint');
lightFigure(gcf);
exportgraphics(gcf, 'a4_q2_price_of_reliability.pdf');

%% ================================================================ %%
%        DO NOT MODIFY BELOW THIS LINE - UNLESS YOU ARE BRAVE        %
%  ================================================================  %

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
