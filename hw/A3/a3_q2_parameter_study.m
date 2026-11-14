%% Clear everything from memory and start from scratch
clear all; close all; clc;

%% ======================================================================= %%
%  A3, TASK 2, SECTION 2.2 -- THE THREE DESIGNS, SIDE BY SIDE
%
%  Engineering Optimisation (M-EOP), Assignment 3, Task 2.
%
%  Runs the baseline and the two parameter studies, then draws all three to
%  ONE scale in ONE frame and prints the comparison table.
%
%      run        what changes      setting
%      baseline   --                (7, 5), 'full'
%      A          connectivity      (7, 5), 'neighbour'
%      B          refinement        (13, 9), 'full'
%
%  WHY THE DRAWING RULE MATTERS.  Bar width is a/aPlotRef with aPlotRef the
%  SAME number in all three panels, and the axes are fixed to the full
%  Lx x Ly domain in all three.  Otherwise a given thickness would mean a
%  different area in each picture, and a design that does not reach the
%  domain edges would be stretched to fill its own panel -- two ways of
%  reading a difference that is not in the areas.
%
%  Requires: Optimization Toolbox.  Run time: a few minutes (run B has 4216
%  bars).
%
%  RUN AND VERIFIED in MATLAB R2026a on 2026-11-14.
%% ======================================================================= %%

%% Define the problem parameters -- identical in all three runs
E     = 2.1e11;      % Young's modulus                                  [Pa]
P     = 1e5;         % external load                                    [N]
amin  = 1e-7;        % smallest area a bar may shrink to                [m^2]
amax  = 1e-2;        % largest area a bar may grow to                   [m^2]
V0    = 0.2;         % volume budget of steel                           [m^3]
Lx    = 6.0;         % domain size                                      [m]
Ly    = 4.0;

%% The three runs.  Everything not named here keeps its baseline value.
runs = struct( ...
    'name',   {'baseline',  'A (neighbour)', 'B (13x9 full)'}, ...
    'nx',     {7,           7,               13}, ...
    'ny',     {5,           5,                9}, ...
    'connectivityType', {'full', 'neighbour', 'full'});

%% Solve all three
results = struct([]);
for k = 1:numel(runs)
    fprintf('==== run %d of %d: %s (%d x %d, %s) ====\n', k, numel(runs), ...
            runs(k).name, runs(k).nx, runs(k).ny, runs(k).connectivityType);
    results(k).out = solveCase(runs(k).nx, runs(k).ny, runs(k).connectivityType, ...
                               Lx, Ly, E, P, amin, amax, V0);
end

%% Solving run B takes minutes, so keep the results: the figures and the
%% table can then be redrawn without optimising again.
save('a3_q2_results.mat', 'results', 'runs', 'Lx', 'Ly', 'V0', 'amin', 'amax');

%% ---- the three structures, same scale, same frame ---------------------
aPlotRef = V0/Lx;                      % the SAME width scale in every panel
figure(1); clf;
set(gcf, 'Position', [100 100 1500 380]);
tl = tiledlayout(1, numel(runs), 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:numel(runs)
    nexttile;
    o = results(k).out;
    plotTrussFixed(o.nodalCoordinates, o.elementConnectivity, o.aStar, aPlotRef, ...
                   o.supportNodes, o.loadNode, Lx, Ly);
    title(sprintf('%s: %d bars\nf^* = %.4e J   (f^*/f_0 = %.4f)', ...
                  runs(k).name, o.numberOfElements, o.fStar, o.fStar/o.f0), ...
          'FontSize', 9);
end
title(tl, sprintf('A3 Task 2.2: three designs, one bar-width scale, one %g x %g m frame', Lx, Ly));
lightFigure(gcf);
exportgraphics(gcf, 'a3_q2_three_designs.png', 'Resolution', 300);

%% ---- the table ---------------------------------------------------------
fprintf('\n================ SECTION 2.2, ITEM 2 -- THE TABLE ================\n');
fprintf('  %-15s %7s %14s %14s %9s %10s %12s %9s\n', ...
        'run', 'bars', 'f0 [J]', 'f* [J]', 'f*/f0', 'bars >1%', 'lambda', 'time [s]');
for k = 1:numel(runs)
    o = results(k).out;
    fprintf('  %-15s %7d %14.6e %14.6e %9.5f %10d %12.6e %9.2f\n', ...
            runs(k).name, o.numberOfElements, o.f0, o.fStar, o.fStar/o.f0, ...
            o.numberKept, o.lambdaVolume, o.elapsed);
end

%% ---- the two comparisons the assignment asks for ----------------------
base = results(1).out;  runA = results(2).out;  runB = results(3).out;

fprintf('\n---- item 3: CONNECTIVITY (baseline ''full'' vs A ''neighbour'') ----\n');
fprintf('  f* baseline = %.6e J,  f* A = %.6e J\n', base.fStar, runA.fStar);
fprintf('  lower f*: %s, by %.1f%%\n', ...
        ternaryString(base.fStar < runA.fStar, 'the BASELINE (full)', 'run A (neighbour)'), ...
        100*abs(runA.fStar - base.fStar)/max(base.fStar, runA.fStar));
fprintf('  It MUST be ''full'': every ''neighbour'' design is a ''full'' design with\n');
fprintf('  the long bars left at amin, so ''neighbour''''s feasible set is a SUBSET.\n');
fprintf('  A subset cannot contain a better minimum -- this is set inclusion, not\n');
fprintf('  an observation about this mesh.\n');
fprintf('  Paid for in bars and time: %d -> %d bars, %.2f -> %.2f s\n', ...
        runA.numberOfElements, base.numberOfElements, runA.elapsed, base.elapsed);

fprintf('\n---- item 4: REFINEMENT (baseline 7x9 cells vs B, both spacings halved) ----\n');
fprintf('  f* baseline = %.6e J,  f* B = %.6e J\n', base.fStar, runB.fStar);
fprintf('  lower f*: %s, by %.1f%%\n', ...
        ternaryString(runB.fStar < base.fStar, 'run B (13x9)', 'the BASELINE (7x5)'), ...
        100*abs(runB.fStar - base.fStar)/max(base.fStar, runB.fStar));
fprintf('  It MUST be B: halving both spacings makes the 7x5 grid a SUBSET of the\n');
fprintf('  13x9 grid, so the coarse feasible set is reproducible on the fine one\n');
fprintf('  (a coarse bar = the collinear fine bars in series, same volume, same\n');
fprintf('  axial stiffness).  The grids NEST, which is what guarantees improvement;\n');
fprintf('  a refinement that does not nest can make f* rise.\n');
fprintf('  Cost: %d -> %d bars (x%.1f), %.2f -> %.2f s (x%.1f)\n', ...
        base.numberOfElements, runB.numberOfElements, ...
        runB.numberOfElements/base.numberOfElements, ...
        base.elapsed, runB.elapsed, runB.elapsed/base.elapsed);

%% ---- assertions: the two inclusions above are predictions, so test them
assert(base.fStar <= runA.fStar*(1 + 1e-6), ...
       '''full'' did not beat ''neighbour'': the set inclusion argument is violated');
assert(runB.fStar <= base.fStar*(1 + 1e-6), ...
       'the nested refinement did not improve f*: check that both spacings halved');
assert(all([base.exitflag runA.exitflag runB.exitflag] == 1), ...
       'a run did not converge (exit flag not 1); its f* must not be compared');
fprintf('\nall assertions passed\n');

%% ================================================================ %%
%  COPIED UNCHANGED from a3_q2_truss_sizeOpt.m -- do not edit them.  %
%  Editing them here would score these designs by a different rule   %
%  than the one that produced the baseline's f*.                     %
%  ================================================================  %

function out = solveCase(nx, ny, connectivityType, Lx, Ly, E, P, amin, amax, V0)
% One complete size optimisation.  Same code as a3_q2_truss_sizeOpt.m.
[nodalCoordinates, elementConnectivity] = groundStructure(nx, ny, Lx, Ly, connectivityType);
numberOfNodes    = size(nodalCoordinates,1);
numberOfDOFs     = 2*numberOfNodes;
numberOfElements = size(elementConnectivity,1);

% supports: two pins, at y = Ly/4 and 3Ly/4, snapped to the nearest node
leftNodes     = find(abs(nodalCoordinates(:,2)) < 1e-12);
supportHeight = [Ly/4; 3*Ly/4];
supportNodes  = zeros(numel(supportHeight),1);
for k = 1:numel(supportHeight)
    [~, nearest]    = min(abs(nodalCoordinates(leftNodes,3) - supportHeight(k)));
    supportNodes(k) = leftNodes(nearest);
end
supportNodes = unique(supportNodes);
assert(numel(supportNodes) == numel(supportHeight), ...
       'ny = %d put both supports on the same node.', ny);
fixedDOFs = sort([2*supportNodes-1; 2*supportNodes])';
freeDOFs  = setdiff(1:numberOfDOFs, fixedDOFs);

% load: P downward at (Lx, Ly/2)
rightNodes = find(abs(nodalCoordinates(:,2) - Lx) < 1e-12);
[~, nearest] = min(abs(nodalCoordinates(rightNodes,3) - Ly/2));
loadNode   = rightNodes(nearest);
forceVector             = zeros(numberOfDOFs,1);
forceVector(2*loadNode) = -P;

elementLength = zeros(numberOfElements,1);
for e = 1:numberOfElements
    nodei = elementConnectivity(e,2);  nodej = elementConnectivity(e,3);
    dx = nodalCoordinates(nodej,2) - nodalCoordinates(nodei,2);
    dy = nodalCoordinates(nodej,3) - nodalCoordinates(nodei,3);
    elementLength(e) = sqrt(dx^2 + dy^2);
end

a0 = V0/sum(elementLength)*ones(numberOfElements,1);
x0 = (a0 - amin)/(amax - amin);
lb = zeros(numberOfElements,1);
ub = ones(numberOfElements,1);
Aineq = (amax - amin)*elementLength'/V0;
bineq = 1 - amin*sum(elementLength)/V0;

[f0, ~] = complianceObjective(x0, nodalCoordinates, elementConnectivity, ...
                              E, forceVector, freeDOFs, amin, amax, 1);
objective = @(x) complianceObjective(x, nodalCoordinates, elementConnectivity, ...
                                     E, forceVector, freeDOFs, amin, amax, f0);
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient', true, ...
    'HessianApproximation', 'lbfgs', ...
    'OptimalityTolerance', 1e-6, ...
    'MaxIterations', 5e3, 'MaxFunctionEvaluations', 5e3, ...
    'Display', 'none');

timer = tic;
[xStar, fhatStar, exitflag, output, lambda] = ...
    fmincon(objective, x0, Aineq, bineq, [], [], lb, ub, [], options);
elapsed = toc(timer);

aStar = amin + (amax - amin)*xStar;
numberVoid = sum(aStar < 0.01*max(aStar));

out = struct('nodalCoordinates', nodalCoordinates, ...
             'elementConnectivity', elementConnectivity, ...
             'numberOfElements', numberOfElements, ...
             'supportNodes', supportNodes, 'loadNode', loadNode, ...
             'aStar', aStar, 'f0', f0, 'fStar', fhatStar*f0, ...
             'exitflag', exitflag, 'iterations', output.iterations, ...
             'firstorderopt', output.firstorderopt, ...
             'lambdaVolume', lambda.ineqlin(1), ...
             'physicalVolume', sum(aStar.*elementLength), ...
             'numberKept', numberOfElements - numberVoid, ...
             'elapsed', elapsed);
fprintf('   %d bars, flag %d, %d iters, f0 = %.4e, f* = %.4e, %.2f s\n', ...
        numberOfElements, exitflag, output.iterations, f0, out.fStar, elapsed);
end

function [fhat, dfhatdx] = complianceObjective(x, nodalCoordinates, ...
                 elementConnectivity, E, forceVector, freeDOFs, amin, amax, f0)
a = amin + (amax - amin)*x;                     % physical areas

% ---- assemble, exactly as in your finite-element course ------------------
numberOfElements = size(elementConnectivity,1);
numberOfDOFs     = 2*size(nodalCoordinates,1);
globalK     = zeros(numberOfDOFs, numberOfDOFs);
Kunit       = cell(numberOfElements,1);         % element matrix at UNIT area
elementDOFs = zeros(numberOfElements,4);
for e = 1:numberOfElements
    nodei = elementConnectivity(e,2);
    nodej = elementConnectivity(e,3);
    dofs  = [2*nodei-1, 2*nodei, 2*nodej-1, 2*nodej];
    xi = nodalCoordinates(nodei,2);  yi = nodalCoordinates(nodei,3);
    xj = nodalCoordinates(nodej,2);  yj = nodalCoordinates(nodej,3);
    Le = sqrt((xj-xi)^2 + (yj-yi)^2);
    l  = (xj-xi)/Le;   m = (yj-yi)/Le;          % direction cosines
    localK = E/Le*[ l^2   l*m  -l^2  -l*m
                    l*m   m^2  -l*m  -m^2
                   -l^2  -l*m   l^2   l*m
                   -l*m  -m^2   l*m   m^2 ];
    elementDOFs(e,:) = dofs;
    Kunit{e}         = localK;                  % = dK_e/da_e, the adjoint needs it
    globalK(dofs,dofs) = globalK(dofs,dofs) + a(e)*localK;
end

% ---- solve, condensed onto the free DOFs ---------------------------------
q = zeros(numberOfDOFs,1);
q(freeDOFs) = globalK(freeDOFs,freeDOFs)\forceVector(freeDOFs);

% ---- compliance ----------------------------------------------------------
f    = forceVector'*q;
fhat = f/f0;                                    % normalised objective

% ---- sensitivity ---------------------------------------------------------
dfda = zeros(numberOfElements,1);
for e = 1:numberOfElements
    qe = q(elementDOFs(e,:));                   % this bar's four displacements
    dfda(e) = -qe' * Kunit{e} * qe;             % self-adjoint: no extra solve
end
dfhatdx = ((amax - amin)/f0)*dfda;              % chain rule d(f/f0)/dx
end

function [nodalCoordinates, elementConnectivity] = groundStructure(nx, ny, Lx, Ly, type)
% Builds the node grid and the candidate bars.
    nodalCoordinates = zeros(nx*ny,3);  n = 0;
    for ix = 1:nx
        for iy = 1:ny
            n = n+1;
            nodalCoordinates(n,:) = [n, (ix-1)*Lx/(nx-1), (iy-1)*Ly/(ny-1)];
        end
    end
    pairs = [];
    for i = 1:nx*ny
        for j = i+1:nx*ny
            dx = round((nodalCoordinates(j,2)-nodalCoordinates(i,2))/(Lx/(nx-1)));
            dy = round((nodalCoordinates(j,3)-nodalCoordinates(i,3))/(Ly/(ny-1)));
            switch type
                case 'neighbour'
                    if max(abs(dx),abs(dy)) == 1, pairs = [pairs; i j]; end
                case 'full'
                    if gcd(abs(dx),abs(dy)) == 1, pairs = [pairs; i j]; end
                otherwise
                    error('connectivityType must be ''neighbour'' or ''full''');
            end
        end
    end
    elementConnectivity = [(1:size(pairs,1))' pairs];
end

function plotTrussFixed(nodalCoordinates, elementConnectivity, a, aPlotRef, ...
                        supportNodes, loadNode, Lx, Ly)
% plotTruss from the lecture, with ONE change: the axes are fixed to the full
% Lx x Ly domain instead of being set from the data, so three panels share a
% frame.  The width rule a/aPlotRef is the lecture's, unchanged.
    hold on;
    aRel = a/max(a);
    for e = 1:size(elementConnectivity,1)
        shade = aRel(e)^0.6;                 % opacity: 1 for the largest bar
        if shade < 0.002, continue; end      % would not register on screen
        ni = elementConnectivity(e,2);  nj = elementConnectivity(e,3);
        xc = [nodalCoordinates(ni,2) nodalCoordinates(nj,2)];
        yc = [nodalCoordinates(ni,3) nodalCoordinates(nj,3)];
        plot(xc, yc, '-', 'Color', [0.15 0.30 0.45 shade], ...
             'LineWidth', 0.05 + 13*a(e)/aPlotRef);
    end
    plot(nodalCoordinates(supportNodes,2), nodalCoordinates(supportNodes,3), '^', ...
         'MarkerSize', 8, 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor', 'none');
    quiver(nodalCoordinates(loadNode,2), nodalCoordinates(loadNode,3), 0, -0.3*Ly, ...
           0, 'LineWidth', 2, 'Color', [0.85 0.20 0.10], 'MaxHeadSize', 3);
    axis equal; axis off;
    xlim([0 Lx] + 0.06*Lx*[-1 1]);           % THE SAME frame in every panel
    ylim([0 Ly] + 0.16*Ly*[-1 1]);
end

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
