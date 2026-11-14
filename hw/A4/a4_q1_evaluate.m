%% Clear everything from memory and start from scratch
clear all; close all; clc;

%% ======================================================================= %%
%  A4, TASK 1, STEPS 2-3 -- PREDICT, RESCALE, ANALYSE, AND COMPARE WITH
%                           THE OPTIMISER
%
%  Engineering Optimisation (M-EOP), Assignment 4, Task 1, sections 1.1-1.3.
%
%  THIS IS L07'S l07_6_use_surrogate.m AND l07_3_truss_sizeopt.m IN ONE FILE.
%  Merged deliberately: f comes from the surrogate's design and f* from the
%  optimiser's, and the ratio f/f* only means something if BOTH were scored
%  by the same analysis code on the same ground structure.  Keeping them in
%  one script makes that impossible to get wrong.  The analysis itself,
%  complianceObjective, is copied from the lecture unchanged.
%
%  It answers section 1.1 items 2-5 and section 1.2 item 1 in one run, by
%  looping over the two load angles:
%      +35 deg   inside the +-75 deg the data covers      (section 1.1)
%      +85 deg   OUTSIDE it -- the extrapolation test     (section 1.2)
%  The network is NOT retrained between them, which is what section 1.2 asks.
%
%  Requires: Deep Learning Toolbox, Optimization Toolbox, and the surrogate
%  file written by a4_q1_train_surrogate.m.
%
%  RUN AND VERIFIED in MATLAB R2026a on 2026-11-14.
%% ======================================================================= %%

%% ---- which network, and which load cases -------------------------------
surrogateFile = 'a4_surrogate.mat';       % section 1.1/1.2: one layer of 80
% surrogateFile = 'a4_surrogate_deep.mat';% section 1.3 [Bonus]: [80 80]

loadAngleDeg = [35, 85];                  % +35: interpolation; +85: extrapolation
angleLabel   = {'section 1.1  (+35 deg, inside the data)', ...
                'section 1.2  (+85 deg, OUTSIDE the data)'};

%% ---- the design, fixed for both cases ----------------------------------
LyQuery = 4.0;       % truss depth                                     [m]
VQuery  = 0.20;      % volume budget                                 [m^3]
HQuery  = 0.75;      % load height, as a fraction of Ly

%% ---- the problem parameters, the SAME truss the network was trained on --
E     = 2.1e11;      % Young's modulus                                  [Pa]
P     = 1e5;         % applied load                                     [N]
amin  = 1e-7;        % smallest area a bar may shrink to                [m^2]
amax  = 1e-2;        % largest area a bar may grow to                   [m^2]
nx = 7;   ny = 5;    % nodes across and up
Lx = 6.0;            % span                                             [m]
connectivityType = 'neighbour';           % 106 bars -- what the data used

%% ---- the trained network -----------------------------------------------
load(surrogateFile, 'net', 'tr');
numberOfOutputs = net.outputs{net.numLayers}.size;
fprintf('surrogate %s: 4 inputs -> %d design variables\n', surrogateFile, numberOfOutputs);
fprintf('  hidden layers: %d,  test MSE %.4e,  stopped on ''%s''\n', ...
        net.numLayers-1, tr.best_tperf, strtrim(tr.stop));

%% ---- the ground structure, built once (Ly is the same in both cases) ---
[nodalCoordinates, elementConnectivity] = groundStructure(nx, ny, Lx, LyQuery, connectivityType);
numberOfNodes    = size(nodalCoordinates,1);
numberOfDOFs     = 2*numberOfNodes;
numberOfElements = size(elementConnectivity,1);

% Constants drifting from the ones that made the training data is the one
% failure that would not announce itself, so check rather than hope.
assert(numberOfOutputs == numberOfElements, ...
       'Network predicts %d values; this ground structure has %d bars.', ...
       numberOfOutputs, numberOfElements);

%% Specify which DOFs are fixed -- the whole left edge, as in L07
leftNodes  = find(abs(nodalCoordinates(:,2)) < 1e-12);
fixedDOFs  = sort([2*leftNodes-1; 2*leftNodes])';
freeDOFs   = setdiff(1:numberOfDOFs, fixedDOFs);

%% Compute element lengths (the volume constraint needs them)
elementLength = zeros(numberOfElements,1);
for e = 1:numberOfElements
    nodei = elementConnectivity(e,2);  nodej = elementConnectivity(e,3);
    dx = nodalCoordinates(nodej,2) - nodalCoordinates(nodei,2);
    dy = nodalCoordinates(nodej,3) - nodalCoordinates(nodei,3);
    elementLength(e) = sqrt(dx^2 + dy^2);
end

aPlotRef = VQuery/Lx;      % a FIXED width scale, so every figure compares

%% ======================= the two load cases ============================
results = struct([]);
for k = 1:numel(loadAngleDeg)

loadAngle = loadAngleDeg(k)*pi/180;       % the data stores RADIANS
fprintf('\n\n########## %s ##########\n', angleLabel{k});

%% Define the DOFs with loads
rightNodes = find(abs(nodalCoordinates(:,2) - Lx) < 1e-12);
[~, nearest] = min(abs(nodalCoordinates(rightNodes,3) - HQuery*LyQuery));
loadNode   = rightNodes(nearest);
forceVector               = zeros(numberOfDOFs,1);
forceVector(2*loadNode-1) =  P*sin(loadAngle);          % horizontal component
forceVector(2*loadNode)   = -P*cos(loadAngle);          % vertical component
fprintf('load: node %d at (%.2f, %.2f) m, (fx,fy) = (%+.3e, %+.3e) N\n', ...
        loadNode, nodalCoordinates(loadNode,2), nodalCoordinates(loadNode,3), ...
        forceVector(2*loadNode-1), forceVector(2*loadNode));

%% ---- 1. ASK THE NETWORK ------------------------------------------------
query = [LyQuery; VQuery; HQuery; loadAngle];      % the four inputs, in order
net(query);                                   % Warm-up: the FIRST call times MATLAB
timer = tic;
xPredicted = net(query);                      % The design variable, straight out
secondsSurrogate = toc(timer);

numberOutside = sum(xPredicted < 0 | xPredicted > 1);
xPredicted = min(max(xPredicted, 0), 1);      % It predicts outside [0,1] sometimes
aRaw = amin + (amax - amin)*xPredicted;       % Recover areas
volumeRaw = aRaw'*elementLength;

fprintf('\nsurrogate answered in %.2e s  (%d of %d outputs fell outside [0,1])\n', ...
        secondsSurrogate, numberOutside, numberOfElements);

%% ---- 2. DOES ITS ANSWER OBEY THE CONSTRAINT? ---------------------------
fprintf('\nITEM 2 -- THE RAW PREDICTION\n');
fprintf('   volume %.4f m^3 of the %.4f allowed  -->  %+.1f%% %s\n', ...
        volumeRaw, VQuery, 100*(volumeRaw/VQuery - 1), ...
        ternaryString(volumeRaw > VQuery, 'OVER BUDGET (violates it)', 'under budget (satisfies it)'));
fRaw = complianceObjective(xPredicted, nodalCoordinates, elementConnectivity, ...
                           E, forceVector, freeDOFs, amin, amax, 1);
fprintf('   compliance as predicted           f = %.2f J\n', fRaw);
fprintf('   -- NOT comparable: the network was never told the budget exists.\n');

%% ---- 3. MAKE IT LEGAL, THEN ANALYSE IT ---------------------------------
a = aRaw * (VQuery/volumeRaw);                % Scale every bar until it fits
xScaled = (a - amin)/(amax - amin);
fScaled = complianceObjective(xScaled, nodalCoordinates, elementConnectivity, ...
                        E, forceVector, freeDOFs, amin, amax, 1);  % f0 = 1 -> joules
barsScaled = sum(a >= 0.01*max(a));

fprintf('\nITEM 3 -- AFTER RESCALING TO EXACTLY THE BUDGET\n');
fprintf('   volume %.4f m^3                    (the constraint is now tight)\n', a'*elementLength);
fprintf('   compliance                        f = %.2f J   <-- the comparable one\n', fScaled);
fprintf('   bars carrying real material: %d of %d\n', barsScaled, numberOfElements);
if any(a < amin) || any(a > amax)
    fprintf('   NOTE: rescaling pushed %d areas outside [amin, amax]\n', ...
            sum(a < amin | a > amax));
end

%% ---- 4. THE GROUND TRUTH: run the optimiser on the same four numbers ---
% This block is l07_3_truss_sizeopt.m's set-up, unchanged.
a0 = VQuery/sum(elementLength)*ones(numberOfElements,1);
x0 = (a0 - amin)/(amax - amin);
lb = zeros(numberOfElements,1);
ub = ones(numberOfElements,1);
Aineq = (amax - amin)*elementLength'/VQuery;
bineq = 1 - amin*sum(elementLength)/VQuery;

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
[xStar, fhatStar, exitflag, output] = ...
    fmincon(objective, x0, Aineq, bineq, [], [], lb, ub, [], options);
secondsOptimiser = toc(timer);
aStar = amin + (amax - amin)*xStar;
fStar = fhatStar*f0;
barsStar = sum(aStar >= 0.01*max(aStar));

fprintf('\nITEM 4 -- THE GROUND TRUTH (the optimiser, on the same four numbers)\n');
fprintf('   exit flag %d, %d iterations, %.2f s\n', exitflag, output.iterations, secondsOptimiser);
fprintf('   compliance                       f* = %.2f J\n', fStar);
fprintf('   bars carrying real material: %d of %d\n', barsStar, numberOfElements);
assert(exitflag == 1, 'the ground-truth optimisation did not converge; f* is not usable');

%% ---- 5. COMPARE --------------------------------------------------------
fprintf('\nITEM 5 -- COMPARE\n');
fprintf('   f  (surrogate, rescaled) = %8.2f J\n', fScaled);
fprintf('   f* (optimiser)           = %8.2f J\n', fStar);
fprintf('   f/f*                     = %8.4f   (%+.1f%% worse than optimal)\n', ...
        fScaled/fStar, 100*(fScaled/fStar - 1));
fprintf('   bars: %d (surrogate) against %d (optimiser)\n', barsScaled, barsStar);
fprintf('   speed: %.2e s against %.2f s  --> %.0fx faster\n', ...
        secondsSurrogate, secondsOptimiser, secondsOptimiser/secondsSurrogate);

%% ---- the two structures, side by side ---------------------------------
figure(k); clf;
set(gcf, 'Position', [100 100 1100 380]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plotTrussFixed(nodalCoordinates, elementConnectivity, a, aPlotRef, ...
               leftNodes, loadNode, forceVector, Lx, LyQuery);
title(sprintf('surrogate, rescaled: f = %.1f J\n%d bars, %.0e s', ...
              fScaled, barsScaled, secondsSurrogate), 'FontSize', 9);
nexttile;
plotTrussFixed(nodalCoordinates, elementConnectivity, aStar, aPlotRef, ...
               leftNodes, loadNode, forceVector, Lx, LyQuery);
title(sprintf('optimiser (ground truth): f^* = %.1f J\n%d bars, %.1f s', ...
              fStar, barsStar, secondsOptimiser), 'FontSize', 9);
title(tl, sprintf('L_y = %.1f m, V = %.2f m^3, H = %.2f, theta = %+d deg   -->   f/f^* = %.3f', ...
                  LyQuery, VQuery, HQuery, loadAngleDeg(k), fScaled/fStar));
lightFigure(gcf);
exportgraphics(gcf, sprintf('a4_q1_compare_%ddeg.pdf', loadAngleDeg(k)));

results(k).angle      = loadAngleDeg(k);
results(k).volumeRaw  = volumeRaw;
results(k).overBudget = 100*(volumeRaw/VQuery - 1);
results(k).fRaw       = fRaw;
results(k).fScaled    = fScaled;
results(k).fStar      = fStar;
results(k).ratio      = fScaled/fStar;
results(k).barsScaled = barsScaled;
results(k).barsStar   = barsStar;
results(k).outside    = numberOutside;

end  % ---------------------------------------------- end of the two cases

%% ---- the summary table, interpolation against extrapolation -----------
fprintf('\n\n================ SUMMARY: +35 deg AGAINST +85 deg ================\n');
fprintf('  %-8s %10s %12s %10s %10s %9s %8s %8s\n', 'theta', 'raw vol %', 'f raw [J]', ...
        'f [J]', 'f* [J]', 'f/f*', 'bars f', 'bars f*');
for k = 1:numel(results)
    fprintf('  %+7d %10.1f %12.1f %10.1f %10.1f %9.4f %8d %8d\n', ...
            results(k).angle, results(k).overBudget, results(k).fRaw, ...
            results(k).fScaled, results(k).fStar, results(k).ratio, ...
            results(k).barsScaled, results(k).barsStar);
end
fprintf('\n  f/f* grew from %.3f to %.3f: %.1fx further from optimal outside the box.\n', ...
        results(1).ratio, results(2).ratio, results(2).ratio/results(1).ratio);
fprintf('  Neither prediction announced a problem -- that is the lesson.\n');

%% ================================================================ %%
%  COPIED UNCHANGED from the optimisation -- do not edit them.       %
%  Editing them here would score this design by a different rule     %
%  than the one that produced f*, and f/f* would mean nothing.       %
%  ================================================================  %

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
                        supportNodes, loadNode, forceVector, Lx, Ly)
% plotTruss from the lecture, with ONE change: the axes are fixed to the full
% Lx x Ly domain, so the two panels share a frame and a thickness means the
% same area in both.  The width rule a/aPlotRef is the lecture's, unchanged.
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
         'MarkerSize', 7, 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor', 'none');
    % Draw the load in its true direction, so +85 deg looks nearly horizontal.
    fx = forceVector(2*loadNode-1);  fy = forceVector(2*loadNode);
    scale = 0.35*Ly/hypot(fx, fy);
    quiver(nodalCoordinates(loadNode,2), nodalCoordinates(loadNode,3), ...
           scale*fx, scale*fy, 0, 'LineWidth', 2, 'Color', [0.85 0.20 0.10], ...
           'MaxHeadSize', 3);
    axis equal; axis off;
    xlim([0 Lx] + 0.10*Lx*[-1 1]);           % THE SAME frame in every panel
    ylim([0 Ly] + 0.18*Ly*[-1 1]);
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
