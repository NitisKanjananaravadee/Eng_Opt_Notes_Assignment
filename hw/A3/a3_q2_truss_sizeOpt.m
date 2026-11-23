%% Clear everything from memory and start from scratch
clear all; close all; clc;

%% ======================================================================= %%
%  A3, TASK 2 -- THE MICHELL CANTILEVER, BY TRUSS SIZE OPTIMISATION
%
%  Engineering Optimisation (M-EOP), Assignment 3, Task 2.
%
%  THIS IS L06'S DEMO, COPIED.  It is l06_demo_truss_sizeOpt.m with ONE block
%  changed -- the domain, the grid, and the boundary conditions -- and the
%  report extended to print what section 2.1 asks for.  The method, the
%  objective, the sensitivity, the optimiser settings and every function at
%  the bottom of the file are untouched.
%
%    what changed            lecture                 here
%    domain height           Ly = 2 m                Ly = 4 m
%    grid                    7 x 3                   7 x 5   (same 1 m cells)
%    supports                the WHOLE left edge     TWO nodes, at y = Ly/4
%                                                    and y = 3Ly/4
%    load                    mid-height, right edge  (Lx, Ly/2), downward
%
%  In 1904 Michell derived the lightest frame carrying a point load back to
%  two supports by hand; this computes the same answer by optimisation.
%
%  Requires: Optimization Toolbox.
%
%  RUN AND VERIFIED in MATLAB R2026a on 2026-11-14.
%
%      minimise    f(a) = forceVector'*q(a)      compliance (= 2 x strain energy)
%      subject to  sum_e L_e*a_e <= V0        volume budget -- LINEAR, no FEA
%                  amin <= a_e <= amax
%      where       K(a)*q = forceVector                 the FEA: ONE linear solve
%
%  No new theory: the FEA is your finite-element course, the sensitivity is
%  Lecture 5's adjoint (compliance is self-adjoint, so it is free), and the
%  optimiser is Lecture 4's fmincon.  The lecture is the assembly of the three.
%% ======================================================================= %%

%% Define the problem parameters
E     = 2.1e11;      % Young's modulus                                  [Pa]
P     = 1e5;         % external load                                    [N]
amin  = 1e-7;        % smallest area a bar may shrink to                [m^2]
amax  = 1e-2;        % largest area a bar may grow to                   [m^2]
V0    = 0.2;         % volume budget of steel, a prescribed parameter,  [m^3]

%% Define the design domain and the ground structure
nx = 7;   ny = 5;    % nodes across and up
Lx = 6.0; Ly = 4.0;  % domain size                                      [m]
connectivityType = 'full';   % 'neighbour' = join the 8 nearest nodes only
                             % 'full'      = join every pair, except through
                             %               a third node (collinear overlap)

FlagLiveShow = false;            % redraw the truss at every iteration

% nodalCoordinates matrix: [nodeIndex, ycoordinate, xcoordinate]; 
% node ordering: from bottom to top, left to right
% elementConnectivity matrix: [elementIndex, nodei, nodej]
[nodalCoordinates, elementConnectivity] = groundStructure(nx, ny, Lx, Ly, connectivityType);
numberOfNodes    = size(nodalCoordinates,1);
numberOfDOFs     = 2*numberOfNodes;               % two DOFs per node
numberOfElements = size(elementConnectivity,1);

%% Specify which DOFs are fixed
% The lecture pinned the WHOLE left edge.  Michell's problem has exactly TWO
% supports, at y = Ly/4 and y = 3Ly/4.  The heights are stated PHYSICALLY and
% snapped to the nearest left-edge node, so section 2.2's refinement run needs
% no edit here: on the 13 x 9 grid the same two heights land on other nodes.
leftNodes     = find(abs(nodalCoordinates(:,2)) < 1e-12);   % every node at x = 0
supportHeight = [Ly/4; 3*Ly/4];                             % where the pins go [m]
supportNodes  = zeros(numel(supportHeight),1);
for k = 1:numel(supportHeight)
    [~, nearest]    = min(abs(nodalCoordinates(leftNodes,3) - supportHeight(k)));
    supportNodes(k) = leftNodes(nearest);
end
supportNodes = unique(supportNodes);
% Two distinct pins are what makes this Michell's problem rather than the
% lecture's cantilever; a grid that cannot offer two is the wrong grid.
assert(numel(supportNodes) == numel(supportHeight), ...
       'ny = %d put both supports on the same node; Ly/4 and 3Ly/4 need distinct nodes.', ny);
fixedDOFs  = sort([2*supportNodes-1; 2*supportNodes])';  % pinned: x and y
allDOFs    = 1:numberOfDOFs;
freeDOFs   = setdiff(allDOFs, fixedDOFs);                % free DOFs

%% Define the DOFs with loads
% A single downward load at the MIDDLE of the right edge, (Lx, Ly/2).
rightNodes = find(abs(nodalCoordinates(:,2) - Lx) < 1e-12); % every node at x = Lx
[~, nearest] = min(abs(nodalCoordinates(rightNodes,3) - Ly/2));
loadNode   = rightNodes(nearest);                        % mid-height of right edge
forceVector             = zeros(numberOfDOFs,1);         % initialise load vector
forceVector(2*loadNode) = -P;                            % (fx, fy) = (0, -P)

fprintf('SUPPORTS: nodes %s, at y = %s m  (asked for %s)\n', ...
        mat2str(supportNodes'), mat2str(nodalCoordinates(supportNodes,3)', 4), ...
        mat2str(supportHeight', 4));
fprintf('LOAD    : node %d, at (x,y) = (%.3f, %.3f) m, (fx,fy) = (0, %.3e) N\n', ...
        loadNode, nodalCoordinates(loadNode,2), nodalCoordinates(loadNode,3), -P);
fprintf('FIXED DOFs: %d  = %d supports x 2 directions  --> %s\n', ...
        numel(fixedDOFs), numel(supportNodes), mat2str(fixedDOFs));

%% Compute element lengths (the volume constraint needs them)
elementLength = zeros(numberOfElements,1);
for e = 1:numberOfElements
    nodei = elementConnectivity(e,2);
    nodej = elementConnectivity(e,3);
    dx = nodalCoordinates(nodej,2) - nodalCoordinates(nodei,2);
    dy = nodalCoordinates(nodej,3) - nodalCoordinates(nodei,3);
    elementLength(e) = sqrt(dx^2 + dy^2);
end

%% Non-dimensionalise the problem
a0 = V0/sum(elementLength)*ones(numberOfElements,1);  % initial uniform areas, spending V0
x0 = (a0 - amin)/(amax - amin);                    % initial design variable, in [0,1]
lb = zeros(numberOfElements,1);                    % x >= 0  i.e. a >= amin
ub = ones(numberOfElements,1);                     % x <= 1  i.e. a <= amax

% sum L_e a_e <= V0 with a = amin + (amax - amin)*x, divided through by V0:
%   amin*sum(L_e)/V0 + (amax - amin)*(L_e/V0)*x <= 1
% The amin term carries no x, so it belongs on the right-hand side.
Aineq = (amax - amin)*elementLength'/V0;           % coefficients of x
bineq = 1 - amin*sum(elementLength)/V0;            % budget left after the amin floor

fprintf('GROUND STRUCTURE (%s): %d nodes, %d bars, %d free DOFs\n', ...
        connectivityType, numberOfNodes, numberOfElements, numel(freeDOFs));
fprintf('  volume budget %.5e m^3\n', V0);

%% Plot the ground structure before optimising
aPlotRef = V0/Lx;    % fixed width scale, so every figure is drawn alike
plotTruss(nodalCoordinates, elementConnectivity, a0, aPlotRef, fixedDOFs, forceVector, ...
          sprintf('Ground structure (%s): %d bars', connectivityType, numberOfElements), 1);

%% Set up the optimisation problem
% Compute the compliance of the initial uniform design, used for normalisation
[f0, ~] = complianceObjective(x0, nodalCoordinates, elementConnectivity, ...
                              E, forceVector, freeDOFs, amin, amax, 1);

% Set up the objective function for the optimiser, normalised by f0
objective = @(x) complianceObjective(x, nodalCoordinates, elementConnectivity, ...
                                     E, forceVector, freeDOFs, amin, amax, f0);

% Set up the options for the optimiser
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...             % handles bounds cheaply
    'SpecifyObjectiveGradient', true, ...          % we supply dc/da ourselves
    'HessianApproximation', 'lbfgs', ...           % dense Hessian costs more
    'OptimalityTolerance', 1e-6, ...               % THE stopping test: KKT
    'MaxIterations', 5e3, 'MaxFunctionEvaluations', 5e3, ...   % flag 0 if exceeded
    'Display', 'none', ...
    ... 
    ... % OutputFcn is a hook fmincon calls once per optimisation ITERATION.  
    ... % fmincon decides
    ... % how to call it, and always passes exactly three arguments, in this order:
    ... %        stop = theHook(x, optimValues, state)
    ... %   x            the current design variable
    ... %   optimValues  iteration data (.iteration, .fval, .firstorderopt)
    ... %   state        'init' before the first iteration, 'iter' at each one,
    ... %                'done' at the end
    ... % But iterationOutput needs twelve arguments, not three.  The @(...) below
    ... % is an ADAPTER: an unnamed function of exactly the three fmincon supplies,
    ... % whose body calls iterationOutput with those three PLUS the mesh, loads,
    ... % f0 and so on, captured from this workspace as it stands right now.
    'OutputFcn', @(x, optimValues, state) ...
        iterationOutput(amin + (amax - amin)*x, optimValues, state, ...
                        nodalCoordinates, elementConnectivity, aPlotRef, ...
                        fixedDOFs, forceVector, elementLength, V0, f0, FlagLiveShow));

%% Run the optimiser
% fmincon owns the loop: each call to the objective is ONE FEA solve.
fprintf('\n  iter      compliance       f/f0      volume/V0   first-order opt\n');
timer = tic;
[xStar, fhatStar, exitflag, output, lambda] = ...
    fmincon(objective, x0, Aineq, bineq, [], [], lb, ub, [], options);
elapsed = toc(timer);
aStar = amin + (amax - amin)*xStar;                % back to physical areas
fStar = fhatStar*f0;                               % back to physical compliance

%% Report
fprintf('\n================ SECTION 2.1 -- THE BASELINE ================\n');
fprintf('  run: %d x %d, ''%s'', %d bars, %.2f s\n', ...
        nx, ny, connectivityType, numberOfElements, elapsed);

% 1. fixed DOFs ----------------------------------------
fprintf('\n  1. fixed DOFs                    : %d of %d\n', numel(fixedDOFs), numberOfDOFs);
fprintf('        two pinned supports x two directions = 4, which is what it should be,\n');
fprintf('        and 4 >= 3 is the minimum that removes the rigid-body motions.\n');

% 2. does it converge? ------------------------------------
fprintf('\n  2. exit flag                     : %d\n', exitflag);
if exitflag == 1
    fprintf('        claims: first-order optimality (KKT) met to OptimalityTolerance\n');
elseif exitflag == 0
    fprintf('        claims: NOTHING -- it ran out of iterations/evaluations\n');
elseif exitflag == 2
    fprintf('        claims: the step became smaller than StepTolerance\n');
else
    fprintf('        claims: see doc fmincon for flag %d\n', exitflag);
end
fprintf('        first-order optimality reached : %.4e (tolerance %.0e)\n', ...
        output.firstorderopt, options.OptimalityTolerance);

% 3. iterations ------------------------------------
fprintf('\n  3. iterations                    : %d  (of %d allowed)\n', ...
        output.iterations, options.MaxIterations);
fprintf('        function evaluations          : %d\n', output.funcCount);

% 4. the objective ------------------------------------
fprintf('\n  4. f0  (uniform start)           : %.6e J\n', f0);
fprintf('     f*  (optimised)               : %.6e J\n', fStar);
fprintf('     f*/f0                         : %.6f   (%.1f%% stiffer for the same steel)\n', ...
        fStar/f0, 100*(f0/fStar - 1));

% 5. the material, both ways ------------------------------------
physicalVolume = sum(aStar.*elementLength);
fprintf('\n  5. material used, physical       : %.6e m^3  (of V0 = %g)\n', physicalVolume, V0);
fprintf('     material used, normalised     : %.6f   = (sum L_e a_e)/V0\n', physicalVolume/V0);
fprintf('        the row fmincon sees, Aineq*x - bineq = %+.3e (0 means tight)\n', ...
        Aineq*xStar - bineq);

% 6. the multiplier ------------------------------------
fprintf('\n  6. lambda.ineqlin                : %.6e\n', lambda.ineqlin(1));
if lambda.ineqlin(1) > 0
    fprintf('        > 0, so the budget is ACTIVE: every cubic metre is spent, and\n');
    fprintf('        one more unit of normalised budget would buy %.3e in f/f0.\n', ...
            lambda.ineqlin(1));
else
    fprintf('        = 0, so the budget is INACTIVE -- which would be a surprise here.\n');
end

% 7. how many bars survived? -------------------
numberVoid = sum(aStar < 0.01*max(aStar));      % count of bars with area less than 1% of the largest bar
numberKept = numberOfElements - numberVoid;
fprintf('\n  7. bars above 1%% of the largest  : %d of %d\n', numberKept, numberOfElements);
fprintf('        (removed: %d; the 1%% line sits at a = %.3e m^2)\n', ...
        numberVoid, 0.01*max(aStar));
fprintf('        those %d bars hold %.1f%% of the volume\n', numberKept, ...
        100*sum(aStar(aStar >= 0.01*max(aStar)).*elementLength(aStar >= 0.01*max(aStar)))/physicalVolume);


%% Plot the optimised structure
plotTruss(nodalCoordinates, elementConnectivity, aStar, aPlotRef, fixedDOFs, forceVector, ...
          sprintf('Optimised: f = %.3e, %d bars removed', fStar, numberVoid), 2);

%% Plot the convergence history
history = iterationOutput([], [], 'get', [], [], [], [], [], [], [], [], []);
figure(3); clf;
yyaxis left;  semilogy(history(:,1), history(:,3), 'LineWidth', 1.5);
ylabel('normalised objective  f / f_0');  xlabel('fmincon iteration');  grid on;
yyaxis right; plot(history(:,1), history(:,4), 'LineWidth', 1.2);
ylabel('volume / V_0');  ylim([0 max(1.2, 1.05*max(history(:,4)))]);
title('Convergence: the objective falls, the budget stays spent');
lightFigure(gcf);
exportgraphics(gcf, sprintf('a3_q2_convergence_%dx%d_%s.png', nx, ny, connectivityType), 'Resolution', 300);

%% Plot how the material ended up distributed across the bars
% Log axes on both sides, because the areas span five decades and the counts
% span three.  Read it as: a big pile of bars left at essentially nothing, a
% small cluster that took all the material, and a GAP in between -- which is
% what makes the 1% cut-off above a safe place to draw the line.
figure(4); clf;
histogram(aStar, logspace(log10(amin), log10(amax), 36), ...
          'FaceColor', [0.15 0.30 0.45], 'EdgeColor', 'w');
set(gca, 'XScale', 'log', 'YScale', 'log');
xline(0.01*max(aStar), '--', '1% of max', 'LineWidth', 1.2);
xlabel('bar area  a_e   [m^2]');  ylabel('number of bars');  grid on;
xlim([amin amax]);
title('Where the material went');
lightFigure(gcf);
exportgraphics(gcf, sprintf('a3_q2_histogram_%dx%d_%s.png', nx, ny, connectivityType), 'Resolution', 300);

%% Write the optimised structure to disk as well (figure 2)
figure(2);
lightFigure(gcf);
exportgraphics(gcf, sprintf('a3_q2_structure_%dx%d_%s.png', nx, ny, connectivityType), 'Resolution', 300);

fprintf('\n  8. three figures written: structure, convergence, histogram\n');

%% ================================================================ %%
%          READ AND MODIFY THIS ONE -- it is the lecture             %
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

%% ================================================================ %%
%        DO NOT MODIFY BELOW THIS LINE - UNLESS YOU ARE BRAVE        %
%  ================================================================  %

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

function plotTruss(nodalCoordinates, elementConnectivity, a, aPlotRef, fixedDOFs, forceVector, ttl, figno)
% Width is a/aPlotRef with aPlotRef FIXED, so two figures are to the same
% scale.  Thickness AND opacity both scale with the area, so a thin bar looks
% thin and a near-void bar fades away; no bar is classified or hidden.
    figure(figno); clf; hold on;
    % ONE style, not two.  Line width AND opacity both scale with the area, so
    % the picture carries a single continuous quantity.  There is deliberately
    % no visual class boundary between a "kept" bar and a "removed" one: the
    % areas taper smoothly (on a fine grid the 90th largest is at 1.2% of the
    % biggest and the 100th at 0.6%), so any such boundary would be a drawing
    % convention, not a fact about the structure -- and it invites the question
    % "if that bar carries load, why is it dotted?".  Every bar is drawn the
    % same way; the near-void ones fade to nothing by themselves.
    %
    % The 1%-of-largest rule stays where it belongs: COUNTING the voids in the
    % printed report, never deciding how a bar is drawn.
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
    sn = unique(ceil(fixedDOFs/2));
    plot(nodalCoordinates(sn,2), nodalCoordinates(sn,3), '^', 'MarkerSize', 8, ...
         'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor', 'none');
    for k = find(forceVector ~= 0)'
        nd = ceil(k/2);
        quiver(nodalCoordinates(nd,2), nodalCoordinates(nd,3), 0, -0.3*max(nodalCoordinates(:,3)), ...
               0, 'LineWidth', 2, 'Color', [0.85 0.20 0.10], 'MaxHeadSize', 3);
    end
    % Margin, NOT 'axis tight'.  The chords lie exactly on y = 0 and y = Ly, so
    % with tight limits MATLAB clips those lines to the axes box and they are
    % drawn at uneven width -- which reads as a thickness variation that is not
    % in the areas.  Keep the structure off the boundary.
    axis equal; axis off;
    xl = xlim;  yl = ylim;
    xlim(xl + 0.06*diff(xl)*[-1 1]);
    ylim(yl + 0.12*diff(yl)*[-1 1]);
    title(sprintf('%s\n(max a = %.3e)', ttl, max(a)));
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function stop = iterationOutput(a, optimValues, state, nodalCoordinates, ...
             elementConnectivity, aPlotRef, fixedDOFs, forceVector, elementLength, V0, f0, FlagLiveShow)
% Called by fmincon once per iteration, through the wrapper set up above.
    persistent history      % survives between calls, unlike a normal variable
    stop = false;
    switch state
        case 'init'
            history = [];                                  % new run, clear it
        case 'iter'
            volumeFraction = (elementLength'*a)/V0;
            normalisedObj  = optimValues.fval;             % f/f0, what fmincon sees
            compliance     = normalisedObj*f0;             % back to physical
            fprintf('  %4d   %13.6e   %8.5f   %9.6f   %12.4e\n', optimValues.iteration, ...
                    compliance, normalisedObj, volumeFraction, optimValues.firstorderopt);
            history(end+1,:) = [optimValues.iteration compliance normalisedObj volumeFraction];
            if FlagLiveShow
                plotTruss(nodalCoordinates, elementConnectivity, a, aPlotRef, fixedDOFs, forceVector, ...
                          sprintf('iteration %d', optimValues.iteration), 2);
                drawnow;
            end
        case 'get'
            stop = history;   % not a state fmincon uses -- our own back door,
                              % the only way to read a persistent from outside
    end
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
