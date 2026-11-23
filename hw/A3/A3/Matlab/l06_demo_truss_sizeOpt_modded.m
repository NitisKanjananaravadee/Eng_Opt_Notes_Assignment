%% Clear everything from memory and start from scratch
clear all; close all; clc;

%% ======================================================================= %%
%  L06 -- TRUSS Size Optimisation: FEA + adjoint sensitivity + optimiser, in one loop
%
%  Engineering Optimisation (M-EOP), Lecture 6.
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

%% Define the design domain and the ground structure(basedline)
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
supportNodes = [1 + (ny-1)/4;  1 + 3*(ny-1)/4];
leftNodes  = find(abs(nodalCoordinates(:,2)) < 1e-12);   % every node at x = 0
fixedDOFs  = sort([2*supportNodes-1; 2*supportNodes])';        % pinned: x and y
allDOFs    = 1:numberOfDOFs;
freeDOFs   = setdiff(allDOFs, fixedDOFs);                % free DOFs

%% Define the DOFs with loads
rightNodes = find(abs(nodalCoordinates(:,2) - Lx) < 1e-12); % every node at x = Lx
loadNode   = rightNodes(ceil(numel(rightNodes)/2));      % middle of right edge
forceVector          = zeros(numberOfDOFs,1);                      % initialise load vector
forceVector(2*loadNode) = -P;                                      % downwards

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

historyBaseline = iterationOutput([], [], 'get', [], [], [], [], [], [], [], [], []);

%% Define the design domain and the ground structure(A)
nxA = 7;   nyA = 5;    % nodes across and up
LxA = 6.0; LyA = 4.0;  % domain size                                      [m]
connectivityTypeA = 'neighbour';   % 'neighbour' = join the 8 nearest nodes only
                             % 'full'      = join every pair, except through
                             %               a third node (collinear overlap)

FlagLiveShowA = false;            % redraw the truss at every iteration

% nodalCoordinates matrix: [nodeIndex, ycoordinate, xcoordinate];
% node ordering: from bottom to top, left to right
% elementConnectivity matrix: [elementIndex, nodei, nodej]
[nodalCoordinatesA, elementConnectivityA] = groundStructure(nxA, nyA, LxA, LyA, connectivityTypeA);
numberOfNodesA    = size(nodalCoordinatesA,1);
numberOfDOFsA     = 2*numberOfNodesA;               % two DOFs per node
numberOfElementsA = size(elementConnectivityA,1);

%% Specify which DOFs are fixed
supportNodesA = [1 + (nyA-1)/4;  1 + 3*(nyA-1)/4];
leftNodesA  = find(abs(nodalCoordinatesA(:,2)) < 1e-12);   % every node at x = 0
fixedDOFsA  = sort([2*supportNodesA-1; 2*supportNodesA])';        % pinned: x and y
allDOFsA    = 1:numberOfDOFsA;
freeDOFsA   = setdiff(allDOFsA, fixedDOFsA);                % free DOFs

%% Define the DOFs with loads
rightNodesA = find(abs(nodalCoordinatesA(:,2) - LxA) < 1e-12); % every node at x = LxA
loadNodeA   = rightNodesA(ceil(numel(rightNodesA)/2));      % middle of right edge
forceVectorA          = zeros(numberOfDOFsA,1);                      % initialise load vector
forceVectorA(2*loadNodeA) = -P;                                      % downwards

%% Compute element lengths (the volume constraint needs them)
elementLengthA = zeros(numberOfElementsA,1);
for e = 1:numberOfElementsA
    nodei = elementConnectivityA(e,2);
    nodej = elementConnectivityA(e,3);
    dx = nodalCoordinatesA(nodej,2) - nodalCoordinatesA(nodei,2);
    dy = nodalCoordinatesA(nodej,3) - nodalCoordinatesA(nodei,3);
    elementLengthA(e) = sqrt(dx^2 + dy^2);
end

%% Non-dimensionalise the problem
a0A = V0/sum(elementLengthA)*ones(numberOfElementsA,1);  % initial uniform areas, spending V0
x0A = (a0A - amin)/(amax - amin);                    % initial design variable, in [0,1]
lbA = zeros(numberOfElementsA,1);                    % x >= 0  i.e. a >= amin
ubA = ones(numberOfElementsA,1);                     % x <= 1  i.e. a <= amax

% sum L_e a_e <= V0 with a = amin + (amax - amin)*x, divided through by V0:
%   amin*sum(L_e)/V0 + (amax - amin)*(L_e/V0)*x <= 1
% The amin term carries no x, so it belongs on the right-hand side.
AineqA = (amax - amin)*elementLengthA'/V0;           % coefficients of x
bineqA = 1 - amin*sum(elementLengthA)/V0;            % budget left after the amin floor

fprintf('GROUND STRUCTURE (%s): %d nodes, %d bars, %d free DOFs\n', ...
        connectivityTypeA, numberOfNodesA, numberOfElementsA, numel(freeDOFsA));
fprintf('  volume budget %.5e m^3\n', V0);

%% Plot the ground structure before optimising
aPlotRefA = V0/LxA;    % fixed width scale, so every figure is drawn alike
plotTruss(nodalCoordinatesA, elementConnectivityA, a0A, aPlotRefA, fixedDOFsA, forceVectorA, ...
          sprintf('Run A ground structure (%s): %d bars', connectivityTypeA, numberOfElementsA), 5);

%% Set up the optimisation problem
% Compute the compliance of the initial uniform design, used for normalisation
[f0A, ~] = complianceObjective(x0A, nodalCoordinatesA, elementConnectivityA, ...
                              E, forceVectorA, freeDOFsA, amin, amax, 1);

% Set up the objective function for the optimiser, normalised by f0
objectiveA = @(x) complianceObjective(x, nodalCoordinatesA, elementConnectivityA, ...
                                     E, forceVectorA, freeDOFsA, amin, amax, f0A);

% Set up the options for the optimiser
optionsA = optimoptions('fmincon', ...
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
                        nodalCoordinatesA, elementConnectivityA, aPlotRefA, ...
                        fixedDOFsA, forceVectorA, elementLengthA, V0, f0A, FlagLiveShowA));

%% Run the optimiser
% fmincon owns the loop: each call to the objective is ONE FEA solve.
fprintf('\n  iter      compliance       f/f0       volume/V0   first-order opt\n');
timerA = tic;
[xStarA, fhatStarA, exitflagA, outputA, lambdaA] = ...
    fmincon(objectiveA, x0A, AineqA, bineqA, [], [], lbA, ubA, [], optionsA);
elapsedA = toc(timerA);
aStarA = amin + (amax - amin)*xStarA;                % back to physical areas
fStarA = fhatStarA*f0A;                               % back to physical compliance
historyA = iterationOutput([], [], 'get', [], [], [], [], [], [], [], [], []);

%% Define the design domain and the ground structure(B)
nxB = 13;  nyB = 9;    % nodes across and up
LxB = 6.0; LyB = 4.0;  % domain size                                      [m]
connectivityTypeB = 'full';   % 'neighbour' = join the 8 nearest nodes only
                             % 'full'      = join every pair, except through
                             %               a third node (collinear overlap)

FlagLiveShowB = false;            % redraw the truss at every iteration

% nodalCoordinates matrix: [nodeIndex, ycoordinate, xcoordinate];
% node ordering: from bottom to top, left to right
% elementConnectivity matrix: [elementIndex, nodei, nodej]
[nodalCoordinatesB, elementConnectivityB] = groundStructure(nxB, nyB, LxB, LyB, connectivityTypeB);
numberOfNodesB    = size(nodalCoordinatesB,1);
numberOfDOFsB     = 2*numberOfNodesB;               % two DOFs per node
numberOfElementsB = size(elementConnectivityB,1);

%% Specify which DOFs are fixed
supportNodesB = [1 + (nyB-1)/4;  1 + 3*(nyB-1)/4];
leftNodesB  = find(abs(nodalCoordinatesB(:,2)) < 1e-12);   % every node at x = 0
fixedDOFsB  = sort([2*supportNodesB-1; 2*supportNodesB])';        % pinned: x and y
allDOFsB    = 1:numberOfDOFsB;
freeDOFsB   = setdiff(allDOFsB, fixedDOFsB);                % free DOFs

%% Define the DOFs with loads
rightNodesB = find(abs(nodalCoordinatesB(:,2) - LxB) < 1e-12); % every node at x = LxB
loadNodeB   = rightNodesB(ceil(numel(rightNodesB)/2));      % middle of right edge
forceVectorB          = zeros(numberOfDOFsB,1);                      % initialise load vector
forceVectorB(2*loadNodeB) = -P;                                      % downwards

%% Compute element lengths (the volume constraint needs them)
elementLengthB = zeros(numberOfElementsB,1);
for e = 1:numberOfElementsB
    nodei = elementConnectivityB(e,2);
    nodej = elementConnectivityB(e,3);
    dx = nodalCoordinatesB(nodej,2) - nodalCoordinatesB(nodei,2);
    dy = nodalCoordinatesB(nodej,3) - nodalCoordinatesB(nodei,3);
    elementLengthB(e) = sqrt(dx^2 + dy^2);
end

%% Non-dimensionalise the problem
a0B = V0/sum(elementLengthB)*ones(numberOfElementsB,1);  % initial uniform areas, spending V0
x0B = (a0B - amin)/(amax - amin);                    % initial design variable, in [0,1]
lbB = zeros(numberOfElementsB,1);                    % x >= 0  i.e. a >= amin
ubB = ones(numberOfElementsB,1);                     % x <= 1  i.e. a <= amax

% sum L_e a_e <= V0 with a = amin + (amax - amin)*x, divided through by V0:
%   amin*sum(L_e)/V0 + (amax - amin)*(L_e/V0)*x <= 1
% The amin term carries no x, so it belongs on the right-hand side.
AineqB = (amax - amin)*elementLengthB'/V0;           % coefficients of x
bineqB = 1 - amin*sum(elementLengthB)/V0;            % budget left after the amin floor

fprintf('GROUND STRUCTURE (%s): %d nodes, %d bars, %d free DOFs\n', ...
        connectivityTypeB, numberOfNodesB, numberOfElementsB, numel(freeDOFsB));
fprintf('  volume budget %.5e m^3\n', V0);

%% Plot the ground structure before optimising
aPlotRefB = V0/LxB;    % fixed width scale, so every figure is drawn alike
plotTruss(nodalCoordinatesB, elementConnectivityB, a0B, aPlotRefB, fixedDOFsB, forceVectorB, ...
          sprintf('Run B ground structure (%s): %d bars', connectivityTypeB, numberOfElementsB), 9);

%% Set up the optimisation problem
% Compute the compliance of the initial uniform design, used for normalisation
[f0B, ~] = complianceObjective(x0B, nodalCoordinatesB, elementConnectivityB, ...
                              E, forceVectorB, freeDOFsB, amin, amax, 1);

% Set up the objective function for the optimiser, normalised by f0
objectiveB = @(x) complianceObjective(x, nodalCoordinatesB, elementConnectivityB, ...
                                     E, forceVectorB, freeDOFsB, amin, amax, f0B);

% Set up the options for the optimiser
optionsB = optimoptions('fmincon', ...
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
                        nodalCoordinatesB, elementConnectivityB, aPlotRefB, ...
                        fixedDOFsB, forceVectorB, elementLengthB, V0, f0B, FlagLiveShowB));

%% Run the optimiser
% fmincon owns the loop: each call to the objective is ONE FEA solve.
fprintf('\n  iter      compliance       f/f0       volume/V0   first-order opt\n');
timerB = tic;
[xStarB, fhatStarB, exitflagB, outputB, lambdaB] = ...
    fmincon(objectiveB, x0B, AineqB, bineqB, [], [], lbB, ubB, [], optionsB);
elapsedB = toc(timerB);
aStarB = amin + (amax - amin)*xStarB;                % back to physical areas
fStarB = fhatStarB*f0B;                               % back to physical compliance
historyB = iterationOutput([], [], 'get', [], [], [], [], [], [], [], [], []);


%% Report
fprintf('------------------------Baseline--------------------------\n')
fprintf('SUPPORTS: nodes %d and %d, %d fixed DOFs of %d\n in Baseline and run A', ...
        supportNodes(1), supportNodes(2), numel(fixedDOFs), numberOfDOFs);
% Does it converge? ------------------------------------
fprintf('  exit flag: %d\n', exitflag);
if exitflag == 1
    fprintf('  status: converged\n');
elseif exitflag == 0
    fprintf('  status: maximum iterations/function evaluations exceeded\n');
elseif exitflag == 2
    fprintf('  status: stopped moving\n');
else
    fprintf('  status: unknown\n');
end
fprintf('iterations: %d   function evaluations: %d\n', ...
        output.iterations, output.funcCount);
fprintf('first-order optimality: %.4e  (tolerance 1e-6)\n', output.firstorderopt);
fprintf('elapsed: %.2f s\n', elapsed);
fprintf('f0 = %.6e J   (uniform start)\n', f0);
fprintf('f* = %.6e J\n', fStar);
fprintf('f*/f0 = %.6f\n', fhatStar);
volumePhysical = sum(aStar.*elementLength);
fprintf('volume(physical): %.6e m^3 (budget V0 = %.6e)\n', volumePhysical, V0);
fprintf('volume/V0: %.3e\n', sum(aStar.*elementLength)/V0);
fprintf('volume constraint Lagrange multiplier: %.3e\n', lambda.ineqlin(1));
numberVoid = sum(aStar < 0.01*max(aStar));
fprintf('bars removed (a < 1%% of max): %d of %d\n', numberVoid, numberOfElements);

fprintf('\n------------------------ RUN A ------------------------\n');
fprintf('SUPPORTS: nodes %d and %d, %d fixed DOFs of %d\n', ...
        supportNodesA(1), supportNodesA(2), numel(fixedDOFsA), numberOfDOFsA);
% Does it converge? ------------------------------------
fprintf('  exit flag: %d\n', exitflagA);
if exitflagA == 1
    fprintf('  status: converged\n');
elseif exitflagA == 0
    fprintf('  status: maximum iterations/function evaluations exceeded\n');
elseif exitflagA == 2
    fprintf('  status: stopped moving\n');
else
    fprintf('  status: unknown\n');
end
fprintf('iterations: %d   function evaluations: %d\n', ...
        outputA.iterations, outputA.funcCount);
fprintf('first-order optimality: %.4e  (tolerance 1e-6)\n', outputA.firstorderopt);
fprintf('elapsed: %.2f s\n', elapsedA);
fprintf('f0 = %.6e J   (uniform start)\n', f0A);
fprintf('f* = %.6e J\n', fStarA);
fprintf('f*/f0 = %.6f\n', fhatStarA);
volumePhysicalA = sum(aStarA.*elementLengthA);
fprintf('volume(physical): %.6e m^3 (budget V0 = %.6e)\n', volumePhysicalA, V0);
fprintf('volume/V0: %.3e\n', sum(aStarA.*elementLengthA)/V0);
fprintf('volume constraint Lagrange multiplier: %.3e\n', lambdaA.ineqlin(1));
numberVoidA = sum(aStarA < 0.01*max(aStarA));
fprintf('bars removed (a < 1%% of max): %d of %d\n', numberVoidA, numberOfElementsA);

fprintf('\n------------------------ RUN B ------------------------\n');
fprintf('SUPPORTS: nodes %d and %d, %d fixed DOFs of %d\n', ...
        supportNodesB(1), supportNodesB(2), numel(fixedDOFsB), numberOfDOFsB);
% Does it converge? ------------------------------------
fprintf('  exit flag: %d\n', exitflagB);
if exitflagB == 1
    fprintf('  status: converged\n');
elseif exitflagB == 0
    fprintf('  status: maximum iterations/function evaluations exceeded\n');
elseif exitflagB == 2
    fprintf('  status: stopped moving\n');
else
    fprintf('  status: unknown\n');
end
fprintf('iterations: %d   function evaluations: %d\n', ...
        outputB.iterations, outputB.funcCount);
fprintf('first-order optimality: %.4e  (tolerance 1e-6)\n', outputB.firstorderopt);
fprintf('elapsed: %.2f s\n', elapsedB);
fprintf('f0 = %.6e J   (uniform start)\n', f0B);
fprintf('f* = %.6e J\n', fStarB);
fprintf('f*/f0 = %.6f\n', fhatStarB);
volumePhysicalB = sum(aStarB.*elementLengthB);
fprintf('volume(physical): %.6e m^3 (budget V0 = %.6e)\n', volumePhysicalB, V0);
fprintf('volume/V0: %.3e\n', sum(aStarB.*elementLengthB)/V0);
fprintf('volume constraint Lagrange multiplier: %.3e\n', lambdaB.ineqlin(1));
numberVoidB = sum(aStarB < 0.01*max(aStarB));
fprintf('bars removed (a < 1%% of max): %d of %d\n', numberVoidB, numberOfElementsB);


%% Plot the optimised structure
plotTruss(nodalCoordinates, elementConnectivity, aStar, aPlotRef, fixedDOFs, forceVector, ...
          sprintf('Optimised: f = %.3e, %d bars removed', fStar, numberVoid), 2);

%% Plot the convergence history
figure(3); clf;
yyaxis left;  semilogy(historyBaseline(:,1), historyBaseline(:,3), 'LineWidth', 1.5);
ylabel('normalised objective  f / f_0');  xlabel('fmincon iteration');  grid on;
yyaxis right; plot(historyBaseline(:,1), historyBaseline(:,4), 'LineWidth', 1.2);
ylabel('volume / V_0');  ylim([0 max(1.2, 1.05*max(historyBaseline(:,4)))]);
title('Convergence: the objective falls, the budget stays spent');

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


%% Plot the optimised structure (run A)
plotTruss(nodalCoordinatesA, elementConnectivityA, aStarA, aPlotRefA, fixedDOFsA, forceVectorA, ...
          sprintf('Run A optimised: f = %.3e, %d bars removed', fStarA, numberVoidA), 6);

%% Plot the convergence history (run A)
figure(7); clf;
yyaxis left;  semilogy(historyA(:,1), historyA(:,3), 'LineWidth', 1.5);
ylabel('normalised objective  f / f_0');  xlabel('fmincon iteration');  grid on;
yyaxis right; plot(historyA(:,1), historyA(:,4), 'LineWidth', 1.2);
ylabel('volume / V_0');  ylim([0 max(1.2, 1.05*max(historyA(:,4)))]);
title('Run A convergence: the objective falls, the budget stays spent');

%% Plot how the material ended up distributed across the bars (run A)
figure(8); clf;
histogram(aStarA, logspace(log10(amin), log10(amax), 36), ...
          'FaceColor', [0.15 0.30 0.45], 'EdgeColor', 'w');
set(gca, 'XScale', 'log', 'YScale', 'log');
xline(0.01*max(aStarA), '--', '1% of max', 'LineWidth', 1.2);
xlabel('bar area  a_e   [m^2]');  ylabel('number of bars');  grid on;
xlim([amin amax]);
title('Run A: where the material went');

%% Plot the optimised structure (run B)
plotTruss(nodalCoordinatesB, elementConnectivityB, aStarB, aPlotRefB, fixedDOFsB, forceVectorB, ...
          sprintf('Run B optimised: f = %.3e, %d bars removed', fStarB, numberVoidB), 10);

%% Plot the convergence history (run B)
figure(11); clf;
yyaxis left;  semilogy(historyB(:,1), historyB(:,3), 'LineWidth', 1.5);
ylabel('normalised objective  f / f_0');  xlabel('fmincon iteration');  grid on;
yyaxis right; plot(historyB(:,1), historyB(:,4), 'LineWidth', 1.2);
ylabel('volume / V_0');  ylim([0 max(1.2, 1.05*max(historyB(:,4)))]);
title('Run B convergence: the objective falls, the budget stays spent');

%% Plot how the material ended up distributed across the bars (run B)
figure(12); clf;
histogram(aStarB, logspace(log10(amin), log10(amax), 36), ...
          'FaceColor', [0.15 0.30 0.45], 'EdgeColor', 'w');
set(gca, 'XScale', 'log', 'YScale', 'log');
xline(0.01*max(aStarB), '--', '1% of max', 'LineWidth', 1.2);
xlabel('bar area  a_e   [m^2]');  ylabel('number of bars');  grid on;
xlim([amin amax]);
title('Run B: where the material went');

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
