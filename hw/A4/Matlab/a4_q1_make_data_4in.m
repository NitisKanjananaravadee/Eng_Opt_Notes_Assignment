%% Clear everything from memory and start from scratch
clear all; close all; clc;
%% ======================================================================= %%
%  A4 TASK 1 -- MAKE THE FOUR-INPUT TRAINING SET
%
%  This is L07's l07_4_make_data.m with all FOUR design parameters varying
%  instead of two.  Every sample is one complete size optimisation.
%
%  YOU DO NOT NEED TO RUN THIS.  truss_data_4in.mat is shipped beside this
%  script.  It is here so that nothing about the data is hidden from you.
%
%  WHAT CHANGED FROM THE LECTURE'S VERSION, AND WHY IT MATTERS:
%    - The depth Ly and the budget V0 now vary too, so the ground structure
%      must be rebuilt INSIDE the loop: the nodes move from sample to
%      sample.  In the lecture they were fixed and it was built once.
%    - Because the nodes move, the element lengths and the volume
%      constraint change every sample as well.
%
%  SELF-CONTAINED.  The optimisation and the ground structure are at the
%  bottom of this file, copied from L07 step 3 unchanged.
%
%  Cost: about 2 s per sample, so the settings below take about half an hour.
%% ======================================================================= %%

%% Define the problem parameters
E     = 2.1e11;      % Young's modulus                                  [Pa]
P     = 1e5;         % applied load                                     [N]
amin  = 1e-7;        % smallest area a bar may shrink to                [m^2]
amax  = 1e-2;        % largest area a bar may grow to                   [m^2]

%% Define the design domain
nx = 7;   ny = 5;    % nodes across and up -> 106 bars
Lx = 6.0;            % span of the design domain                        [m]
connectivityType = 'neighbour';

%% What varies from sample to sample, and over what range
%  ALL FOUR.  p = [Ly, V0, loadHeight, loadAngle]
rangeLow  = [2.0, 0.06, 0.00, -75*pi/180];
rangeHigh = [6.0, 0.40, 1.00, +75*pi/180];
% The load can only sit on a NODE.  The right edge has ny of them, at
% fractions 0, 1/4, ... 1 of Ly, so the height input takes exactly those
% values -- sampled, not snapped.  Everything else is genuinely continuous.
loadHeightLevels = (0:ny-1)/(ny-1);
numberOfSamples = 2000;

%% How many bars?  The connectivity does not depend on Ly, so one probe
%  tells us, and the answer is the number of network outputs.
[probeNodes, probeElements] = groundStructure(nx, ny, Lx, 4.0, connectivityType);
numberOfElements = size(probeElements,1);
numberOfNodes    = size(probeNodes,1);
numberOfDOFs     = 2*numberOfNodes;
fprintf('%d bars, so the network will have %d outputs\n', ...
        numberOfElements, numberOfElements);

options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient', true, ...
    'HessianApproximation', 'lbfgs', ...
    'OptimalityTolerance', 1e-6, ...
    'MaxIterations', 5e3, 'MaxFunctionEvaluations', 5e3, ...
    'Display', 'none');

%% Sample, optimise, record
rng(21);                                  % fixed, so the file is reproducible
Depth  = zeros(numberOfSamples, 1);       % truss depth Ly                   [m]
Budget = zeros(numberOfSamples, 1);       % volume budget V0                [m^3]
H      = zeros(numberOfSamples, 1);       % load height, as a fraction of Ly
Angle  = zeros(numberOfSamples, 1);       % load angle from vertical        [rad]
X      = zeros(numberOfSamples, numberOfElements);  % the DESIGN VARIABLE, in [0,1]
fStar  = zeros(numberOfSamples, 1);       % the compliance reached            [J]
keep   = false(numberOfSamples, 1);

timer = tic;
for k = 1:numberOfSamples
    p = rangeLow + (rangeHigh - rangeLow).*rand(1,4);
    p(3) = loadHeightLevels(randi(numel(loadHeightLevels)));   % on a node
    Ly = p(1);  V0 = p(2);  loadHeight = p(3);  loadAngle = p(4);

    % The nodes MOVE with Ly, so the whole ground structure is rebuilt here
    [nodalCoordinates, elementConnectivity] = ...
        groundStructure(nx, ny, Lx, Ly, connectivityType);

    leftNodes  = find(abs(nodalCoordinates(:,2)) < 1e-12);
    fixedDOFs  = sort([2*leftNodes-1; 2*leftNodes])';
    freeDOFs   = setdiff(1:numberOfDOFs, fixedDOFs);

    rightNodes = find(abs(nodalCoordinates(:,2) - Lx) < 1e-12);
    [~, nearest] = min(abs(nodalCoordinates(rightNodes,3) - loadHeight*Ly));
    loadNode = rightNodes(nearest);
    forceVector               = zeros(numberOfDOFs,1);
    forceVector(2*loadNode-1) =  P*sin(loadAngle);
    forceVector(2*loadNode)   = -P*cos(loadAngle);

    elementLength = zeros(numberOfElements,1);
    for e = 1:numberOfElements
        nodei = elementConnectivity(e,2);  nodej = elementConnectivity(e,3);
        dx = nodalCoordinates(nodej,2) - nodalCoordinates(nodei,2);
        dy = nodalCoordinates(nodej,3) - nodalCoordinates(nodei,3);
        elementLength(e) = sqrt(dx^2 + dy^2);
    end

    a0    = V0/sum(elementLength)*ones(numberOfElements,1);
    x0    = (a0 - amin)/(amax - amin);
    lb    = zeros(numberOfElements,1);
    ub    = ones(numberOfElements,1);
    Aineq = (amax - amin)*elementLength'/V0;
    bineq = 1 - amin*sum(elementLength)/V0;

    [f0, ~] = complianceObjective(x0, nodalCoordinates, elementConnectivity, ...
                                  E, forceVector, freeDOFs, amin, amax, 1);
    objective = @(x) complianceObjective(x, nodalCoordinates, elementConnectivity, ...
                                         E, forceVector, freeDOFs, amin, amax, f0);
    [xStar, fhatStar, exitflag] = ...
        fmincon(objective, x0, Aineq, bineq, [], [], lb, ub, [], options);
    if exitflag <= 0, continue; end

    Depth(k)  = Ly;
    Budget(k) = V0;
    H(k)      = loadHeight;
    Angle(k)  = loadAngle;
    X(k,:)    = xStar';       % the DESIGN VARIABLE itself, not the areas
    fStar(k)  = fhatStar*f0;
    keep(k)   = true;

    if mod(k,100) == 0
        fprintf('   %d/%d kept %d, %.1f min elapsed\n', ...
                k, numberOfSamples, sum(keep), toc(timer)/60);
    end
end
Depth = Depth(keep);  Budget = Budget(keep);  H = H(keep);
Angle = Angle(keep);  X = X(keep,:);  fStar = fStar(keep);

%% Save -- THE DATA AND NOTHING ELSE
%  Depth, Budget, H and Angle are the four inputs; X is the design variable the
%  optimiser returned for them; fStar is the compliance it reached.  The truss
%  itself is defined by the constants at the top of this file, and every script
%  that reads this data states those constants again for itself.
save('truss_data_4in.mat', 'Depth', 'Budget', 'H', 'Angle', 'X', 'fStar');

fprintf('\nkept %d of %d samples in %.1f min -> truss_data_4in.mat\n', ...
        numel(fStar), numberOfSamples, toc(timer)/60);
A = amin + (amax - amin)*X;      % recover areas, for the summary below only
% Count unused bars by L06's own criterion: below 1% of the largest area in
% THAT design.  Robust -- whether a bar lands exactly on amin or a hair above
% depends on the stopping tolerance and should not change what we report.
unused = mean(arrayfun(@(i) mean(A(i,:) <  0.01*max(A(i,:))), 1:size(A,1)));
used   = mean(arrayfun(@(i) sum( A(i,:) >= 0.01*max(A(i,:))), 1:size(A,1)));
fprintf('%.0f%% of bars carry no real material: a median of %.0f of %d do the work\n', ...
        100*unused, used, size(A,2));
fprintf('compliance ranged from %.1f to %.1f J\n', min(fStar), max(fStar));

%% ================================================================ %%
%   The two functions below are COPIED FROM STEP 3, unchanged.       %
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
