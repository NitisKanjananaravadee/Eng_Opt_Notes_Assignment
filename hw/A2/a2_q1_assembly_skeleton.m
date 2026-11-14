%% A2 Q1 SKELETON -- FEA ASSEMBLY
%
% Engineering Optimisation (M-EOP), Assignment 2, Task 1, part 4.
%
% No optimisation appears anywhere in this script. 

clear; clc

%% =======================================================================
%%  1 -- Define the problem parameters
%% =======================================================================
% The eight spring rates k1 ... k8 in N/m, and the three loads P1, P2, P3
k1 = 5000;  k2 = 1500;  k3 = 2000;  k4 = 1000;
k5 = 2500;  k6 =  500;  k7 = 3000;  k8 = 3500;
P1 = 1000;  P2 = 2000;  P3 = 3000;

%% Define the number of elements and nodes
numberOfElements = 8;
numberOfNodes    = 5;

%% =======================================================================
%%  2 -- Set up the element connectivity table
%% =======================================================================
% One row per spring, four columns:
%
%        [ element number , node i , node j , spring stiffness ]
%
% Nodes 1, 2, 3 are the three carts; node 4 is the LEFT wall and node 5 the
% RIGHT wall -- two walls, so two distinct nodes.  k7 and k8 both run from
% cart 3 to the right wall, so they are two separate elements on (3,5).
elementConnectivity = [ 1   1   4   k1       % left wall  -- cart 1
                        2   2   4   k2       % left wall  -- cart 2
                        3   3   4   k3       % left wall  -- cart 3
                        4   1   2   k4       % cart 1     -- cart 2
                        5   1   3   k5       % cart 1     -- cart 3
                        6   2   3   k6       % cart 2     -- cart 3
                        7   3   5   k7       % cart 3     -- right wall
                        8   3   5   k8 ];    % cart 3     -- right wall

%% =======================================================================
%%  3 -- Specify which nodes are fixed
%% =======================================================================
fixedNodes = [4 5];          % neither wall moves
allNodes   = 1:numberOfNodes;
freeNodes  = setdiff(allNodes, fixedNodes);   % everything that is not fixed

%% =======================================================================
%%  4 -- Assemble the global stiffness matrix
%% =======================================================================
localK  = [ 1 -1;
           -1  1];
globalK = zeros(numberOfNodes, numberOfNodes);
for e = 1:numberOfElements
    nodei = elementConnectivity(e,2);
    nodej = elementConnectivity(e,3);
    elementStiffness = elementConnectivity(e,4);
    globalK(nodei,nodei) = globalK(nodei,nodei) + elementStiffness*localK(1,1);
    globalK(nodei,nodej) = globalK(nodei,nodej) + elementStiffness*localK(1,2);
    globalK(nodej,nodei) = globalK(nodej,nodei) + elementStiffness*localK(2,1);
    globalK(nodej,nodej) = globalK(nodej,nodej) + elementStiffness*localK(2,2);
end

%% Define the global force vector
f = zeros(numberOfNodes,1);                   % zero at every node to start
f(1) = P1;  f(2) = P2;  f(3) = P3;            % loads on the three carts
                                              % nodes 4, 5 (walls) carry none

%% Show the assembled system, before any boundary condition
fprintf('globalK (5 x 5, both walls included):\n');
fprintf('   %9.0f %9.0f %9.0f %9.0f %9.0f\n', globalK');
fprintf('\nf'' = [%.0f %.0f %.0f %.0f %.0f]\n', f);
fprintf('\n  row sums of globalK: [%.0f %.0f %.0f %.0f %.0f]  <- all zero, because a\n', sum(globalK,2));
fprintf('  rigid-body translation of everything stretches no spring.  That is\n');
fprintf('  why globalK is SINGULAR and must be condensed before it can be solved.\n');
fprintf('  rank(globalK) = %d of %d\n\n', rank(globalK), numberOfNodes);

%% =======================================================================
%%  5 -- Condense, and solve for the free nodes
%% =======================================================================
% Striking the wall rows and columns IS the boundary condition: a prescribed
% dof carries no equilibrium equation.
freeK = globalK(freeNodes, freeNodes);
u = zeros(numberOfNodes,1);                 % zero at the fixed nodes
u(freeNodes) = freeK \ f(freeNodes);        % one linear solve

%% ---- report -------------------------------------------------------------
fprintf('freeK (3 x 3, after applying u(4) = u(5) = 0):\n');
fprintf('   %9.0f %9.0f %9.0f\n', freeK');
fprintf('  symmetric %d, positive definite %d\n\n', issymmetric(freeK), all(eig(freeK)>0));

fprintf('displacements:\n');
fprintf('   u(%d) = %+.8f m\n', [allNodes; u']);

fprintf('\nchecks:\n');
fprintf('  residual ||freeK*u - f|| = %.2e N   (the out-of-balance force)\n', ...
        norm(freeK*u(freeNodes) - f(freeNodes)));
fprintf('  potential energy 0.5*u''Ku - u''f = %.6f J\n', ...
        0.5*u(freeNodes)'*freeK*u(freeNodes) - u(freeNodes)'*f(freeNodes));
R = -(globalK(fixedNodes,:)*u);               % one reaction per fixed node
fprintf('  reactions: left wall %.1f N, right wall %.1f N, sum %.1f N\n', R, sum(R));
fprintf('  against an applied %.1f N\n', P1+P2+P3);

%% ---- now write your answer --------------------------------------------
% (a) Report freeK and the displacements.
% (b) Compare freeK with the Hessian of V your energy script printed in
%     part 3.  Are they the same matrix, and why?  