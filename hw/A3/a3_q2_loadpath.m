%% Clear everything from memory and start from scratch
clear; clc; close all

%% ======================================================================= %%
%  A3, Question 2 -- check 4 of the Lecture 6 workflow: NAME THE LOAD PATH
%
%  The structure plot of a3_q2_truss_sizeOpt.m draws every bar in one
%  colour, so it shows WHERE the material went but not WHAT each member
%  does.  This script recovers the axial force of every bar of the
%  optimised baseline and redraws it coloured by sign, so the claim
%  "upper chord in tension, lower chord in compression" is read off the
%  numbers rather than asserted.
%
%  For bar e with direction cosines (l, m) and end displacements u_i, u_j,
%
%      N_e = (E A_e / L_e) * ( (u_j - u_i) . (l, m) )      [N]
%
%  i.e. axial stiffness times the elongation.  N_e > 0 is TENSION.
%  This is the same K_e^unit that supplied the sensitivity, used forwards
%  instead of backwards.
%
%  Requires: base MATLAB only (a3_q2_results.mat must exist -- run
%            a3_q2_parameter_study.m first).
%
%  RUN AND VERIFIED in MATLAB R2026a on 2026-11-14.
%% ======================================================================= %%

S = load('a3_q2_results.mat');          % written by a3_q2_parameter_study.m
E  = 2.1e11;                            % Young's modulus [Pa]
P  = 1e5;                               % applied load    [N]

r  = S.results(1).out;                  % run 1 = the baseline, 7 x 5 'full'
nodalCoordinates    = r.nodalCoordinates;
elementConnectivity = r.elementConnectivity;
aStar               = r.aStar;
ne                  = r.numberOfElements;
nNode               = size(nodalCoordinates,1);

%% ---- 1. rebuild the analysis of the optimised design -------------------
numberOfDOFs = 2*nNode;
globalK      = zeros(numberOfDOFs);
elementLength = zeros(ne,1);
dirCos        = zeros(ne,2);
dofsOf        = zeros(ne,4);

for e = 1:ne
    i = elementConnectivity(e,2);   j = elementConnectivity(e,3);
    dx = nodalCoordinates(j,2) - nodalCoordinates(i,2);
    dy = nodalCoordinates(j,3) - nodalCoordinates(i,3);
    Le = hypot(dx, dy);
    l  = dx/Le;   m = dy/Le;

    elementLength(e) = Le;
    dirCos(e,:)      = [l m];
    dofs             = [2*i-1, 2*i, 2*j-1, 2*j];
    dofsOf(e,:)      = dofs;

    unitK = (E/Le) * [ l*l  l*m -l*l -l*m
                       l*m  m*m -l*m -m*m
                      -l*l -l*m  l*l  l*m
                      -l*m -m*m  l*m  m*m ];
    globalK(dofs,dofs) = globalK(dofs,dofs) + aStar(e)*unitK;
end

forceVector = zeros(numberOfDOFs,1);
forceVector(2*r.loadNode) = -P;

fixedDOFs = sort([2*r.supportNodes-1; 2*r.supportNodes])';
freeDOFs  = setdiff(1:numberOfDOFs, fixedDOFs);

u = zeros(numberOfDOFs,1);
u(freeDOFs) = globalK(freeDOFs,freeDOFs) \ forceVector(freeDOFs);

compliance = forceVector' * u;

%% ---- 2. axial force in every bar ---------------------------------------
axialForce = zeros(ne,1);
for e = 1:ne
    ue = u(dofsOf(e,:));
    elongation    = (ue(3) - ue(1))*dirCos(e,1) + (ue(4) - ue(2))*dirCos(e,2);
    axialForce(e) = E*aStar(e)/elementLength(e) * elongation;   % + = tension
end

%% ---- 3. keep only the structural bars ----------------------------------
keep = aStar > 0.01*max(aStar);          % the same 1% rule as the histogram

%% ---- 4. report ---------------------------------------------------------
fprintf('================ CHECK 4 -- THE LOAD PATH ================\n');
fprintf('baseline 7 x 5 ''full'': %d bars, %d above the 1%% line\n', ne, sum(keep));
fprintf('compliance recovered from the saved design: %.6e J\n', compliance);
fprintf('  (a3_q2_parameter_study.m reported          %.6e J)\n\n', r.fStar);

fprintf('of the %d structural bars:\n', sum(keep));
fprintf('  in TENSION     (N > 0): %d\n', sum(keep & axialForce > 0));
fprintf('  in COMPRESSION (N < 0): %d\n', sum(keep & axialForce < 0));
fprintf('  largest tension      : %+.4e N\n', max(axialForce(keep)));
fprintf('  largest compression  : %+.4e N\n', min(axialForce(keep)));

%% the support reactions -- the couple that resists the overall moment
reactionFull = globalK*u - forceVector;
fprintf('\nTHE SUPPORT REACTIONS (the couple)\n');
for k = 1:numel(r.supportNodes)
    nd = r.supportNodes(k);
    fprintf('  support at y = %.2f m : Rx = %+.4e N,  Ry = %+.4e N\n', ...
            nodalCoordinates(nd,3), reactionFull(2*nd-1), reactionFull(2*nd));
end
fprintf('  sum of Ry = %+.4e N  (must balance the applied %.4e N)\n', ...
        sum(reactionFull(2*r.supportNodes)), P);
fprintf('  sum of Rx = %+.4e N  (must vanish: no applied horizontal load)\n', ...
        sum(reactionFull(2*r.supportNodes-1)));
leverArm     = abs(diff(nodalCoordinates(r.supportNodes,3)));
coupleMoment = abs(reactionFull(2*r.supportNodes(1)-1)) * leverArm;
fprintf('  the horizontal pair is a COUPLE: %.4e N x %.2f m = %.4e N m\n', ...
        abs(reactionFull(2*r.supportNodes(1)-1)), leverArm, coupleMoment);
fprintf('  the applied load makes P x Lx          = %.4e N m  (they must agree)\n', P*S.Lx);

%% what the bars touching each support do
fprintf('\nTHE BARS MEETING EACH SUPPORT\n');
for k = 1:numel(r.supportNodes)
    nd  = r.supportNodes(k);
    at  = keep & (elementConnectivity(:,2) == nd | elementConnectivity(:,3) == nd);
    fprintf('  support at y = %.2f m : %d bars, %d tension, %d compression,\n', ...
            nodalCoordinates(nd,3), sum(at), ...
            sum(at & axialForce > 0), sum(at & axialForce < 0));
    fprintf('      thickest of them carries %+.4e N\n', ...
            axialForce(find(at & aStar == max(aStar(at)), 1)));
end

%% antisymmetry: the geometry is mirror-symmetric about y = Ly/2 and the load
%% sits on that line, so mirroring maps the structure to itself but flips the
%% load.  Every bar therefore has a mirror twin carrying exactly MINUS its force.
mirrorOf = zeros(ne,1);
coord    = nodalCoordinates(:,2:3);
mirrored = [coord(:,1), S.Ly - coord(:,2)];
nodeTwin = zeros(nNode,1);
for nd = 1:nNode
    [d, w] = min(hypot(coord(:,1)-mirrored(nd,1), coord(:,2)-mirrored(nd,2)));
    if d < 1e-9, nodeTwin(nd) = w; end
end
barKey = @(a,b) min(a,b)*1e6 + max(a,b);
keyList = arrayfun(@(e) barKey(elementConnectivity(e,2), elementConnectivity(e,3)), (1:ne)');
for e = 1:ne
    tw = barKey(nodeTwin(elementConnectivity(e,2)), nodeTwin(elementConnectivity(e,3)));
    w  = find(keyList == tw, 1);
    if ~isempty(w), mirrorOf(e) = w; end
end
hasTwin   = mirrorOf > 0;
antiError = max(abs(axialForce(hasTwin) + axialForce(mirrorOf(hasTwin)))) / ...
            max(abs(axialForce));
fprintf('\nMIRROR ANTISYMMETRY\n');
fprintf('  every bar has a twin reflected in y = Ly/2: %d of %d\n', sum(hasTwin), ne);
fprintf('  max |N_e + N_twin| / max|N| = %.1e\n', antiError);
fprintf('  (not machine precision: the OPTIMISED DESIGN is itself only symmetric\n');
fprintf('   to the optimality tolerance, firstorderopt = %.2e)\n', r.firstorderopt);

%% ---- 5. the figure -----------------------------------------------------
figure('Color','w','Position',[100 100 760 520]); hold on
aPlotRef = S.V0/S.Lx;
for e = find(keep)'
    i = elementConnectivity(e,2);  j = elementConnectivity(e,3);
    lw = max(0.5, 7*aStar(e)/max(aStar));
    if axialForce(e) > 0
        col = [0.80 0.15 0.15];      % tension  -- red
    else
        col = [0.10 0.30 0.70];      % compression -- blue
    end
    plot(nodalCoordinates([i j],2), nodalCoordinates([i j],3), '-', ...
         'Color', col, 'LineWidth', lw);
end
plot(nodalCoordinates(r.supportNodes,2), nodalCoordinates(r.supportNodes,3), ...
     '^', 'MarkerSize', 9, 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor','none');
quiver(nodalCoordinates(r.loadNode,2), nodalCoordinates(r.loadNode,3), ...
       0, -0.35*S.Ly, 0, 'LineWidth', 2, 'Color', [0.1 0.1 0.1], 'MaxHeadSize', 3);
axis equal; xlim([-0.3 S.Lx+0.5]); ylim([-0.3 S.Ly+0.3]);
xlabel('x [m]'); ylabel('y [m]');
title(sprintf('Baseline load path: red = tension, blue = compression (%d bars)', sum(keep)));
box on; hold off
lightFigure(gcf);
% PNG, not vector PDF: the report embeds raster figures, so 'ContentType'
% is dropped (it is a vector-only option) and the raster density set instead.
exportgraphics(gcf, 'a3_q2_loadpath.png', 'Resolution', 300);

fprintf('\nfigure written: a3_q2_loadpath.png\n');

%% ---- 6. assertions -- this script's test suite -------------------------
assert(abs(compliance - r.fStar)/r.fStar < 1e-10, ...
       'recovered compliance disagrees with the saved f*');
assert(all(isfinite(axialForce)), 'non-finite axial force');
% global equilibrium: the vertical reactions must carry the applied load
reactions = globalK(fixedDOFs,:)*u - forceVector(fixedDOFs);
assert(abs(sum(reactions(2:2:end)) - P) < 1e-6*P, ...
       'vertical reactions do not balance the applied load');
assert(antiError < 1e-3, 'axial forces are not mirror-antisymmetric');
assert(abs(coupleMoment - P*S.Lx) < 1e-4*P*S.Lx, ...
       'the support couple does not match the applied moment');
fprintf('all assertions passed\n');

%% ================================================================ %%
%  Local functions
%% ================================================================ %%
function lightFigure(fig)
    try
        theme(fig, 'light');
    catch
        set(fig, 'Color', 'w');
        set(findall(fig, 'Type', 'axes'), 'Color', 'w');
    end
    drawnow;
end
