%% Clear everything from memory and start from scratch
clear all; close all; clc;
%% ======================================================================= %%
%  A4, TASK 1, STEP 1 -- TRAIN THE FOUR-INPUT SURROGATE
%
%  Engineering Optimisation (M-EOP), Assignment 4, Task 1.
%
%  THIS IS L07'S l07_5_train_surrogate.m, COPIED, with the three changes the
%  assignment lists and no others:
%
%    1. load truss_data_4in.mat, taking Depth, Budget, H, Angle and X;
%    2. build the input matrix from all four: P = [Depth Budget H Angle]';
%    3. numberOfNeurons = 80, saved to a4_surrogate.mat.
%
%      inputs     4 : depth Ly, budget V, load height H, load angle theta
%      outputs  106 : the design variable for every bar
%
%  SECTION 1.3 (BONUS) is the ONE line marked below: make networkShape
%  [80 80] instead of 80.  The file name follows the shape automatically, so
%  the shallow network is not overwritten and the two can be compared.
%
%  Requires: Deep Learning Toolbox, and truss_data_4in.mat in Matlab/.
%
%  RUN AND VERIFIED in MATLAB R2026a on 2026-11-14.
%% ======================================================================= %%

%% ---- THE ONE SETTING THAT CHANGES BETWEEN SECTION 1.1 AND SECTION 1.3 --
networkShape = 80;            % section 1.1 and 1.2: one hidden layer of 80
% networkShape = [80 80];     % section 1.3 [Bonus]: two hidden layers of 80

maximumEpochs = 1000;         % section 1.3 asks for this to be RAISED if the
                              % deeper net stops on 'epochs' rather than on
                              % the validation criterion

%% ---- the data ----------------------------------------------------------
% truss_data_4in.mat holds one COMPLETED size optimisation per row:
%   Depth  n x 1    the truss depth Ly, 2 m to 6 m                       [m]
%   Budget n x 1    the volume budget V, 0.06 to 0.40                  [m^3]
%   H      n x 1    load height, as a fraction of the truss depth
%   Angle  n x 1    load angle from vertical, -75 to +75 deg            [rad]
%   X      n x 106  the DESIGN VARIABLE the optimiser returned, in [0,1]
%   fStar  n x 1    the compliance reached                                [J]
% NOTE Angle is stored in RADIANS, so a query in degrees must be converted.
load(fullfile('Matlab','truss_data_4in.mat'), 'Depth', 'Budget', 'H', 'Angle', 'X');

P = [Depth Budget H Angle]';   %   4 x n -- one COLUMN per sample
T = X';                        % 106 x n -- predict the design variable itself

% NOTHING IS SCALED BY HAND.  The design variable is bounded in [0,1] by the
% optimiser's own lower and upper bounds, so the targets arrive O(1) already.
% Predicting areas instead would mean targets of ~1e-3 m^2, small enough that
% the gradient trips min_grad below after 8 epochs and training quits.

fprintf('%d completed size optimisations\n', numel(H));
fprintf('inputs : depth, budget, load height, load angle\n');
fprintf('outputs: %d design variables\n', size(X,2));
fprintf('ranges : Ly %.2f..%.2f m, V %.3f..%.3f m^3, H %s, theta %+.1f..%+.1f deg\n', ...
        min(Depth), max(Depth), min(Budget), max(Budget), ...
        mat2str(unique(H)', 3), min(Angle)*180/pi, max(Angle)*180/pi);

%% ---- 1. CONSTRUCT ------------------------------------------------------
rng(1);                                      % Fix the starting weights and the split
numberOfNeurons = networkShape;              % 80, or [80 80] for section 1.3
net = fitnet(numberOfNeurons, 'trainscg');   % 4 inputs -> hidden -> 106 outputs

% A VECTOR gives MORE HIDDEN LAYERS -- the only change needed to go deeper.
fprintf('network: 4 -> %s -> %d, %d hidden layer(s)\n', ...
        mat2str(numberOfNeurons), size(T,1), numel(numberOfNeurons));

% Set the activation functions for all hidden layers
for k = 1:net.numLayers-1                            % Every HIDDEN layer, however many
    net.layers{k}.transferFcn = 'tansig';            % tanh in every hidden layer
end
net.layers{net.numLayers}.transferFcn = 'purelin';   % The OUTPUT layer, left linear

%% ---- 2. TRAIN ----------------------------------------------------------
net.divideFcn = 'dividerand';               % Divide data randomly into training, validation, and test sets
net.divideParam.trainRatio = 0.70;          % Fitted: the error actually minimised
net.divideParam.valRatio   = 0.15;          % Never fitted; decides WHEN TO STOP
net.divideParam.testRatio  = 0.15;          % Touched once, at the end, as the score

% Training terminates at whichever of these five stopping criteria is met
net.trainParam.max_fail = 6;            % Stop once validation error rose 6 epochs
net.trainParam.epochs   = maximumEpochs;% Hard cap on epochs
net.trainParam.min_grad = 1e-6;         % Stop once the gradient gets this flat
net.trainParam.time     = Inf;          % No wall-clock limit
net.trainParam.goal     = 0;            % Target error; 0, so it never fires

net.trainParam.showWindow      = false; % false so the script also runs in -batch
net.trainParam.showCommandLine = true;  % Also print in command line

[net, tr] = train(net, P, T);           % The whole optimisation; P and T already columns

fprintf('kept epoch %d of %d run\n', tr.best_epoch, numel(tr.perf)-1);  % net holds THIS epoch
fprintf('train MSE %.4e\n', tr.best_perf);    % Fitted to these; always the optimistic estimate
fprintf('valid MSE %.4e\n', tr.best_vperf);   % Chose the epoch, so not a fair score
fprintf('TEST  MSE %.4e\n', tr.best_tperf);   % On the design variable, in [0,1]
fprintf('stopping reason: %s\n', tr.stop);    % Why training ended

%% ---- did it converge, or just run out of patience? --------------------
% Section 1.3 item 1 turns on this distinction, so it is tested, not eyeballed.
stoppedOnEpochs = contains(lower(tr.stop), 'epoch');
if stoppedOnEpochs
    fprintf('\n  *** STOPPED ON THE EPOCH CAP (%d) -- NOT CONVERGED ***\n', maximumEpochs);
    fprintf('  Its MSE cannot fairly be compared with a run that stopped on the\n');
    fprintf('  validation criterion.  Raise maximumEpochs and train again.\n');
else
    fprintf('\n  converged: stopped on ''%s'', after %d epochs, keeping epoch %d\n', ...
            strtrim(tr.stop), numel(tr.perf)-1, tr.best_epoch);
    fprintf('  the %d epochs after the kept one are max_fail = %d validation failures\n', ...
            numel(tr.perf)-1 - tr.best_epoch, net.trainParam.max_fail);
end

%% ---- 3. SAVE -----------------------------------------------------------
% The network predicts the design variable, which needs no units and no
% rescaling, so the artefact is the network and nothing else.  The training
% record tr goes with it: section 1.3 compares the two runs' test MSE.
if isscalar(numberOfNeurons)
    surrogateFile = 'a4_surrogate.mat';            % section 1.1 / 1.2
else
    surrogateFile = 'a4_surrogate_deep.mat';       % section 1.3 [Bonus]
end
save(surrogateFile, 'net', 'tr');
fprintf('\nsaved %s\n', surrogateFile);

%% ---- Visualise the network and its performance ------------------------
figure; plotperform(tr)               % train/validation/test vs epoch
lightFigure(gcf);
exportgraphics(gcf, sprintf('a4_q1_plotperform_%s.pdf', ...
               strjoin(string(numberOfNeurons), 'x')));

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
