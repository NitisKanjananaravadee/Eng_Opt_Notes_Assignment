%% A1 Q3 SKELETON -- buying five vegetables from four farms (Rao 3.101)
%
% Engineering Optimisation (M-EOP), Assignment 1.
%
% Most of this script is written for you.  There are FOUR short gaps, each
% marked  <<< TODO n >>>, and each is a single line.  The comment just above
% every gap says exactly what is wanted.
%
% TODO 2 is already DONE FOR YOU as a worked example -- read it, then TODO 3 is
% the same idea with rows and columns swapped.
%
% If you get stuck, run the script anyway: it will fail at the first unfilled
% gap, which tells you where you are.
%
% Requires: Optimization Toolbox (linprog).

clear; clc

%% ---- the data (given -- do not change) ---------------------------------
%        potato tomato  okra  eggpl spinach
price = [  200    600   1600    800   1200      % farm 1     [$ per tonne]
           300    550   1400    850   1100      % farm 2
           250    650   1500    700   1000      % farm 3
           150    500   1700    900   1300];    % farm 4
supply = [180; 200; 100; 120];       % max tonnes per farm, ALL TYPES COMBINED
demand = [100;  60;  20;  80;  40];  % min tonnes of each vegetable required
veg    = ["potato" "tomato" "okra" "eggplant" "spinach"];

[nF, nV] = size(price);              % nF = 4 farms, nV = 5 vegetables

%% =======================================================================
%%  Build the LP and solve it
%% =======================================================================
% x will hold the 4x5 grid of purchases, flattened into one column of 20.
% MATLAB flattens column by column, so price(:) and x line up automatically.

% <<< TODO 1 >>>  The cost vector.  One line, using price(:).
% Remember x is flattened the same way, so the two line up automatically.
% c = ...

% --- farm capacity rows:  for each farm f,  sum over v of x(f,v) <= supply(f)
% The mask trick: make a 4x5 matrix of zeros, put 1s where you want to add up,
% then flatten it with M(:)' so it lines up with x.
Afarm = zeros(nF, nF*nV);
for f = 1:nF
    M = zeros(nF, nV);
    M(f, :) = 1;               % <<< TODO 2 -- DONE FOR YOU >>>
    %                            1s across row f: "everything bought FROM farm f"
    Afarm(f, :) = M(:)';       % flatten, so it lines up with x
end

% --- store requirement rows: for each vegetable v, sum over f of x(f,v) >= demand(v)
Aveg = zeros(nV, nF*nV);
for v = 1:nV
    M = zeros(nF, nV);
    % <<< TODO 3 >>>  Exactly like TODO 2 above, but down a COLUMN instead of
    % across a row: put 1s in column v of M, so this row adds up everything
    % bought OF vegetable v.  One line.

    Aveg(v, :) = M(:)';
end

% linprog only accepts rows of the form  (row) * x <= (value).  The capacity
% rows are already <=.  The demand rows are >=, so multiply both sides by -1,
% which flips the inequality.
A = [Afarm; -Aveg];
b = [supply; -demand];

% <<< TODO 4a >>>  Lower bounds: you cannot buy a negative quantity.  One line.
% There are nF*nV variables in total.
% lb = ...

opts = optimoptions('linprog', 'Display', 'none');

% <<< TODO 4b >>>  Solve it.  Ask for the solution and the objective value.
% The call is  linprog(c, A, b, Aeq, beq, lb, ub, opts)  -- there are no
% equality constraints and no upper bounds here, so pass [] for Aeq, beq and ub.
% [x, cost] = linprog( ... );

X = reshape(x, nF, nV);        % unflatten the answer back into a 4x5 grid

%% =======================================================================
%%  Reporting -- written for you.  Check it against your hand reasoning.
%% =======================================================================
fprintf('\n  buying plan (tonnes)\n');
fprintf('  farm ');  fprintf('%10s', veg);  fprintf('%10s%7s\n', 'total', 'cap');
for f = 1:nF
    fprintf('   %d  ', f);  fprintf('%10.1f', X(f,:));
    tag = "";
    if abs(sum(X(f,:)) - supply(f)) < 1e-6, tag = "  ACTIVE"; end
    fprintf('%10.1f%7.0f%s\n', sum(X(f,:)), supply(f), tag);
end
fprintf('  need');  fprintf('%10.0f', demand);  fprintf('\n');
fprintf('\n  minimum cost = $%.0f\n', cost);

fprintf('\n  vegetables bought from more than one farm:\n');
for v = 1:nV
    farms = find(X(:,v) > 1e-6);
    if numel(farms) > 1
        fprintf('     %-9s from farms %s\n', veg(v), mat2str(farms'));
    end
end

%% ---- now write your answer --------------------------------x------------
% Using the output above, answer in your report:
%   * Which farms end up at capacity (marked ACTIVE)?
%   * Which vegetables are split between two farms, and WHY?  Compare each
%     farm's price for that vegetable against what else that farm is supplying.
