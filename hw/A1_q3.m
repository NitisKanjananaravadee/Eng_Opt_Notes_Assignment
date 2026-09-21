% Assignment 1, Question 3 - Rao 3.101
% Grocery store buys 5 vegetables from 4 farms for one month.
% Variables x(f,v) = tonnes of vegetable v bought from farm f, so 20 of them.
% Requires: Optimization Toolbox (linprog).

clear; clc

% price [$/t], one row per farm, one column per vegetable
% columns: potato tomato okra eggplant spinach
P = [200 600 1600 800 1200      % farm 1
     300 550 1400 850 1100      % farm 2
     250 650 1500 700 1000      % farm 3
     150 500 1700 900 1300];    % farm 4

supply = [180; 200; 100; 120];  % max total tonnes each farm can deliver
demand = [100; 60; 20; 80; 40]; % min tonnes of each vegetable the store needs

[nf, nv] = size(P);
n = nf*nv;

% x is ordered farm by farm: x11 x12 ... x15 x21 ... x45
c = reshape(P', n, 1);

% one row per farm: everything bought from that farm <= its supply
A_supply = kron(eye(nf), ones(1, nv));

% one row per vegetable: bought from all farms >= demand, so flip the sign
A_demand = -repmat(eye(nv), 1, nf);

A = [A_supply; A_demand];
b = [supply; -demand];
lb = zeros(n, 1);

opts = optimoptions('linprog', 'Display', 'none');
[x, fval] = linprog(c, A, b, [], [], lb, [], opts);

X = reshape(x, nv, nf)';   % back into the shape of the price table

fprintf('buying plan [t], rows = farms, columns = potato tomato okra eggplant spinach\n');
disp(round(X, 4));
fprintf('total cost = $%.0f\n\n', fval);

fprintf('farm      bought   supply\n');
for f = 1:nf
    fprintf('farm %d   %6.1f   %6.1f\n', f, sum(X(f,:)), supply(f));
end

fprintf('\nvegetable totals [t]: ');
fprintf('%.1f ', sum(X, 1));
fprintf('\nrequired            : ');
fprintf('%.1f ', demand);
fprintf('\n');
