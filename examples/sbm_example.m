% Stochastic block model example (SBM)

% Seed number
rng(2021);

%% Setup graphs
% Generate G1 from SBM
[s,A1] = stochastic_block_model([10 10 10 10], [0.9 0.1 0.1 0.1; 0.1 0.9 0.1 0.1; 0.1 0.1 0.9 0.1; 0.1 0.1 0.1 0.9]);
n = length(A1);

A2 = A1;
A2(2,14) = 0;
A2(14,2) = 0;
A2(7,16) = 0;
A2(16,7) = 0;

A3 = A1;
A3(2,3) = 0;
A3(3,2) = 0;
A3(5,6) = 0;
A3(6,5) = 0;

% Convert to transition matrix
P1 = adj_to_trans(A1);
P2 = adj_to_trans(A2);
P3 = adj_to_trans(A3);

%% Define cost matrices 
% 0-1 cost
c1 = ones(n,n) - eye(n);

% Degree cost
c22 = get_degree_cost(A1, A2);
c23 = get_degree_cost(A1, A3);

%% Compute OT costs
result = zeros(5,2);

% Entropic OTC
[result(1,1), P12, stat_dist12] = entropic_otc(P1, P2, c1, 100, 1000, 50, 1000, 1);
[result(1,2), P13, stat_dist13] = entropic_otc(P1, P3, c1, 100, 1000, 50, 1000, 1);

[result(2,1), P22, stat_dist22] = entropic_otc(P1, P2, c22, 100, 1000, 50, 1000, 1);
[result(2,2), P23, stat_dist23] = entropic_otc(P1, P3, c23, 100, 1000, 50, 1000, 1);

% GOT distance
result(3,1) = got_dist(get_lap(A1), get_lap(A2));
result(3,2) = got_dist(get_lap(A1), get_lap(A3));

% GW distance
[result(4,1), ~] = fgw_dist(c1, A1, A2, ones(n, 1)/(n), ones(n, 1)/(n), 1, 1.0);
[result(4,2), ~] = fgw_dist(c1, A1, A3, ones(n, 1)/(n), ones(n, 1)/(n), 1, 1.0);

% Frobenius norm
result(5,1) = norm(get_lap(A1)-get_lap(A2), 'fro');
result(5,2) = norm(get_lap(A1)-get_lap(A3), 'fro');

% Organize results into table
rowNames = {'0-1','degree','GOT','GW', 'Frobenius'};
colNames = {'G1vsG2','G1vsG3'};
ResultTable = array2table(result,'RowNames',rowNames,'VariableNames',colNames);
disp(ResultTable);

%% Plot graphs G1, G2, G3
expid = strrep(num2str(now, '%f'), '.', '_');
format longG;
datadir = ['C:\Users\oconn\Documents\Research\OptimalJoinings\GraphOTC\Data\Examples\SBM\'];
savedir = [datadir expid '\'];
mkdir(savedir);

g1 = graph(A1);
h = plot(g1, 'NodeLabel',{});
title('G1')
saveas(gcf,[savedir 'SBM1.png'])

g1 = graph(A1);
h = plot(g1, 'NodeLabel',{});
highlight(h, [4 10],[13 17],'EdgeColor','r','LineWidth',1.5,'LineStyle','--')
title('G2')
saveas(gcf,[savedir 'SBM2.png'])

g1 = graph(A1);
h = plot(g1, 'NodeLabel',{});
highlight(h, [2 4],[3 6],'EdgeColor','r','LineWidth',1.5,'LineStyle','--')
title('G3')
saveas(gcf,[savedir 'SBM3.png'])
