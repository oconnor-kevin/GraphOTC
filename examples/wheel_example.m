% Wheel graph example

% Seed number
rng(2021);

%% Setup graphs
% Number of spokes
n = 15;

% Construct adjacency matrices
A1 = zeros(n+1);
A1(2:(n+1)) = 1;
A1(3:(n+1),2:n) = eye(n-1);
A1((n+1),2) = 1;
A1 = A1 + A1';

A2 = A1;
A2(2,3) = 0;
A2(3,2) = 0;

A3 = A1;
A3(1,2) = 0;
A3(2,1) = 0;

% Convert to transition matrix
P1 = adj_to_trans(A1);
P2 = adj_to_trans(A2);
P3 = adj_to_trans(A3);

%% Define cost matrices 
% 0-1 cost
c1 = ones(n+1,n+1) - eye(n+1);

% Degree cost
c22 = get_degree_cost(A1, A2);
c23 = get_degree_cost(A1, A3);

%% Compute OT costs
result = zeros(5,2);

% Exact OTC 
[result(1,1), P12, stat_dist12] = exact_otc(P1, P2, c1);
[result(1,2), P13, stat_dist13] = exact_otc(P1, P3, c1);

[result(2,1), P22, stat_dist22] = exact_otc(P1, P2, c22);
[result(2,2), P23, stat_dist23] = exact_otc(P1, P3, c23);

% GOT distance
result(3,1) = got_dist(get_lap(A1), get_lap(A2));
result(3,2) = got_dist(get_lap(A1), get_lap(A3));

% GW distance
[result(4,1), ~] = fgw_dist(c1, A1, A2, ones(n+1, 1)/(n+1), ones(n+1, 1)/(n+1), 1, 1.0);
[result(4,2), ~] = fgw_dist(c1, A1, A3, ones(n+1, 1)/(n+1), ones(n+1, 1)/(n+1), 1, 1.0);

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
datadir = ['C:\Users\oconn\Documents\Research\OptimalJoinings\GraphOTC\Data\Examples\Wheel\'];
savedir = [datadir expid '\'];
mkdir(savedir);

g1 = graph(A1);
h = plot(g1, 'NodeLabel',{});
title('G1')
saveas(gcf,[savedir 'hub_n15_G1.png'])

g1 = graph(A1);
h = plot(g1, 'NodeLabel',{});
highlight(h, [1],[2],'EdgeColor','r','LineWidth',1.5,'LineStyle','--')
title('G2')
saveas(gcf,[savedir 'hub_n15_G2.png'])

g1 = graph(A1);
h = plot(g1, 'NodeLabel',{});
highlight(h, [2],[3],'EdgeColor','r','LineWidth',1.5,'LineStyle','--')
title('G3')
saveas(gcf,[savedir 'hub_n15_G3.png'])
