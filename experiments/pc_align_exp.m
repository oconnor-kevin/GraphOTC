%%
% pointcloud_alignment_experiment.m
%%

clear
rng(2021);

% Setup directories
datadir = 'C:\Users\oconn\Documents\Research\OTC\GraphOTC\Data\PointCloudAlignment\';
format longG;
expid = strrep(num2str(now, '%f'), '.', '_');
expdir = [datadir expid '\'];
mkdir(expdir);

% Experiment parameters
n_iter = 5;
global n_points1 n_points2 n_clouds point_sigma dimension;
n_points1 = 10;
n_points2 = 5;
n_clouds = 4;
point_sigma = 1;
dimension = 3;
lambda_vec = [1e-2];
mean_sigma_vec = [1 2 3];

% Cross-validation parameters
do_cv = 1;
cv_iter = 10;
alpha_grid = 0:0.1:1;

% Setup data tables
values = table([], [], [], []);
values.Properties.VariableNames = {'Algorithm' 'Sigma' 'Lambda' 'Accuracy'};
edge_values = table([], [], [], []);
edge_values.Properties.VariableNames = {'Algorithm' 'Sigma' 'Lambda' 'Accuracy'};

%parpool(5);
%parfor mean_sigma=mean_sigma_vec
for mean_sigma=mean_sigma_vec
    disp(['mean sigma' num2str(mean_sigma)]);
    for lambda=lambda_vec
        disp(['lambda' num2str(lambda)]);
        if do_cv
            disp('Doing cross-validation');
            fgw_params = pc_align_cv(mean_sigma, lambda, cv_iter, alpha_grid);
            disp('Best FGW params:');
            disp(fgw_params);
        end
        
        parpool(5);
        parfor iter=1:n_iter
        %for iter=1:n_iter            
            disp(['Iteration' num2str(iter)]);          
            
            % Generate point clouds
            [points1, points2] = get_pointclouds(mean_sigma);
            
            % Construct graphs
            distmat1 = squareform(pdist(points1));
            distmat2 = squareform(pdist(points2));
            A1 = -lambda*distmat1;
            A2 = -lambda*distmat2;
            % Rescale
            A1 = (A1-min(A1,[],'all'))/(max(A1,[],'all') - min(A1,[],'all'));
            A2 = (A2-min(A2,[],'all'))/(max(A2,[],'all') - min(A2,[],'all'));
            % Threshold
            thresh = 0.1;
            A1(A1 < thresh) = 0;
            A2(A2 < thresh) = 0;

            n1 = size(A1,1);
            n2 = size(A2,1);

            % Compute transition matrices and stationary distributions
            P1 = adj_to_trans(A1);
            P2 = adj_to_trans(A2);
            %P1 = adj_to_trans(A1+1e-4);
            %P2 = adj_to_trans(A2+1e-4);
            stat_dist1 = approx_stat_dist(P1, 100)';
            stat_dist2 = approx_stat_dist(P2, 100)';
            
            % Construct cost function
            c = zeros(size(A1,1), size(A2,1));
            for i=1:size(A1,1)
                for j=1:size(A2,1)
                    c(i,j) = sum((points1(i,:) - points2(j,:)).^2);
                end
            end
            c = c ./ max(c, [], 'all');

            % Run algorithms
            [~, otc_edge_alignment, otc_alignment] = exact_otc(P1, P2, c);
            [~, fgw_alignment] = fgw_dist(c, A1, A2, stat_dist1, stat_dist2, 1, fgw_params(1));
            [~, gw_alignment] = fgw_dist(c, A1, A2, stat_dist1, stat_dist2, 1, 1);
            [otsd_alignment, ~] = computeot_lp(c', stat_dist1, stat_dist2');
            otsd_alignment = reshape(otsd_alignment, n2, n1)';
            [~, copt_alignment] = copt_dist(A1, A2, 10);
            copt_alignment = copt_alignment ./ sum(copt_alignment, 'all');
            
            %% Evaluate the alignments
            % Vertex alignment
            aligned_mass_otc = eval_alignment(otc_alignment);
            aligned_mass_fgw = eval_alignment(fgw_alignment);
            aligned_mass_gw = eval_alignment(gw_alignment);
            aligned_mass_otsd = eval_alignment(otsd_alignment);
            aligned_mass_copt = eval_alignment(copt_alignment);
            
            disp(['OTC Accuracy: ' num2str(aligned_mass_otc)]);
            disp(['FGW Accuracy: ' num2str(aligned_mass_fgw)]);
            disp(['GW Accuracy: ' num2str(aligned_mass_gw)]);            
            disp(['OT-SD Accuracy: ' num2str(aligned_mass_otsd)]);
            disp(['COPT Accuracy: ' num2str(aligned_mass_copt)]);

            values = [values; {'OTC' mean_sigma lambda aligned_mass_otc}];
            values = [values; {'FGW' mean_sigma lambda aligned_mass_fgw}];
            values = [values; {'GW' mean_sigma lambda aligned_mass_gw}];            
            values = [values; {'OT-SD' mean_sigma lambda aligned_mass_otsd}];
            values = [values; {'COPT' mean_sigma lambda aligned_mass_copt}];
            
            % Edge alignment
            % Refit FGW with edge-tuned parameters
            [~, fgw_alignment] = fgw_dist(c, A1, A2, stat_dist1, stat_dist2, 1, fgw_params(2));            
            
            edge_mass_otc = 0;
            for x_row=1:n1
                for x_col=1:n1
                    for y_row=1:n2
                        for y_col=1:n2
                            idx1 = n2*(x_row-1)+y_row;
                            idx2 = n2*(x_col-1)+y_col;
                            if fix((x_row-1)/n_points1) == fix((x_col-1)/n_points1) && fix((y_row-1)/n_points2) == fix((y_col-1)/n_points2)
                                edge_mass_otc = edge_mass_otc + otc_edge_alignment(idx1, idx2)*otc_alignment(x_row, y_row);
                            end
                        end
                    end
                end
            end

            edge_mass_fgw = 0;
            for x_row=1:n1
                for x_col=1:n1
                    for y_row=1:n2
                        for y_col=1:n2
                            idx1 = n2*(x_row-1)+y_row;
                            idx2 = n2*(x_col-1)+y_col;
                            if fix((x_row-1)/n_points1) == fix((x_col-1)/n_points1) && fix((y_row-1)/n_points2) == fix((y_col-1)/n_points2)
                                edge_mass_fgw = edge_mass_fgw + fgw_alignment(x_row, y_row)*fgw_alignment(x_col, y_col);
                            end
                        end
                    end
                end
            end
            
            
            edge_mass_gw = 0;
            for x_row=1:n1
                for x_col=1:n1
                    for y_row=1:n2
                        for y_col=1:n2
                            idx1 = n2*(x_row-1)+y_row;
                            idx2 = n2*(x_col-1)+y_col;
                            if fix((x_row-1)/n_points1) == fix((x_col-1)/n_points1) && fix((y_row-1)/n_points2) == fix((y_col-1)/n_points2)
                                edge_mass_gw = edge_mass_gw + gw_alignment(x_row, y_row)*gw_alignment(x_col, y_col);
                            end
                        end
                    end
                end
            end
            
            disp(['OTC Edge Accuracy: ' num2str(edge_mass_otc)]);
            disp(['FGW Edge Accuracy: ' num2str(edge_mass_fgw)]);
            disp(['GW Edge Accuracy: ' num2str(edge_mass_gw)]);          

            edge_values = [edge_values; {'OTC' mean_sigma lambda edge_mass_otc}];
            edge_values = [edge_values; {'FGW' mean_sigma lambda edge_mass_fgw}];
            edge_values = [edge_values; {'GW' mean_sigma lambda edge_mass_gw}];
        end
        delete(gcp('nocreate'));
    end
end

% Write accuracies
disp(values);
writetable(values, [expdir 'alignment_accuracies.txt']);
disp(edge_values);
writetable(edge_values, [expdir 'edge_alignment_accuracies.txt']);

% Vertex alignment evaluation function
function score = eval_alignment(alignment)
    global n_clouds n_points1 n_points2;
    score = 0;
    for cl = 1:n_clouds
        score = score + sum(alignment(((cl-1)*n_points1 + 1):cl*n_points1, ((cl-1)*n_points2 + 1):cl*n_points2), 'all');
    end            
end

% Point cloud sampling function
function [points1, points2] = get_pointclouds(mean_sigma)
    global dimension n_clouds n_points1 n_points2 point_sigma;
    n_points1 = 10;
    n_points2 = 5;
    n_clouds = 4;
    point_sigma = 1;
    dimension = 3;
    
    % Randomly draw means of each Gaussian
    means = normrnd(0, mean_sigma, n_clouds, dimension);

    % Draw point clouds
    points1 = zeros(n_points1*n_clouds, dimension);
    for c = 1:n_clouds
        for dim=1:dimension
            points1(((c-1)*n_points1+1):c*n_points1, dim) = normrnd(means(c,dim), point_sigma, [n_points1, 1]);
        end
    end

    points2 = zeros(n_points2*n_clouds, 2);
    for c = 1:n_clouds
        for dim=1:dimension
            points2(((c-1)*n_points2+1):c*n_points2, dim) = normrnd(means(c,dim), point_sigma, [n_points2, 1]);
        end
    end
end