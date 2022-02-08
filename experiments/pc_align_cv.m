%%
% pc_align_cv.m
%%

function fgw_params = pc_align_cv(mean_sigma, lambda, tm_scaling, cv_iter, alpha_grid)
    dimension = 3;
    n_points1 = 10;
    n_points2 = 5;
    n_clouds = 4;
    point_sigma = 1;

    % Setup variables
    fgw_acc = zeros(length(alpha_grid), cv_iter);
    fgw_edge_acc = zeros(length(alpha_grid), cv_iter);

    parpool(5);    
    parfor iter=1:cv_iter
    %for iter=1:cv_iter
        disp(['CV Iteration' num2str(iter)]);
        
        % Iteration accuracy
        iter_fgw_acc = zeros(length(alpha_grid), 1);
        iter_fgw_edge_acc = zeros(length(alpha_grid), 1);
        
        % Generate point clouds
        means = normrnd(0, mean_sigma, [n_clouds, dimension]);

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

        % Construct graphs
        distmat1 = squareform(pdist(points1));
        distmat2 = squareform(pdist(points2));
        A1 = -distmat1;
        A2 = -distmat2;
        % Rescale
        A1 = (A1-min(A1,[],'all'))/(max(A1,[],'all') - min(A1,[],'all'));
        A2 = (A2-min(A2,[],'all'))/(max(A2,[],'all') - min(A2,[],'all'));
        % Threshold
        thresh = 0.1;
        A1(A1 < thresh) = 0;
        A2(A2 < thresh) = 0;

        n1 = size(A1,1);
        n2 = size(A2,1);
        
        % Construct transition matrices and stationary distributions
        if strcmp(tm_scaling, 'linear')
            P1 = A1./sum(A1,2);
            P2 = A2./sum(A2,2);
        elseif strcmp(tm_scaling, 'exponential')
            P1 = adj_to_trans(lambda*A1);
            P2 = adj_to_trans(lambda*A2);
        end

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

        for alpha_idx = 1:length(alpha_grid)
            alpha = alpha_grid(alpha_idx);
            
            % Run FGW and evaluate
            [~, fgw_alignment] = fgw_dist(c, A1, A2, stat_dist1, stat_dist2, 1, alpha);

            aligned_mass_fgw = 0;
            for cl = 1:n_clouds
                aligned_mass_fgw = aligned_mass_fgw + sum(fgw_alignment(((cl-1)*n_points1 + 1):cl*n_points1, ((cl-1)*n_points2 + 1):cl*n_points2), 'all');
            end            
            iter_fgw_acc(alpha_idx) = aligned_mass_fgw;
            disp(['FGW Accuracy: ' num2str(aligned_mass_fgw)]);

            edge_mass_fgw = 0;
            for x_row=1:n1
                for x_col=1:n1
                    for y_row=1:n2
                        for y_col=1:n2
                            if fix((x_row-1)/n_points1) == fix((x_col-1)/n_points1) && fix((y_row-1)/n_points2) == fix((y_col-1)/n_points2)
                                edge_mass_fgw = edge_mass_fgw + fgw_alignment(x_row, y_row)*fgw_alignment(x_col, y_col);
                            end
                        end
                    end
                end
            end
            iter_fgw_edge_acc(alpha_idx) = edge_mass_fgw;
            disp(['FGW Edge Accuracy: ' num2str(edge_mass_fgw)]);
        end
        
        fgw_acc(:,iter) = iter_fgw_acc;
        fgw_edge_acc(:,iter) = iter_fgw_edge_acc;
    end
    delete(gcp('nocreate'));
    
    fgw_acc = mean(fgw_acc, 2);
    disp('FGW mean accuracy');
    disp(fgw_acc);
    
    fgw_edge_acc = mean(fgw_edge_acc, 2);
    disp('FGW mean edge accuracy');
    disp(fgw_edge_acc);
  
    alpha_idx1 = find(fgw_acc == max(max(fgw_acc)));
    alpha_idx2 = find(fgw_edge_acc == max(max(fgw_edge_acc)));
    fgw_params = [alpha_grid(alpha_idx1(1)); alpha_grid(alpha_idx2(1))];
end