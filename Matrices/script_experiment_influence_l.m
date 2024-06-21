clear;

%%% Script to compute the influence of the sampling l on the absolute %%%
%%% error for the matrix case %%%

%%% Parameters
tau = 20; %Number of iterations
r = 10; % Rank of the matrix to initialize (see the create_initialization file)
n = 3000; % Size of the matrix (square matrix)
s = 2; % Parameter used in the construction of the examples
example = 'f_1'; %Type of example
init_vec = ["Sampling"]; %The different initializations
coeff_l = 20:20:140; %Coefficients of the sampling size
param_plot_vec = ["-d","-v","-o","-+","-^","-s","-*"]; % Shape of the different points
% Color of the different plots
marker_face_color_vect = ["#4DBEEE","#0072BD","#77AC30","#7E2F8E","#EDB120","#D95319","#FF0000"];
repet = 50; %Number to repeat the experiments
iter_CUR = 3; %Number of iterations for the initial CUR approximation

size_l = size(coeff_l,2);
abs_error = zeros([size_l,tau+1]);

%Fix the seed%
seed = 50;
rng_source = 'twister';
rng(seed,rng_source);

%%% Create the example
[M, matrix_name] = create_example(example,n,r,s);

%%% Get the best low rank approximation
[U,Sigma,V] = svd(M);
U_r = U(:,1:r);
V_r = V(:,1:r);

Sigma_r_M = Sigma(r,r);
Sigma_r_plus_1_M = Sigma(r+1,r+1);

Best_approx = U_r*Sigma(1:r,1:r)*V_r';

best_error = norm(M - Best_approx,'fro');

% Create arrays to get the max and min number for the y-axis
max_abs_error = best_error*ones([1,size(init_vec)]);
min_abs_error = best_error*ones([1,size(init_vec)]);


fprintf('--- Test matrix : %s, Size : %d x %d , rank : %d, repetition : %d ---\n',matrix_name,n,n,r,repet)
fprintf('--- Sigma_{r+1}(M)/Sigma_{r}(M) = %f,   theta = %f --- \n',Sigma_r_plus_1_M/Sigma_r_M,best_error/Sigma_r_plus_1_M)

%Start the loops for the different plots
for j = 1:size(init_vec,2)
    
    init = init_vec(j);
        
    %%% Form the initial LRA
    rng(seed,rng_source);
    [A_0,B_0] = create_initialization(init,M,U_r,Sigma(1:r,1:r)*V_r',r,iter_CUR);

    fprintf('--- %s --- \n',init)
    fprintf('\n')

    % Create intermediate arrays for the min and max of the y-axis (for
    % parfor loops, need to take extra care)
    max_abs_error_fixed_l = best_error*ones([1,size_l]);
    min_abs_error_fixed_l = best_error*ones([1,size_l]);
    
    %Use Parfor to accelerate the code
    parfor i = 1:size_l
        
        alpha = coeff_l(i);

        %Fix the seed so that the experiment is reproductible even with
        %parfor loops (they are executed independently)
        rng(seed,rng_source);
        abs_error_fixed_l = zeros([1,tau+1]);
        
        %Perform the iterative refinement algorithm
        for k = 1:repet
            [A,B,error] = iterative_refinement_fast_influence_l(M,A_0,B_0,tau,r,alpha*r);
            abs_error_fixed_l = abs_error_fixed_l + error;
        
        end
        
        %Print the relative error obtained
        abs_error_fixed_l = abs_error_fixed_l/(repet);
        max_abs_error_fixed_l(i) = max(max_abs_error_fixed_l(i),max(abs_error_fixed_l));
        min_abs_error_fixed_l(i) = min(min_abs_error_fixed_l(i),min(abs_error_fixed_l));
        
        fprintf('Experiment for l = %dr \n',alpha)
        fprintf('\n')
        fprintf('     Step N°       Relative Error \n')
        for k=0:tau
            fprintf('        %d           %e    \n',k,abs_error_fixed_l(k+1))
        end
        fprintf('\n')

        abs_error(i,:) = abs_error_fixed_l;
    
    end 
    
    %Update the min and max y-axis
    max_abs_error(j) = max(max_abs_error(j),max(max_abs_error_fixed_l));
    min_abs_error(j) = min(min_abs_error(j),min(min_abs_error_fixed_l));

    %%% Plot the error
    for i = 1:size_l

        alpha = coeff_l(i);
        param_plot = param_plot_vec(i);
        color = marker_face_color_vect(i);
    
        figure(j);
        semilogy(0:tau,abs_error(i,:),char(param_plot),'DisplayName',sprintf('l = %dr',alpha),...
            'MarkerEdgeColor',color,'MarkerFaceColor',color,'Color',color);
        hold on;

    end

end

%Adjust the different parameters to get a clean plot
for j = 1:size(init_vec,2)

    figure(j);
    set(gcf,'Position',[100,100,800,400])
    semilogy(0:tau,best_error*ones([1,tau+1]),'--k','DisplayName','Truncated SVD');
    title(sprintf('%s, absolute error with %s initialization',matrix_name,init_vec(j)),'Interpreter','latex')
    ylabel('$\|M-AB\|_{F}$','Interpreter','latex')
    xlabel('N° Step');
    xlim([0 tau+1])
    ylim([10^(floor(log10(min_abs_error(j)))) 10^(ceil(log10(max_abs_error(j))))])
    legend('Location','eastoutside');
    legend show;

end
