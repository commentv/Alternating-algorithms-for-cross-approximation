clear;

%%% Script to compute the influence of the sampling l on the absolute %%%
%%% error for the tensor case %%%

%%% Parameters
tau = 20; %Number of iterations
rank_tensor = 10; % Rank of the tensor to initialize (see the create_initialization file)
R = 10; % Rank of the tensor to approximate
s = 1; % Parameter used in the construction of the examples
dims = [300,300,300]; % Dimension of the tensors
fmt=['[' repmat(' %d ',1,numel(dims)) ']']; %to print the size in the terminal
example = 'f_1'; % %Type of example
init_vec = ["Sampling"]; %The different initializations
coeff_l = 2:2:14; %Coefficients of the sampling size
param_plot_vec = ["-d","-v","-o","-+","-^","-s","-*"]; % Shape of the different points
% Color of the different plots
marker_face_color_vect = ["#4DBEEE","#0072BD","#77AC30","#7E2F8E","#EDB120","#D95319","#FF0000"];
repet = 20; %Number to repeat the experiments

%Fix the seed%
seed = 50;
rng_source = 'twister';
rng(seed,rng_source);

%%% Create the example
[M,tensor_name] = create_example(example,dims,s,rank_tensor);

%%% Get the THOSVD low rank approximation

[G_THOSVD,List_Matrix_THOSVD] = create_initialization("THOSVD",M,R);
[A_THOSVD,B_THOSVD] = matricization_Tucker(G_THOSVD,List_Matrix_THOSVD,1);
error = norm(matricization(M,1)-A_THOSVD*B_THOSVD,'fro');

error_THOSVD = error*ones([1,tau+1]);

% Create arrays to get the max and min number for the y-axis
max_rel_error = error*ones([1,size(init_vec)]);
min_rel_error = error*ones([1,size(init_vec)]);

size_l = size(coeff_l,2);


absolute_error = zeros([size_l,tau+1]);


fprintf(strcat('--- Test tensor : %s, Size : ', fmt ,' , rank : %d, repetition : %d ---\n'),tensor_name,dims,R,repet);
for j = 1:size(init_vec,2)
    
    init = init_vec(j);
        
    %%% Form the initial LRA
    rng(seed,rng_source);
    [G_init,Matrix_List_init] = create_initialization(init,M,R);

    fprintf('--- %s --- \n',init)
    fprintf('\n')
    
    % Create intermediate arrays for the min and max of the y-axis (for
    % parfor loops, need to take extra care)
    max_rel_error_fixed_l = error*ones([1,size_l]);
    min_rel_error_fixed_l = error*ones([1,size_l]);
    
    %Use Parfor to accelerate the code
    parfor i = 1:size_l

        alpha = coeff_l(i);
        l = alpha*R;
        l_tilde = 4*alpha*R^2;

        %Fix the seed so that the experiment is reproductible even with
        %parfor loops (they are executed independently)
        rng(seed,rng_source);
        absolute_error_fixed_l = zeros([1,tau+1]);

        %Perform the iterative refinement algorithm
        for k = 1:repet
            [~,~,error_algo] = sublinear_iterative_refinement_rescaled_fast_infl_l(M,G_init,Matrix_List_init,tau,R,l,l_tilde);

            absolute_error_fixed_l = absolute_error_fixed_l + error_algo;
            
        end
        
        %%% Form the relative error and print the results
        absolute_error_fixed_l = absolute_error_fixed_l/repet;

        max_rel_error_fixed_l(i) = max(max_rel_error_fixed_l(i),max(absolute_error_fixed_l));

        min_rel_error_fixed_l(i) = min(min_rel_error_fixed_l(i),min(absolute_error_fixed_l));

        absolute_error(i,:) = absolute_error_fixed_l;
        
        fprintf('Experiment for l = %dR \n',alpha)
        fprintf('\n')
        fprintf('---Sublinear Algorithm ---\n')
        fprintf('     Step N°       Absolute Error \n')
        for m=0:tau
            fprintf('        %d           %e    \n',m,absolute_error_fixed_l(m+1))
        end

    end
    
    %Update the min and max y-axis
    max_rel_error(j) = max(max_rel_error_fixed_l);
    min_rel_error(j) = min(min_rel_error_fixed_l);
    
    %%% Plot the error
    for i = 1:size_l
        alpha = coeff_l(i);
        param_plot = param_plot_vec(i);
        color = marker_face_color_vect(i);

        figure(j);
        semilogy(0:tau,absolute_error(i,:),char(param_plot),'DisplayName',sprintf('l = %dR',alpha),...
            'MarkerEdgeColor',color,'MarkerFaceColor',color,'Color',color);
        hold on;
    end
    
end

%Adjust the different parameters to get a clean plot
for j = 1:size(init_vec,2)

    figure(j);
    set(gcf,'Position',[100,100,800,400])
    semilogy(0:tau,error_THOSVD,'--k','DisplayName','THOSVD');
    title(sprintf('%s, absolute error with %s initialization',tensor_name,init_vec(j)),'Interpreter','latex')
    ylabel('$\|M-\hat{M}\|_{F}$','Interpreter','latex')
    xlabel('N° Step');
    xlim([0 tau+1])
    ylim([10^(floor(log10(min_rel_error(j)))) 10^(ceil(log10(max_rel_error(j))))])
    legend('Location','eastoutside');
    legend show;

end