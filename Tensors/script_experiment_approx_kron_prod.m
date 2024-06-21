clear;

%%% Script to compute the influence of the presampling l_tilde on the %%%
%%% quality of the approximation of the chain of kronecker products %%%

%%% Parameters
tau = 20; %Number of iterations
rank_tensor = 10; % Rank of the tensor to initialize (see the create_initialization file)
R = 10; % Rank of approximation
s = 1; % Parameter used in the construction of the examples
dims = [300,300,300]; % Dimensions of the tensor
fmt=['[' repmat(' %d ',1,numel(dims)) ']']; %to print the size in the terminal
example = 'f_1'; %Type of example
init_vec = ["STHOSVD","Sampling","Range Finder"]; %The different initializations
coeff_l_tilde = [3/4, 2:2:14]; %Coefficients of the sampling size
param_plot_vec = ["-pentagram","-d","-v","-o","-+","-^","-s","-*"]; % Shape of the different points
% Color of the different plots
marker_face_color_vect = ["#00FFFF","#4DBEEE","#0072BD","#77AC30","#7E2F8E","#EDB120","#D95319","#FF0000"];
repet = 20; %Number to repeat the experiments

%Fix the seed%
seed = 50;
rng_source = 'twister';
rng(seed,rng_source);

%%% Create the example
[M,tensor_name] = create_example(example,dims,s,rank_tensor);

size_l_tilde = size(coeff_l_tilde,2);
d = size(dims,2);

% Create arrays to get the max and min number for the y-axis
max_rel_error = 1e-16*ones([d,size(init_vec)]);
min_rel_error = 1e16*ones([d,size(init_vec)]);

absolute_error = zeros([size_l_tilde,tau,d]);


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
    max_rel_error_fixed_l = 1e-16*ones([d,size_l_tilde]);
    min_rel_error_fixed_l = 1e16*ones([d,size_l_tilde]);

    %Use Parfor to accelerate the code
    parfor i = 1:size_l_tilde

        alpha = coeff_l_tilde(i);
        l = 14*R;
        l_tilde = 4*alpha*R^2;

        %Fix the seed so that the experiment is reproductible even with
        %parfor loops (they are executed independently)
        rng(seed,rng_source);
        absolute_error_fixed_l = {};

        %Perform the iterative refinement algorithm
        for mu = 1:d
            absolute_error_fixed_l{mu} = zeros([1,tau]);
        end
        for k = 1:repet
            [~,~,error] = sublinear_iterative_refinement_rescaled_fast_distance_col(M,G_init,Matrix_List_init,tau,R,l,l_tilde);
            
            % Get the error for every matricization modes
            for mu = 1:d
                absolute_error_fixed_l{mu} = absolute_error_fixed_l{mu} + error{mu};
            end
            
        end
        
        %%% Form the relative error and print the results
        for mu = 1:d
            absolute_error_fixed_l{mu} = absolute_error_fixed_l{mu}/repet;

            max_rel_error_fixed_l(mu,i) = max(max_rel_error_fixed_l(mu,i),max(absolute_error_fixed_l{mu}));
    
            min_rel_error_fixed_l(mu,i) = min(min_rel_error_fixed_l(mu,i),min(absolute_error_fixed_l{mu}));
    
            absolute_error(i,:,mu) = absolute_error_fixed_l{mu};
        end
        
        fprintf('Experiment for l = %dR \n',alpha)
        fprintf('\n')
        fprintf('     Step N°       Error of presampling \n')
        for mu = 1:d
            fprintf('\n')
            fprintf('mu = %d\n',mu)
            for m=1:tau
                fprintf('        %d           %e    \n',m,absolute_error_fixed_l{mu}(m))
            end
        end

    end

    for mu = 1:d
        max_rel_error(mu,j) = max(max_rel_error_fixed_l(mu,:));
        min_rel_error(mu,j) = min(min_rel_error_fixed_l(mu,:));
    end
    
    %%% Plot the error
    for mu = 1:d
        for i = 1:size_l_tilde
            alpha = coeff_l_tilde(i);
            param_plot = param_plot_vec(i);
            color = marker_face_color_vect(i);
    
            figure((j-1)*d +mu);
            semilogy(1:tau,absolute_error(i,:,mu),char(param_plot),'DisplayName',sprintf('l tilde = 4*%dR^2',alpha),...
                'MarkerEdgeColor',color,'MarkerFaceColor',color,'Color',color);
            hold on;
        end
    end
    
end

%Adjust the different parameters to get a clean plot
for j = 1:size(init_vec,2)
    for mu = 1:d

        figure((j-1)*d+mu);
        set(gcf,'Position',[100,100,800,400])
        title(sprintf('%s, Error of presampling with %s initialization, for mu = %d',tensor_name,init_vec(j),mu),'Interpreter','latex')
        ylabel('$\|A-CC^{+}A\|_{F}$','Interpreter','latex')
        xlabel('N° Step');
        xlim([1 tau+1])
        ylim([10^(floor(log10(min_rel_error(mu,j)))) 10^(ceil(log10(max_rel_error(mu,j))))])
        legend('Location','eastoutside');
        legend show;
    end

end