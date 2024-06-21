clear;

%%% Script to compute the transition diagrams for the matrix case %%%

%%% Parameters
tau = 20; %Number of iterations
r = 10; % Rank of the matrix to initialize (see the create_initialization file)
n = 3000; % Size of the matrix (square matrix)
s = 1;
example = 'Fast decay'; %Type of example
init_vec = ["Truncation"]; %The different initializations
coeff_l = 2:2:20; %Coefficients of the sampling size
r_vect = 6:6:60; %The different ranks tested
repet = 20; %Number to repeat the experiments
iter_CUR = 3; %Number of iterations for the initial CUR approximation

%Fix the seed%
seed = 50;
rng_source = 'twister';
rng(seed,rng_source);

%%% Create the example
[M,matrix_name] = create_example(example,n,r,s);

[U,Sigma,V] = svd(M);

fprintf('--- Test matrix : %s, Size : %d x %d , rank : %d, repetition : %d ---\n',matrix_name,n,n,r,repet)

L = length(coeff_l);

%Set the y axis labels
y_axis_labels = {};
for i = 1:L
    y_axis_labels{i} = strcat(num2str(coeff_l(L+1-i)),'R');
end

%Start the loops for the different plots
for i = 1:size(init_vec,2)
    
    init = init_vec(i);

    fprintf('--- %s --- \n',init)
    fprintf('\n')

    phase_transition_error = zeros([length(coeff_l),length(r_vect)]);

    for j = 1:length(r_vect)
        
        %Get the best approximation depending on the considered rank r
        r = r_vect(j);

        U_r = U(:,1:r);
        Sigma_r_V_r_T = Sigma(1:r,1:r)*V(:,1:r)';
        
        Best_approx = U_r*Sigma_r_V_r_T;
        
        best_error = norm(M - Best_approx,'fro');

        %%% Form the initial LRA
        rng(seed,rng_source);
        [A_0,B_0] = create_initialization(init,M,U_r,Sigma_r_V_r_T,r,iter_CUR);

        fprintf('\n');
        fprintf('Experiment for r = %d \n',r);
        fprintf('\n');
        
        %Use Parfor to accelerate the code
        parfor k = 1:L

            alpha = coeff_l(k);

            
            %Fix the seed so that the experiment is reproductible even with
            %parfor loops (they are executed independently)
            relative_error = 0;
            rng(seed,rng_source);

            %Perform the iterative refinement algorithm
            for p = 1:repet
                [~,~,error] = iterative_refinement_fast_transition_diagram(M,A_0,B_0,tau,r,alpha*r);
                relative_error = relative_error + error;
            end
            
            %Print the relative error obtained
            relative_error = relative_error/(best_error*repet);

            fprintf('Experiment for l = %dr \n',alpha)
            fprintf('\n')
            fprintf('     Step N°       Relative Error \n')
            fprintf('        20           %f    \n',relative_error)
            fprintf('\n')
            
            %Recuperate the relative error (parfor loops need special care)
            phase_transition_error(k,j) = relative_error;

        end

    end
    
    %Plot the obtaine phase transition diagram for one initialization
    phase_transition_error = flip(phase_transition_error,1);

    figure(i);
    h = heatmap(r_vect,y_axis_labels,phase_transition_error,'CellLabelColor','none',...
        'ColorLimits',[1.10 1.20],'ColorbarVisible','off');
    h.NodeChildren(3).Title.Interpreter = 'latex';
    title(sprintf('%s, relative error with %s initialization',matrix_name,init_vec(i)));
    ylabel('Size of the sampling l');
    xlabel('Rank r');
    color = colormap("gray");
    colormap(flipud(color));
end