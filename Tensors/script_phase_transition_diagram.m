clear;

%%% Script to compute the transition diagrams for the tensor case %%%

%%% Parameters
tau = 20; %Number of iterations
rank_tensor = 10; % Rank of the tensor to initialize (see the create_initialization file)
s = 2; % Parameter used in the construction of the examples
dims = [300,300,300]; % Dimensions of the tensor
d = length(dims);
fmt=['[' repmat(' %d ',1,numel(dims)) ']']; %to print the size in the terminal
example = 'f_1'; %Type of example
init_vec = ["Sampling","Range Finder"]; %The different initializations
coeff_l = 2:2:20; %Coefficients of the sampling size
r_vect = 4:4:40; % Different ranks used
repet = 5; %Number to repeat the experiments

%Fix the seed%
seed = 50;
rng_source = 'twister';
rng(seed,rng_source);

%%% Create the example
[M,tensor_name] = create_example(example,dims,s,rank_tensor);

M_1 = matricization(M,1);

best_error_vect = [];

%%% Get the best low rank approximation for every ranks
for i = 1:length(r_vect)
    R = r_vect(i);
    [G_THOSVD,List_Matrix_THOSVD] = create_initialization("THOSVD",M,R);
    [A_1,B_1] = matricization_Tucker(G_THOSVD,List_Matrix_THOSVD,1);
    best_error_vect = [best_error_vect norm(M_1-A_1*B_1,"fro")];
end

fprintf(strcat('--- Test tensor : %s, Size : ', fmt ,' repetition : %d ---\n'),tensor_name,dims,repet);
fprintf('--- Seed : %d, Random Generator Source : %s ---\n',seed,rng_source)

L = length(coeff_l);
y_axis_labels = {};

% Set the y_axis labels
parfor i = 1:L
    y_axis_labels{i} = strcat(num2str(coeff_l(L+1-i)),'R');
end

for i = 1:size(init_vec,2)
    
    init = init_vec(i);

    fprintf('--- %s --- \n',init)
    fprintf('\n')

    phase_transition_error = zeros([length(coeff_l),length(r_vect)]);

    for j = 1:length(r_vect)

        R = r_vect(j);

        best_error = best_error_vect(j);

        %%% Form the initial LRA
        rng(seed,rng_source);
        [G_init,Matrix_List_init] = create_initialization(init,M,R);


        fprintf('\n');
        fprintf('Experiment for R = %d \n',R);
        fprintf('\n');
        
        %Use Parfor to accelerate the code
        parfor k = 1:L

            alpha = coeff_l(k);

            %Fix the seed so that the experiment is reproductible even with
            %parfor loops (they are executed independently)
            relative_error_array = zeros([1,repet]);
            rng(seed,rng_source);

            %Perform the iterative refinement algorithm
            for p = 1:repet
                [~,~,relative_error_array(p)] = sublinear_iterative_refinement_rescaled_fast(M,G_init,Matrix_List_init,tau,R,alpha*R^2,4*alpha*(R^(d-1)));
            end

            relative_error_array = relative_error_array / best_error;
            
            relative_error = mean(relative_error_array);

            fprintf('Experiment for l = %dr \n',alpha)

            fprintf('     Step N°       Relative Error \n')
            
            fprintf('        20           %f    \n',relative_error)
            
            fprintf('\n')

            phase_transition_error(k,j) = relative_error;

        end

    end
    
    %%% Plot the error
    phase_transition_error = flip(phase_transition_error,1);

    figure(i);
    h = heatmap(r_vect,y_axis_labels,phase_transition_error,'CellLabelColor','none',...
        'ColorLimits',[4 5],'ColorbarVisible','off');
    h.NodeChildren(3).Title.Interpreter = 'latex';
    title(sprintf('%s, relative error with %s initialization',tensor_name,init_vec(i)))
    ylabel('Size of the sampling l');
    xlabel('Rank R');
    color = colormap("gray");
    colormap(flipud(color));

end