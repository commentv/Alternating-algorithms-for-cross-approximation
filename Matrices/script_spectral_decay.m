clear;

%%% Parameters

%examples used
example_vect = ["Slow decay", "Fast decay", "Gaussian Bump", "Cauchy","f_1"];
%colors used in the graph
marker_face_color_vect = ["#4DBEEE","#0072BD","#77AC30","#EDB120","#FF0000"];
%Sizes of the matrices of the different examples
sizes = [3000,3000,3000,2000,3000];
s = 1;
r = 10;

%%% Fix the seed
seed = 50;
rng_source = 'twister';
rng(seed,rng_source);

max_value = 1e-16;

for i = 1:size(example_vect,2)
    
    %Get the different parameters for the loop
    example = example_vect(i);
    color = marker_face_color_vect(i);
    n = sizes(i);
    
    %Create the example
    [M,matrix_name] = create_example(example,n,r,s);
    
    %Get the singular value decay
    Sigma_vect = svd(M);
    
    %Get the maximum value of the singular values for the plots
    max_value = max(max_value,Sigma_vect(1));
    
    %Plot the corresponding singular decay
    figure(1);
    semilogy(1:n,Sigma_vect,'-','DisplayName',matrix_name,...
            'MarkerEdgeColor',color,'MarkerFaceColor',color,'Color',color);
    hold on;

end

%Adjust the different parameters to get a clean plot
figure(1);
set(gcf,'Position',[100,100,800,400])
title('Spectral decay','Interpreter','latex')
ylabel('$\sigma_{i}(M)$','Interpreter','latex')
xlabel('i');
xlim([1 1500])
ylim([1e-16 10^(ceil(log10(max_value)))])
legend('Interpreter','latex')
legend('Location','eastoutside');
legend show;