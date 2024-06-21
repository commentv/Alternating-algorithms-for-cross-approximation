function [M, matrix_name] = create_example(example,n,r,s)

    %%%
    % Inputs : A string "example" denoting the example to create, the size
    % of the matrix n (square matrix), the parameter r denoting the rank of
    % some examples, and the parameter s for the f_1 example.
    %
    % Outputs : The matrix created from the type of example, as well as its 
    % name for plotting purposes.
    %%%

    switch example
        case 'Fast decay'

            Sigma = diag([ones([1,r]) 2.^(-(1:(n-r)))]);
            [U,~] = qr(randn([n,n]));
            [V,~] = qr(randn([n,n]));
            M = U*Sigma*V';
            matrix_name = 'Fast decay';

        case 'Slow decay'

            Sigma = diag([ones([1,r]) (2:(n+1-r)).^-2]);
            [U,~] = qr(randn([n,n]));
            [V,~] = qr(randn([n,n]));
            M = U*Sigma*V';
            matrix_name = 'Slow decay';

        case 'Cauchy'
            X = unifrnd(0,100,[n 1]);
            Y = -unifrnd(100,200,[n 1]);
            M = gallery('cauchy',X,Y);
            matrix_name = "Cauchy";

        case 'Shaw'
            % Use the file "shaw.m"
            M = shaw(n);
            matrix_name = 'Shaw';

        case 'Gaussian Bump'

            x = linspace(-1,1,n);
            radius = sqrt(x'.^2+x.^2);
            radius_indicator = ceil(1-radius);
            M = exp(-1./(1-radius.^2));
            M(M==inf)=0;
            M = M.*radius_indicator;
            matrix_name = 'Gaussian Bump';

        case 'Exponential'

            f = @(x,y) -exp(-1/2*(bsxfun(@plus,x'.^2,y.^2)));
            M = f(linspace(-1,1,n),linspace(-1,1,n));
            matrix_name = 'Exponential';

        case 'f_1'

            f = @(x,y,s) 1./(bsxfun(@plus,x'.^s,y.^s).^(1/s));
            M = f(1:n,1:n,s);
            matrix_name = sprintf('$f_1(x,y,s)$, with s = %d',s);

    end