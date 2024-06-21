function [M, tensor_name] = create_example(example,dims,s,rank_tensor)

    %%%
    % Inputs : A string "example" denoting the example to create, the
    % dimension of the tensor "dims", an integer "s" representing a
    % parameter of some examples and an integer "rank_tensor" representing
    % the multilinear rank of the tensor to initialize.
    %
    % Outputs : The tensor created from the type of example, as well as its 
    % name for plotting purposes.
    %%%

    switch example
        case 'Gaussian'

            d = length(dims);

            Matrix_List = {};
            G = randn(rank_tensor*ones([1,d]));

            for mu = 1:d
                Matrix_List{mu} = randn([dims(1),rank_tensor]);
            end

            [A,B] = matricization_Tucker(G,Matrix_List,1);
            Matricization_1 = A*B;

            M = tensorization(Matricization_1,1,dims);

            norm_M = tensornorm(M);

            M = M/norm_M;

            tensor_name = 'Gaussian';

           
        case 'f_1'
            
            % This type of functions are used for vectorization in matlab
            f = @(x,y,z,s) 1./(bsxfun(@plus,bsxfun(@plus,x'.^s,y.^s),reshape(z.^s,[1,1,length(z)])).^(1/s));
            M = f(1:dims(1),1:dims(2),1:dims(3),s);
            tensor_name = sprintf('$f_1(x,y,z,s)$, with s = %d',s);

         case 'Gaussian Bump'
            
            d = length(dims);

            radius = zeros(dims);

            for mu = 1:d
                n_mu = dims(mu);
                vect = linspace(-1,1,n_mu);
                Tensor = ones(dims);

                matricization_mu = vect'.*matricization(Tensor,mu);
                Tensor = tensorization(matricization_mu,mu,dims);
                radius = radius + Tensor.^2;
            end

            radius = sqrt(radius);
            radius_indicator = ceil(1-radius);
            M = exp(-1./(1-radius.^2));
            M(M==inf)=0;
            M = M.*radius_indicator;
            tensor_name = 'Gaussian Bump';

        case 'f_2'
            
            % This type of functions are used for vectorization in matlab
            f = @(x,y,z) exp(-sqrt(bsxfun(@plus,bsxfun(@plus,(x'-1).^2,(y-1).^2),reshape((z-1).^2,[1,1,length(z)]))));
            M = f(linspace(-1,1,dims(1)),linspace(-1,1,dims(2)),linspace(-1,1,dims(3)));
            tensor_name = '$f_2(x,y,z)$';

        case 'f_3'
            
            % This type of functions are used for vectorization in matlab
            f = @(x,y,z) cosh(3*bsxfun(@plus,bsxfun(@plus,x',y),reshape(z,[1,1,length(z)]))).^(-2);
            M = f(linspace(-1,1,dims(1)),linspace(-1,1,dims(2)),linspace(-1,1,dims(3)));
            tensor_name = '$f_3(x,y,z)$';

        case 'f_4'

            % This type of functions are used for vectorization in matlab
            f = @(x,y,z) 1e5/(1+1e5*(bsxfun(@plus,bsxfun(@plus,x'.^2,y.^2),reshape(z.^2,[1,1,length(z)]))));
            M = f(linspace(-1,1,dims(1)),linspace(-1,1,dims(2)),linspace(-1,1,dims(3)));
            tensor_name = '$f_4(x,y,z)$';

        case 'Exponential'

            % This type of functions are used for vectorization in matlab
            f = @(x,y,z) -exp(-1/2*(bsxfun(@plus,bsxfun(@plus,x'.^2,y.^2),reshape(z.^2,[1,1,length(z)]))));
            M = f(linspace(-1,1,dims(1)),linspace(-1,1,dims(2)),linspace(-1,1,dims(3)));
            tensor_name = 'Exponential';

        case 'Alpine'

            % This type of functions are used for vectorization in matlab
            f = @(x,y,z) bsxfun(@plus,bsxfun(@plus,abs(10*x'.*sin(10*x')+x')...
                ,abs(10*y.*sin(10*y)+y)),reshape(abs(10*z.*sin(10*z)+z),[1,1,length(z)]));
            M = f(linspace(-1,1,dims(1)),linspace(-1,1,dims(2)),linspace(-1,1,dims(3)));
            tensor_name = 'Alpine';

    end
return