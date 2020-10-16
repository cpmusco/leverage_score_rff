function F = gaussianKernelMRFF(D,gamma,s,f)
if nargin<4
  f = 3;
end
    %% Modifier Random Fourier Feature Embedding for Gaussian Kernel
    % Outputs an embedding f_i for each row d_i of the data matrix D such 
    % that the the inner product between embeddings approximates the Gaussian 
    % kernel inner product. I.e. f_i^Tf_j ~= exp(-gamma*||d_i - d_j||^2).
    %
    % usage : 
    %
    % input:
    %
    %  * D : A matrix with n rows (data points) and d columns (features)
    %  
    %  * gamma : kernel variance parameter. The kernel function approximated
    %  is k(d_i,d_j) = exp(-gamma*||d_i - d_j||^2.
    %
    %  * s : how many random features to use for each f_i
    %
    %  * f : width parameter. Random Fourier features samples will be drawn
    %  from the range [-f*sqrt(gamma/(2*pi^2)),f*sqrt(gamma/(2*pi^2))]. 
    %
    % output:
    %
    %  * F : A matrix with n rows (data embeddings) and s columns (random Fourier features)
    
    [n,d] = size(D);
    d = 2;
    sig2 = gamma/(2*pi^2);
    U = zeros(d,s);
    for j = 1:s
        done = 0;
        while(~done)
            eta = f*sqrt(sig2)*(rand(d,1)*2 - 1);
            if(norm(eta) <= f*sqrt(sig2))
                U(:,j) = eta;
                done = 1;
            end
        end
    end
    w = (2*pi*sig2)^(-d/2)*exp(-sum(U.^2)/(2*sig2))*pi*(f*sqrt(sig2))^2;  
    b = 2*pi*rand(1,s);
    F = sqrt(2)*cos(2*pi*D*U + b).*sqrt(w) / sqrt(s);
end