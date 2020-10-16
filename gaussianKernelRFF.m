function F = gaussianKernelRFF(D,gamma,s)
    %% Original Random Fourier Feature Embedding for Gaussian Kernel
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
    % output:
    %
    %  * F : A matrix with n rows (data embeddings) and s columns (random Fourier features)
    
   [n,d] = size(D);
    g = sqrt(gamma/(2*pi^2))*randn(d,s);
    b = 2*pi*rand(1,s);
    F = sqrt(2)*cos(2*pi*D*g + b) / sqrt(s);
end 