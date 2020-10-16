function F = cauchyKernelRFF(D,gamma,s)
    %% Original Random Fourier Feature Embedding for Cauchy Kernel
    % Outputs an embedding f_i for each row d_i of the data matrix D such 
    % that the the inner product between embeddings approximates the Cauchy 
    % kernel inner product. I.e. f_i^Tf_j ~= 1/(1+2*pi*gamma^2*||d_i-d_j||^2).
    %
    % usage : 
    %
    % input:
    %
    %  * D : A matrix with n rows (data points) and d columns (features)
    %  
    %  * sigma : kernel variance parameter. The kernel function approximated
    %  is k(d_i,d_j) = 1/(1+2*pi*gamma^2*||d_i-d_j||^2).
    %
    %  * s : how many random features to use for each f_i
    %
    % output:
    %
    %  * F : A matrix with n rows (data embeddings) and s columns (random Fourier features)
    
   [n,d] = size(D);
   
   gamma = gamma/2;
   
    % set up PDF for 1D sampling
    lim = sqrt(20*gamma^2);
    limstep = lim/1000;
    dgrid = -lim:limstep:lim;
    o = length(dgrid);
    pdf = abs(dgrid/gamma).^(d-1).*exp(-sqrt(dgrid.*dgrid)*sqrt(2)/gamma);
    pdf = pdf/sum(pdf);
    basepdf = abs(dgrid/gamma).^(d-1).*exp(-sqrt(dgrid.*dgrid)*sqrt(2)/gamma);
    basepdf = basepdf/sum(basepdf);
    
    w = zeros(1,s);
    % select lengths of etas
    etas = zeros(1,s);
    for j = 1:s
        y = randsample(o,1,true,pdf);
        etas(j) = dgrid(y)+limstep*(rand(1,1)-.5);
        w(j) = basepdf(y)/pdf(y);
    end
    b = 2*pi*rand(1,s);
    
    % choose random directions for etas
    etavec = randn(d,s);
    etavec = etavec./sqrt(sum(etavec.*etavec,1));
    
    F = sqrt(2)*cos(2*pi*D*(etavec.*etas) + b).*sqrt(w) / sqrt(s);
end 