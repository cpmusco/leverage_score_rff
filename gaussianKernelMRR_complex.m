function F = gaussianKernelMRR(D,gamma,s)
    %% Guassian kernel generator
    % Outputs a submatrix of the Gaussian kernel with variance paramater 
    % gamma for the data rows of D. 
    %
    % usage : 
    %
    % input:
    %
    %  * D : A matrix with n rows (data points) and d columns (features)
    %  
    %  * gamma : kernel variance parameter
    %
    %  * s : how many RR samples to use
    %
    % output:
    %
    %  * F : 
    
    [n,d] = size(D);
    d = 2;
    f = 3;
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