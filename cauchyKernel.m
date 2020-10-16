function Ksub = cauchyKernel(D,rowInd,colInd,gamma)
    %% Cauchy kernel generator
    % Outputs a submatrix of the Cauchy kernel with width paramater gamma 
    % for the data rows of D. 
    %
    % usage : 
    %
    % input:
    %
    %  * D : A matrix with n rows (data points) and d columns (features)
    %
    %  * rowInd, colInd : Lists of indices between 1 and n. 
    %
    %  NOTE: colInd can be an empty list, in which case the **diagonal** 
    %  entries of the kernel will be output for the indices in rowInd.
    %  
    %  * sigma : kernel variance parameter
    %
    % output:
    %
    %  * Ksub : Let K(i,j) = 1/(1+2*pi*gamma^2*||D(i,:)-D(j,:)||^2). Then Ksub = 
    %  K(rowInd,colInd). Or if colInd = [] then Ksub = diag(K)(rowInd).
    
    if(isempty(colInd))
        Ksub = ones(length(rowInd),1);
    else
        nsqRows = sum(D(rowInd,:).^2,2);
        nsqCols = sum(D(colInd,:).^2,2);
        Ksub = bsxfun(@minus,nsqRows,D(rowInd,:)*(2*D(colInd,:))');
        Ksub = bsxfun(@plus,nsqCols',Ksub);
        Ksub = 1./(1+2*pi*gamma.^2*Ksub);         
    end
end 