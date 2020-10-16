%% computing upper bound for leverage scores under Gaussian measure
% discretization step size for approximation k-sparse Fourier operator
step = .01;
% range to approximation leverage scores over
r = 8;
x = [-r:step:r];
% sparsity of function class
k = 5;
% defining gaussian distribution g
sigma = 1;
g = (1/sqrt(2*pi*sigma^2))*exp(-x.^2/(2*sigma^2))*step;

% the total number of iterations run is outerloops*n. In the paper
% outerloops was set to 1000, so we ran for 10 million iterations.
outerloops = 5;
n = 10000;

% frequencies will be randomly sampled, with variance randomly chosen on a
% geometric scale for the width vector.
width = .01*2.^(0:10);
nw = length(width);

% ubtight will maintain current estimates for the leverage scores at x
% after running the loop below it should be a pretty close estimate
ubtight = zeros(size(x));
% maxlines tracks the Fourier function found with largest magnitude at x
maxlines = zeros(length(x),length(x));
for o = 1:outerloops
    % print out to keep track of what iteration we're on
    o
for i = 1:n
    w = randi(nw);
    % randomly chosen frequencies
    lambs = width(w)*randn(k,1);
    % discretized approximation to the k-sparse Fourier operator
    A = exp(-sqrt(-1)*(lambs*x)').*(sqrt(g)'*ones(size(lambs))');
    coeffs = conj(A)*pinv(A'*A);
    % exact leverage scores for A
    levs = sum(A.*coeffs,2)';
    slevs = levs/step;
    ubtight = max(ubtight,slevs);   
    for j = find(ubtight-slevs == 0)
        % for any x where we found a higher leverage score estimate, we
        % reconstruct and store the function f that realizes this score
        f = A*coeffs(j,:)';
        sf = (1/step)*abs(f).^2/norm(f)^2;
        maxlines(:,j) = sf';
    end
end
end


%% computing upper bound for leverage scores under Laplace measure
% (same code as above, with some different parameters)
step = .01;
rz = 16;
xz = [-rz:step:rz];
sigma = 1;
k = 5;
z = (1/sqrt(2*sigma^2))*exp(-sqrt(2)*abs(xz)/sigma)*step;

outerloops = 5;
n = 10000;
width = .01*2.^(0:10);
nw = length(width);

ubtightz = zeros(size(xz));
maxlinesz = zeros(length(xz),length(xz));
for o = 1:outerloops
    % print out to keep track of what iteration we're on
    o
for i = 1:n
    w = randi(nw);
    lambs = width(w)*randn(k,1);
    A = exp(-sqrt(-1)*(lambs*xz)').*(sqrt(z)'*ones(size(lambs))');
    coeffs = conj(A)*pinv(A'*A);
    levs = sum(A.*coeffs,2)';
    slevs = levs/step;
    ubtightz = max(ubtightz,slevs);   
    for j = find(ubtightz-slevs == 0)
        f = A*coeffs(j,:)';
        sf = (1/step)*abs(f).^2/norm(f)^2;
        maxlinesz(:,j) = sf';
    end
end
end

%% reconstruct Figure 1 from the paper for Gaussian and Laplace leverage scores
% gaussian
figure();
hold;
num_dark_plots = 1;
num_light_plots = 200;
for i = 1:num_light_plots
    plot(x,real(maxlines(:,randi(length(x)))),'Color',[0.6 0.8 1],'linewidth',.25);
end
pl(2) = plot(x,ubtight,'k','linewidth',2);
pl(1) = plot(x,real(maxlines(:,900)),'Color',[0 0 1],'linewidth',2);
gauss_up = zeros(size(x));
gauss_up = (1/sqrt(4*sigma^2))*exp(-x.^2/(4*sigma^2))/step;
m = abs(x) < 4.5;
gauss_up(m) = 2.1;
pl(3) = plot(x,gauss_up, 'r-.', 'linewidth',2);

set(gca,'fontsize',16)
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$','FontSize',20,'interpreter','latex');
legend(pl)
legend(pl,'Example of $\frac{|f(x)|^2\cdot g(x)}{\|f\|_g^2}$ for $f\in \mathcal{T}_5$.','Empirical estimate for $\tau_{5,g}(x)$.','Possible closed form upper bound $\bar{\tau}_{5,g}(x)$.', 'FontSize',16,'interpreter','latex','Location', 'north');
ylim([0,1.4*max(abs(gauss_up))])
xlim([-r,r])

% laplace
figure();
hold;
num_dark_plots = 1;
num_light_plots = 200;
for i = 1:num_light_plots
    plot(xz,real(maxlinesz(:,randi(length(xz)))),'Color',[0.6 0.8 1],'linewidth',.25);
end
pl(2) = plot(xz,ubtightz,'k','linewidth',2);
pl(1) = plot(xz,real(maxlinesz(:,1700)),'Color',[0 0 1],'linewidth',2);
laplace_up = 1/sqrt(2*sqrt(8)*sigma)^2*exp(-sqrt(2)*abs(xz)/(sqrt(8)*sigma))/step;
m = abs(xz) < 8;
laplace_up(m) = 4./(1 + abs(xz(m))*1);
pl(3) = plot(xz,laplace_up, 'r-.', 'linewidth',2);

set(gca,'fontsize',16)
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$','FontSize',20,'interpreter','latex');
legend(pl)
legend(pl,'Example of $\frac{|f(x)|^2\cdot z(x)}{\|f\|_z^2}$ for $f\in \mathcal{T}_5$.','Empirical estimate for $\tau_{5,z}(x)$.','Possible closed form upper bound $\bar{\tau}_{5,z}(x)$.','FontSize',16,'interpreter','latex','Location', 'north');
ylim([0,1.35*max(abs(laplace_up))])
xlim([-rz,rz])