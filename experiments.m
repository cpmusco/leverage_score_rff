% external MATLAB function used for plotting geopolotical borders
% available directly from https://www.mathworks.com/matlabcentral/fileexchange/50390-borders
addpath('./borders_v3.1.2/borders/')
% precipitation data availabe from https://grass.osgeo.org/download/sample-data/
% the include MATLAB data file was exported using GrassGIS (grass.osgeo.org)
load('slovakia_rain.mat')
% shape file for Slovakian border usedful for plotting
% available directly from https://gadm.org/maps.html
slov = shaperead('./gadm36_SVK_0.shp', 'UseGeoCoords', true);

% clean up data
% fix quirk where first entry in data matrix has spurious value
map_data(1,1) = NaN; 
% remove all NaN rows and columns
good_cols = ~isnan(max(map_data));
good_rows = ~isnan(max(map_data'));
map_data = map_data(find(good_rows),find(good_cols));
map_data = flip(map_data',2);

% move data from grid format to standard data matrix format with latitude
% and longitude values as predictors
[nr, nc] = size(map_data);
lon = linspace(min(slov.Lon), max(slov.Lon), nr);
lat = linspace(min(slov.Lat), max(slov.Lat), nc);
[lat, lon] = meshgrid(lat, lon);
rawDat = zeros(1,5);
c = 1;
for i = 1:nr
    for j = 1:nc
        if(~isnan(map_data(i,j)))
            rawDat(c,:) = [i,j,lat(i,j),lon(i,j),map_data(i,j)];
            c = c + 1;
        end
    end
end
X = rawDat(:,3:4); y = rawDat(:,5);
n = size(X,1);
% mean center target
my = mean(y);
y = y-my;

% split data into training and testing data
rind = randperm(size(rawDat,1));
ntrain = 6400;
rtrain = rind(1:ntrain);
rtest = rind(ntrain+1:end);
rtrain = sort(rtrain);
rtest = sort(rtest);
Xtrain = X(rtrain,:);
ytrain = y(rtrain,:);
Xtest = X(rtest,:);
ytest = X(rtest,:);

% Plot 1: Initial visualization of precipitation data.
figure(); hold();
latlim = ([47.5 50]);
lonlim = ([16 23]);
ax = worldmap(latlim, lonlim);
precip = map_data;
geoshow(lat, lon, precip, 'displaytype', 'texturemap');
colg = [114,149,87]/255;
bordersm('Poland','k','LineWidth',1,'FaceColor',colg);
bordersm('Hungary','k','LineWidth',1,'FaceColor',colg);
bordersm('Austria','k','LineWidth',1,'FaceColor',colg);
bordersm('Czech Republic','k','LineWidth',1,'FaceColor',colg);
bordersm('Ukraine','k','LineWidth',1,'FaceColor',colg);
bordersm('Romania','k','LineWidth',1,'FaceColor',colg);
bordersm('Slovakia','k','LineWidth',4);
colorbar
scatterm(X(rind(1:200),1),X(rind(1:200),2), 'k.')
title('Slovakia Precipitation Data','FontSize',16,'interpreter','latex');
exportgraphics(gca,'visualize_data.png','Resolution',600) 



%% Plot 2 Spectrum Comparison
% construct baseline kernel matrix for xtrain
gamma = 64;
lambda = .008;
K = gaussianKernel(X,rtrain,rtrain,gamma);
Kpred = gaussianKernel(X,1:n,rtrain,gamma);
% compute top 
eigsK = eigs(K,400);
nk = max(abs(eigsK));

% compute approximations using the classical Rahimi-Recht random Fourier 
% features (RR) and our modified random Fourier features with 200, 400, and
% 800 samples
s = 200;
FRFF = gaussianKernelRFF(Xtrain,gamma,s);
FMRFF = gaussianKernelMRFF(Xtrain,gamma,s);
eigsRFF200 = flip(eig(FRFF'*FRFF));
specRFF200 = abs(eigs(K - FRFF*FRFF',1))/nk;
eigsMRFF200 = flip(eig(FMRFF'*FMRFF));
specMRFF200 = abs(eigs(K - FMRFF*FMRFF',1))/nk;

s = 400;
FRFF = gaussianKernelRFF(Xtrain,gamma,s);
FMRFF = gaussianKernelMRFF(Xtrain,gamma,s);
eigsRFF400 = flip(eig(FRFF'*FRFF));
specRFF400 = abs(eigs(K - FRFF*FRFF',1))/nk;
eigsMRFF400 = flip(eig(FMRFF'*FMRFF));
specMRFF400 = abs(eigs(K - FMRFF*FMRFF',1))/nk;

s = 800;
FRFF = gaussianKernelRFF(Xtrain,gamma,s);
FMRFF = gaussianKernelMRFF(Xtrain,gamma,s);
eigsRFF800 = flip(eig(FRFF'*FRFF));
specRFF800 = abs(eigs(K - FRFF*FRFF',1))/nk;
eigsMRFF800 = flip(eig(FMRFF'*FMRFF));
specMRFF800 = abs(eigs(K - FMRFF*FMRFF',1))/nk;

lim = 400;
figure();
p5 = semilogy(1:lim, eigsRFF400(1:lim),'Linewidth',3,'Color',[0.8500, 0.3250, 0.0980]);
hold;
p6 = semilogy(1:lim, eigsRFF800(1:lim),'k:','Linewidth',3,'Color',[0.8500, 0.3250, 0.0980]);
p3 = semilogy(1:lim, eigsMRFF400(1:lim),'Linewidth',3,'Color',[0, 0.4470, 0.7410]);
p4 = semilogy(1:lim, eigsMRFF800(1:lim),'k:','Linewidth',3,'Color',[0, 0.4470, 0.7410]);
p1 = semilogy(1:lim,eigsK(1:lim),'k','Linewidth',1);
legend([p1,p5,p6,p3,p4],'Eigenvalues of Original Kernel', 'Eigenvalues of Classical RFF Approx., 400 Samples','Eigenvalues of Classical RFF Approx., 800 Samples', 'Eigenvalues of Modified RFF Approx., 400 Samples', 'Eigenvalues of Modified RFF Approx., 800 Samples', 'FontSize',14,'interpreter','latex','Location', 'southwest');
set(gca,'fontsize',16)
set(gca,'TickLabelInterpreter','latex');
pbaspect([2 1 1])
xlabel('$i$','FontSize',20,'interpreter','latex');
ylabel('Eigenvalue $\lambda_i$','FontSize',20,'interpreter','latex');
ylim([100*min([eigsRFF400(1:lim);eigsMRFF400(1:lim);eigsK]),3*max([eigsRFF400;eigsMRFF400;eigsK])])
exportgraphics(gca,'spectra.png','Resolution',600) 


%% Plot 3 Spectral Norm Error
z = [specRFF200 specMRFF200; specRFF400 specMRFF400; specRFF800 specMRFF800;];
figure();
bX = categorical({'200 Samples','400 Samples','800 Samples'});
bX = reordercats(bX,{'200 Samples','400 Samples','800 Samples'});
b = bar(bX,z);
for i = 1:size(z,2)
    b(i).FaceColor = 'flat';
end
for i = 1:size(z,1)
    b(1).CData(i,:) = [0.8500 0.3250 0.0980];
    b(2).CData(i,:) = [0 0.4470 0.7410];
end
legend('Classical RFF Approximation', 'Modified RFF Approximation', 'FontSize',16,'interpreter','latex','Location', 'northeast');
set(gca,'fontsize',16)
set(gca,'TickLabelInterpreter','latex');
ylabel('2-Norm Error $\frac{\|K - \tilde{K}\|_2}{\|K\|_2}$','FontSize',20,'interpreter','latex');
pbaspect([2 1 1])
exportgraphics(gca,'spectral_norm_error.png','Resolution',600)  

%% Plot 4 Preconditioning
gamma = 64;
lambda = .008;
K = gaussianKernel(X,rtrain,rtrain,gamma);
Kpred = gaussianKernel(X,1:n,rtrain,gamma);

maxit = 500;
% convergence without preconditioning
[alpha,flag,relres,iter,resvec] = pcg(K + lambda*eye(size(K)),ytrain,1e-14,maxit);

s = 800;
FRFF = gaussianKernelRFF(Xtrain,gamma,s);
FMRFF = gaussianKernelMRFF(Xtrain,gamma,s);

% random Fourier features preconditioner
[U,S,V] = svd(FRFF);
slong = zeros(size(K,1),1);
slong(1:size(S,2)) = diag(S);
iRFF = U*diag(1./(slong.^2+lambda))*U';
RFFfun = @(x) iRFF*x;

% modified random Fourier features preconditioner
[UM,SM,VM] = svd(FMRFF);
slongM = zeros(size(K,1),1);
slongM(1:size(SM,2)) = diag(SM);
iMRFF = UM*diag(1./(slongM.^2+lambda))*UM';
MRFFfun = @(x) iMRFF*x;

[alphaRFF800,flag,relres,iter,resvecRFF800] = pcg(K + lambda*eye(size(K)),ytrain,1e-14,maxit,RFFfun);
[alphaMRFF800,flag,relres,iter, resvecMRFF800] = pcg(K + lambda*eye(size(K)),ytrain,1e-14,maxit,MRFFfun);

s = 400;
FRFF = gaussianKernelRFF(Xtrain,gamma,s);
FMRFF = gaussianKernelMRFF(Xtrain,gamma,s,2);

% random Fourier features preconditioner
[U,S,V] = svd(FRFF);
slong = zeros(size(K,1),1);
slong(1:size(S,2)) = diag(S);
iRFF = U*diag(1./(slong.^2+lambda))*U';
RFFfun = @(x) iRFF*x;

% modified random Fourier features preconditioner
[UM,SM,VM] = svd(FMRFF);
slongM = zeros(size(K,1),1);
slongM(1:size(SM,2)) = diag(SM);
iMRFF = UM*diag(1./(slongM.^2+lambda))*UM';
MRFFfun = @(x) iMRFF*x;

[alphaRFF400,flag,relres,iter,resvecRFF400] = pcg(K + lambda*eye(size(K)),ytrain,1e-14,maxit,RFFfun);
[alphaMRFF400,flag,relres,iter, resvecMRFF400] = pcg(K + lambda*eye(size(K)),ytrain,1e-14,maxit,MRFFfun);

% convergence comparison
miny = min([resvecRFF800;resvecRFF400;resvecMRFF800;resvecMRFF400])/norm(ytrain);
maxy = max([resvecRFF800;resvecRFF400;resvecMRFF800;resvecMRFF400])/norm(ytrain);
figure();
p2 = semilogy(1:length(resvecRFF400),resvecRFF400/norm(ytrain),'Linewidth',3,'Color',[0.8500, 0.3250, 0.0980]);
hold;
p3 = semilogy(1:length(resvecRFF800),resvecRFF800/norm(ytrain),'k:','Linewidth',3,'Color',[0.8500, 0.3250, 0.0980]);
p4 = semilogy(1:length(resvecMRFF400),resvecMRFF400/norm(ytrain),'Linewidth',3,'Color',[0, 0.4470, 0.7410]);
p5 = semilogy(1:length(resvecMRFF800),resvecMRFF800/norm(ytrain),'k:','Linewidth',3,'Color',[0, 0.4470, 0.7410]);
p1 = semilogy(1:length(resvec),resvec/norm(ytrain),'k','Linewidth',2);
legend([p1,p2,p3,p4,p5],'No Preconditioner', 'Classical RFF Precond., 400 Samples', 'Classical RFF Precond., 800 Samples', 'Modified RFF Precond., 400 Samples', 'Modified RFF Precond., 800 Samples', 'FontSize',16,'interpreter','latex','Location', 'southwest');
set(gca,'fontsize',16)
set(gca,'TickLabelInterpreter','latex');
xlabel('Iteration','FontSize',18,'interpreter','latex');
ylabel('Residual Error $\frac{\|(K+\lambda I)x - y\|_2}{\|y\|_2}$','FontSize',18,'interpreter','latex');
ylim([miny,maxy])
xlim([0,500])
pbaspect([1.75 1 1])
exportgraphics(gca,'precond.png','Resolution',600) 


%% Plot 5 Downstream Prediction Error
s = 400;
FRFF = gaussianKernelRFF(Xtrain,gamma,s);
FMRFF = gaussianKernelMRFF(Xtrain,gamma,s,2);

% random Fourier features preconditioner
[U,S,V] = svd(FRFF);
slong = zeros(size(K,1),1);
slong(1:size(S,2)) = diag(S);
iRFF = U*diag(1./(slong.^2+lambda))*U';
RFFfun = @(x) iRFF*x;

% modified random Fourier features preconditioner
[UM,SM,VM] = svd(FMRFF);
slongM = zeros(size(K,1),1);
slongM(1:size(SM,2)) = diag(SM);
iMRFF = UM*diag(1./(slongM.^2+lambda))*UM';
MRFFfun = @(x) iMRFF*x;

ypredBest = Kpred*alphaMRFF800;
best = norm(ypredBest-y)/norm(y)
errsNo = zeros(1,3);
errsRFF = zeros(1,3);
errsMRFF = zeros(1,3);
c = 1;
for it = [50,75,100]
    [alpha,flag,relres,iter,resvec] = pcg(K + lambda*eye(size(K)),ytrain,1e-14,it);
    [alphaRFF400,flag,relres,iter,resvecRFF400] = pcg(K + lambda*eye(size(K)),ytrain,1e-14,it,RFFfun);
    [alphaMRFF400,flag,relres,iter, resvecMRFF400] = pcg(K + lambda*eye(size(K)),ytrain,1e-14,it,MRFFfun);

    ypredNo = Kpred*alpha;
    ypredRFF = Kpred*alphaRFF400;
    ypredMRFF = Kpred*alphaMRFF400;
    
    errsNo(c) = norm(ypredNo-y)/norm(y);
    errsRFF(c) = norm(ypredRFF-y)/norm(y);
    errsMRFF(c) = norm(ypredMRFF-y)/norm(y);
    c = c+1;
end

errsNo(2:3) = errsNo(1:2)
errsRFF(2:3) = errsRFF(1:2)
errsMRFF(2:3) = errsMRFF(1:2)

z = [errsNo(1) errsRFF(1) errsMRFF(1) best; errsNo(2) errsRFF(2) errsMRFF(2) best; errsNo(3) errsRFF(3) errsMRFF(3) best;].^2;
figure();
bX = categorical({'25 Iterations','50 Iterations','75 Iterations'});
bX = reordercats(bX,{'25 Iterations','50 Iterations','75 Iterations'});
b = bar(bX,z);
for i = 1:size(z,2)
    b(i).FaceColor = 'flat';
end
for i = 1:size(z,1)
	b(1).CData(i,:) = [.5 .5 .5];
    b(2).CData(i,:) = [0.8500 0.3250 0.0980];
    b(3).CData(i,:) = [0 0.4470 0.7410];
    b(4).CData(i,:) = [0 0 0];
end

legend('No Preconditioner', 'Classical RFF Precond., 400 Samples', 'Modified RFF Precond., 400 Samples', 'Best Possible Error', 'FontSize',16,'interpreter','latex','Location', 'northeast');
set(gca,'fontsize',16)
set(gca,'TickLabelInterpreter','latex');
ylabel('Prediction Error $\frac{\|y_{pred} - y_{test}\|^2_2}{\|y_{test}\|^2_2}$','FontSize',18,'interpreter','latex');
ylim([0,max(max(z))*1.1])
pbaspect([1.75 1 1])
exportgraphics(gca,'regression_error.png','Resolution',600) 


%% Plot 6 Visualization of Final Interpolation Result
map_pred = map_data;
for i = 1:length(ypredBest)
    map_pred(rawDat(i,1),rawDat(i,2)) = ypredBest(i);
end
figure(); hold();
latlim = ([47.5 50]);
lonlim = ([16 23]);
ax = worldmap(latlim, lonlim);
precip = map_pred + my;
geoshow(lat, lon, precip, 'displaytype', 'texturemap');
colg = [114,149,87]/255;
bordersm('Poland','k','LineWidth',1,'FaceColor',colg);
bordersm('Hungary','k','LineWidth',1,'FaceColor',colg);
bordersm('Austria','k','LineWidth',1,'FaceColor',colg);
bordersm('Czech Republic','k','LineWidth',1,'FaceColor',colg);
bordersm('Ukraine','k','LineWidth',1,'FaceColor',colg);
bordersm('Romania','k','LineWidth',1,'FaceColor',colg);
bordersm('Slovakia','k','LineWidth',4);
colorbar
title('Slovakia Precipitation, Predicted','FontSize',16,'interpreter','latex');
exportgraphics(gca,'visualize_data_predict.png','Resolution',600) 

%% Plot 7 Visualization of Modified Sampling Distributions
% gaussian
lim = 5
step = .1
[Xsamp,Ysamp] = meshgrid(-lim:step:lim);
Z = sqrt(Xsamp.^2 + Ysamp.^2);
f = 40;
Z = 1 - exp(f*(Z-3))./(exp(f*(Z-3)) + 1);
figure();
sPlot = surf(Xsamp,Ysamp,Z,'FaceAlpha',.5,'FaceColor',[0, 0.4470, 0.7410])
ylim([-lim,lim])
xlim([-lim,lim])
zlim([0,max(max(Z))*1.1])
exportgraphics(gca,'mrr_gauss_dist.png','Resolution',600) 

[Xsamp,Ysamp] = meshgrid(-5:.1:5);
Z = exp(-(Xsamp.^2 + Ysamp.^2));
figure();
sPlot = surf(Xsamp,Ysamp,Z,'FaceAlpha',.5,'FaceColor',[0.8500, 0.3250, 0.0980])
ylim([-lim,lim])
xlim([-lim,lim])
zlim([0,max(max(Z))*1.1])
exportgraphics(gca,'rr_gauss_dist.png','Resolution',600) 

% cauchy
lim = 8
step = .15
[Xsamp,Ysamp] = meshgrid(-lim:step:lim);
Z = 1./(1 + sqrt(Xsamp.^2 + Ysamp.^2))
M = sqrt(Xsamp.^2 + Ysamp.^2) < 7.5;
Z = Z.*M
figure();
sPlot = surf(Xsamp,Ysamp,Z,'FaceAlpha',.5,'FaceColor',[0, 0.4470, 0.7410])
ylim([-lim,lim])
xlim([-lim,lim])
zlim([0,max(max(Z))*1.1])
exportgraphics(gca,'mrr_exp_dist.png','Resolution',600) 

Z = exp(-sqrt(Xsamp.^2 + Ysamp.^2));
figure();
sPlot = surf(Xsamp,Ysamp,Z,'FaceAlpha',.5,'FaceColor',[0.8500, 0.3250, 0.0980])
ylim([-lim,lim])
xlim([-lim,lim])
zlim([0,max(max(Z))*1.1])
exportgraphics(gca,'rr_exp_dist.png','Resolution',600) 