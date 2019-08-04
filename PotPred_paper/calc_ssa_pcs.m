%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    compute a single spectral analysis on PCs time series
%    determined after applying and EOF analysis (using cdo) 
%    on the SSTa 
%
%    Following the idea in Monselesan et al 2015 and O'Kane et al 2016
%
%    With the help of -- copied from!:
%    https://dept.atmos.ucla.edu/tcd/ssa-tutorial-matlab
%
%    Author: Eva C.
%    Created: Oct 2018
%    Last Modif:
%
%    TESTING phase!
%    What about comparing, writing it and using the mcssa toolbox 
%    for SSA analysis i python, also include onte Carlo SSA 
%    (https://github.com/VSainteuf/mcssa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
%addpath(genpath('/data_scratch/matlab_scripts'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF on the yearly mean and 1x1 deg resolution grid
%{
fname_pc = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/yearly/PC_yearlyONmonthly_1deg.nc'
M=12 % dimension of the desired embedding (window length)
%}

%figfile = '/v_Munk_Drive/ecougnon/ana/PotPred/vSZ/SSA/'
figfile = '/home/ecougnon/Documents/MarineX_Project_2017-2019/scripts/tmp_data2plot/'

% EOF on the unfiltered (daily) SSTa and (1/4)x(1/4) deg resolution grid
fname_pc = '/home/ecougnon/Documents/MarineX_Project_2017-2019/scripts/cdo_tool/PCs_daily_trend_1deg_NoWeight.nc'
M=365 % dimension of the desired embedding (window length)
time = ncread(fname_pc,'time');
fname_eigva = '/home/ecougnon/Documents/MarineX_Project_2017-2019/scripts/cdo_tool/eigval_daily_1deg_NoWeight.nc'
eof_pcs_ = ncread(fname_pc,'SSTa');
% Normalizing each Principal Component (dividing by max amplitude of that PC)
amp = max(abs(eof_pcs_),[],2); % max abs value of each principal component
amp2 = repmat(amp,[1,length(eof_pcs_(1,:))]);
eof_pcs_norm = eof_pcs_./amp2;

%{
fname_pc = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/yearly/PCs_YearlyOnMonthly_trend_NoWeight.nc'
M=12 % dimension of the desired embedding (window length)
time = ncread(fname_pc,'time');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose a mode to do the SSA analysis on:
RC_allPCs = nan(8,length(time),3); %12*35
PCs_all = nan(8,length(time)-M+1,3); %12*35-M+1
Lambda_all = nan(8,M);
for mode=1:8
    hmode=mode; %12
    eof_pcs = eof_pcs_norm(mode,:);
%ncread(fname_pc,'SSTa',[hmode,1],[1,Inf]) * ...
%              sqrt(ncread(fname_eigva,'SSTa',[1 1 hmode],[1 1 1]));
    N = length(eof_pcs);
    t = [1:N];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes SSA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the covariance matrix
% with Y the time-delayed embedding of the time series
    Y=zeros(N-M+1,M);
    for m=1:M
        Y(:,m) = eof_pcs((1:N-M+1)+m-1);
    end;
    Cemb=Y'*Y / (N-M+1);

% Calc eigenvalues (LAMBDA) and eigeinvectors (RHO)
    [RHO,LAMBDA] = eig(Cemb);
    LAMBDA = diag(LAMBDA);               % extract the diagonal elements
    [LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
    RHO = RHO(:,ind);                    % and eigenvectors

% percentage of explained variance
    LAMBDA_var = LAMBDA/sum(LAMBDA) * 100;

% Calc PCs
    PC = Y*RHO;

% calc reconstructed components RC
    RC=zeros(N,M);
    for m=1:M
        buf=PC(:,m)*RHO(:,m)'; % invert projection
        buf=buf(end:-1:1,:);
        for n=1:N % anti-diagonal averaging
            RC(n,m)=mean( diag(buf,-(N-M+1)+n) );
        end
    end;

    RC_allPCs(mode,:,:) = RC(:,1:3);
    PCs_all(mode,:,:) = PC(:,1:3);
    Lambda_all(mode,:) = LAMBDA_var;

end
%{
nccreate([figfile 'All_SSA_RCs_YearlyOnMonthly_NoWeight.nc'],'RC_allPCs', ...
         'Dimension',{'modes',8,'time',N,'slow_RC',3},'Format','netcdf4');
nccreate([figfile 'All_SSA_RCs_YearlyOnMonthly_NoWeight.nc'],'PCs_all', ...
         'Dimension',{'modes',8,'time',N,'slow_RC',3},'Format','netcdf4');
nccreate([figfile 'All_SSA_RCs_YearlyOnMonthly_NoWeight.nc'],'Var_expl', ...
         'Dimension',{'modes',8,'Components',M},'Format','netcdf4');
ncwrite([figfile 'All_SSA_RCs_YearlyOnMonthly_NoWeight.nc'],'RC_allPCs',RC_allPCs);
ncwrite([figfile 'All_SSA_RCs_YearlyOnMonthly_NoWeight.nc'],'PCs_all',PCs_all);
ncwrite([figfile 'All_SSA_RCs_YearlyOnMonthly_NoWeight.nc'],'Var_expl',Lambda_all);
%}
%{
nccreate([figfile 'All_SSA_RCs_daily_NoWeight_365.nc'],'RC_allPCs', ...
         'Dimension',{'modes',8,'time',N,'slow_RC',3},'Format','netcdf4');
nccreate([figfile 'All_SSA_RCs_daily_NoWeight_365.nc'],'PCs_all', ...
         'Dimension',{'modes',8,'time',N,'slow_RC',3},'Format','netcdf4');
nccreate([figfile 'All_SSA_RCs_daily_NoWeight_365.nc'],'Var_expl', ...
         'Dimension',{'modes',8,'Components',M},'Format','netcdf4');
ncwrite([figfile 'All_SSA_RCs_daily_NoWeight_365.nc'],'RC_allPCs',RC_allPCs);
ncwrite([figfile 'All_SSA_RCs_daily_NoWeight_365.nc'],'PCs_all',PCs_all);
ncwrite([figfile 'All_SSA_RCs_daily_NoWeight_365.nc'],'Var_expl',Lambda_all);
%}
nccreate([figfile 'All_SSA_RCs_daily_NoWeight_365.nc'],'RC_allPCs', ...
         'Dimension',{'modes',8,'time',N,'slow_RC',3},'Format','netcdf4');
nccreate([figfile 'All_SSA_RCs_daily_NoWeight_365.nc'],'PCs_all', ...
         'Dimension',{'modes',8,'time',N,'slow_RC',3},'Format','netcdf4');
nccreate([figfile 'All_SSA_RCs_daily_NoWeight_365.nc'],'Var_expl', ...
         'Dimension',{'modes',8,'Components',M},'Format','netcdf4');
ncwrite([figfile 'All_SSA_RCs_daily_NoWeight_365.nc'],'RC_allPCs',RC_allPCs);
ncwrite([figfile 'All_SSA_RCs_daily_NoWeight_365.nc'],'PCs_all',PCs_all);
ncwrite([figfile 'All_SSA_RCs_daily_NoWeight_365.nc'],'Var_expl',Lambda_all);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxwidth=800;
maxheight=1100;
C = linspecer(8);
X_axis = [1982:1/365:2017-(1/365)];

%{
for mode=1:8

eof_pcs = ncread(fname_pc,'SSTa',[mode,1],[1,Inf]);

fig=figure('Position',[1 1 maxwidth maxheight-300]),
set(gcf,'name','Variance explained -- eigenvalues')
plot(Lambda_all(mode,:),'o-');
title(['Variance explained for each mode, with mode1: ' num2str(Lambda_all(mode,1)) ' ;mode2: ' num2str(Lambda_all(mode,2)) ' ;mode3: ' num2str(Lambda_all(mode,3)) ''] )
xlim([0,15]);
grid
export_fig([figfile 'Var_explained_Daily_PC' num2str(mode) ''],'-eps','-transparent','-painters')
end
%}


fig=figure('Position',[1 1 maxwidth maxheight]),
set(gcf,'name','Original PCs and their 3 least oscillating RC')
clf;
for mode = 5:5
%    eof_pcs = ncread(fname_pc,'SSTa',[hmode,1],[1,Inf]) * ...
%              sqrt(ncread(fname_eigva,'SSTa',[1 1 hmode],[1 1 1]));

    subplot(4,1,mode-4) % subplot(4,1,mode)
    plot(X_axis,eof_pcs_norm(mode,:),'k-','LineWidth',2);
    hold on,
    plot(X_axis,RC_allPCs(mode,:,1),'-','Color',C(1,:),'LineWidth',2);
    plot(X_axis,RC_allPCs(mode,:,2),'-','Color',C(2,:),'LineWidth',2);
    plot(X_axis,RC_allPCs(mode,:,3),'-','Color',C(4,:),'LineWidth',2);
    legend(['PC' num2str(mode) ''],'RC1','RC2','RC3','Location','SouthEast');
    xlim([1982 2017]);
    title(['PC' num2str(mode) ' time series and 3 least oscillating RCs']);
end
%export_fig([figfile 'RC1-3_daily_SumPC5-6_NoWeight_365'],'-eps','-transparent','-painters')


%{
fig=figure('Position',[1 1 maxwidth maxheight]),
set(gcf,'name','Original PCs and the 3 least oscillating RC combined')
clf;
for mode = 1:4
    subplot(4,1,mode)
%    eof_pcs = ncread(fname_pc,'SSTa',[hmode,1],[1,Inf]) * ...
%              sqrt(ncread(fname_eigva,'SSTa',[1 1 hmode],[1 1 1]));
    plot(X_axis,eof_pcs_norm(mode,:),'-k','LineWidth',2);
    hold on,
    plot(X_axis,sum(RC_allPCs(mode,:,1:3),3),'r-','LineWidth',2);
    legend(['PC' num2str(mode) ''],'Reconstruction with RCs 1-3','Location','SouthEast');
    xlim([1982 2017]);
    title(['PC' num2str(mode) ' time series with the sum of 3 least oscillating RCs']);
end
export_fig([figfile 'RC1-3sum_daily_allPC_NoWeight_365'],'-eps','-transparent','-painters')
%}

%{
figure(3);
set(gcf,'name','Eigenvectors RHO and eigenvalues LAMBDA')
clf;
subplot(3,1,1);
plot(LAMBDA,'o-');
subplot(3,1,2);
plot(RHO(:,1:2), '-');
legend('1', '2');
subplot(3,1,3);
plot(RHO(:,3:4), '-');
legend('3', '4');


figure(4);
set(gcf,'name','Principal components PCs')
clf;
for m=1:4
  subplot(4,1,m);
  plot(t(1:N-M+1),PC(:,m),'k-');
  ylabel(sprintf('PC %d',m));
  ylim([-10000 10000]);
end;

figure(5);
set(gcf,'name','Reconstructed components RCs')
clf;
for m=1:4
  subplot(4,1,m);
  plot(t,RC(:,m),'r-');
  ylabel(sprintf('RC %d',m));
  ylim([-1500 1500]);
end;

figure(6);
set(gcf,'name','Original time series and reconstruction RC')
clf;
subplot(2,1,1)
plot(t,eof_pcs,'b-',t,sum(RC(:,:),2),'r-');
legend('Original','Complete reconstruction');

subplot(2,1,2)
plot(t,eof_pcs,'b','LineWidth',2);
plot(t,eof_pcs,'b-',t,sum(RC(:,1:2),2),'r-');
legend('Original','Reconstruction with RCs 1-2');
%}


