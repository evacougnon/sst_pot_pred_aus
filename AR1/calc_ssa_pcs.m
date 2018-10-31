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
addpath(genpath('/data_scratch/matlab_scripts'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF on the yearly mean and 1x1 deg resolution grid
%{
fname_pc = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/yearly/PC_yearlyONmonthly_1deg.nc'
M=12 % dimension of the desired embedding (window length)
%}

figfile = '/v_Munk_Drive/ecougnon/ana/PotPred/vSZ/SSA/'


% EOF on the unfiltered (daily) SSTa and (1/4)x(1/4) deg resolution grid
fname_pc = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/PCs_daily_trend_1deg.nc'
M=360 % dimension of the desired embedding (window length)
time = ncread(fname_pc,'time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose a mode to do the SSA analysis on:
RC_allPCs = nan(8,length(time),3); %12*35
PCs_all = nan(8,length(time)-M+1,3); %12*35-M+1
Lambda_all = nan(8,M);
for mode=1:8
hmode=mode; %12
eof_pcs = ncread(fname_pc,'SSTa',[hmode,1],[1,Inf]);
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
nccreate([figfile 'All_SSA_RCs_Daily.nc'],'RC_allPCs', ...
         'Dimension',{'modes',8,'time',N,'slow_RC',3},'Format','netcdf4');
nccreate([figfile 'All_SSA_RCs_Daily.nc'],'PCs_all', ...
         'Dimension',{'modes',8,'time',N,'slow_RC',3},'Format','netcdf4');
nccreate([figfile 'All_SSA_RCs_Daily.nc'],'Var_expl', ...
         'Dimension',{'modes',8,'Components',M},'Format','netcdf4');
ncwrite([figfile 'All_SSA_RCs_Daily.nc'],'RC_allPCs',RC_allPCs);
ncwrite([figfile 'All_SSA_RCs_Daily.nc'],'PCs_all',PCs_all);
ncwrite([figfile 'All_SSA_RCs_Daily.nc'],'Var_expl',Lambda_all);




%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxwidth=1300;
maxheight=600;
C = linspecer(8);

fig=figure('Position',[1 1 maxwidth maxheight-300]),
set(gcf,'name','Variance explained -- eigenvalues')
plot(LAMBDA_var,'o-');
title(['Variance explained for each mode, with mode1: ' num2str(LAMBDA_var(1)) ' ;mode2: ' num2str(LAMBDA_var(2)) ' ;mode3: ' num2str(LAMBDA_var(3)) ''] )
xlim([0,15]);
grid
export_fig([figfile 'Var_explained_Daily_PC' num2str(hmode) ''],'-eps','-transparent','-painters')

fig=figure('Position',[1 1 maxwidth maxheight]),
set(gcf,'name','Original time series and 3 least oscillating RC')
clf;
subplot(2,1,1)
plot(t,eof_pcs,'k-','LineWidth',2);
hold on,
plot(t,RC(:,1),'-','Color',C(1,:),'LineWidth',2);
plot(t,RC(:,2),'-','Color',C(2,:),'LineWidth',2);
plot(t,RC(:,3),'-','Color',C(4,:),'LineWidth',2);
legend('Original','RC1','RC2','RC3','Location','SouthEast');
title('Original time series and 3 least oscillating RCs');
subplot(2,1,2)
plot(t,eof_pcs,'-k','LineWidth',2);
hold on,
plot(t,sum(RC(:,1:3),2),'r-','LineWidth',2);
legend('Original','Reconstruction with RCs 1-3','Location','SouthEast');
title('Original time series with the sum of 3 least oscillating RCs');
export_fig([figfile 'RC1-3_Original_Daily_PC' num2str(hmode) ''],'-eps','-transparent','-painters')
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


