function [cov_F,mn_F,cov_X,mn_X]=mc_WC(Nc,mu_vec,sig_vec,tau_vec,rv_vec,sp_vec,Gm,CinMat)
%
%Simulate full network of Nc (input) WC, coupled/corrNoise, etc
%!mu_vec, sig_vec, tau_vec are all Nc x 1 column vectors!
%F is sigmoidal; 1/2*(1+tanh((Inp-rv_vec)./sp_vec))
%Gm is an NcxNc coupling matrix, with no autaptic (diag(Gm)=0)
%CinMat is an NcxNc correlation matrix; ones on diag and PSD!

%--- should check: (mu_vec,sig_vec,tau_vec,rv_vec,sp_vec) are all Nc x 1
%---    CinMat is PSD, Gm is Nc x Nc ---

rng('shuffle') %seed random number generator

dt=0.01; % in msec
t_end=500; % in msec
Lt=floor(t_end/dt)+1; %total num time-steps
sq_dt=1/sqrt(dt);

N_relz=5000;

%state variables (variables needed to do sim)
mu_Mat=repmat(mu_vec',N_relz,1);
xV=randn(N_relz,Nc)*diag(sig_vec./sqrt(2))+mu_Mat; %# of cells
%sig_Mat=repmat(sig_vec',N_relz,1);
rv_Mat=repmat(rv_vec',N_relz,1);
sp_Mat=repmat(sp_vec',N_relz,1);
tauInv_Mat=repmat(1./tau_vec',N_relz,1);
xi=zeros(N_relz,Nc); %white noise forcing
covM_xi=CinMat.*(sig_vec*sig_vec'); %pure cov matrix
R_cor=chol(covM_xi); %R_cor'*randn(Nc,1) gives correct noise term (RHS)

%OUTPUTS
mn_X=zeros(Nc,1);    %mean of Xj
cov_X=zeros(Nc,Nc);  %entire cov matrix of Xj (includes var)
mn_F=zeros(Nc,1);    %mean of F(Xj)
cov_F=zeros(Nc,Nc);  %entire cov matrix of Fj (includes var)

%reset vals
mu_rsV=zeros(N_relz,Nc);
mu_Fv=zeros(N_relz,Nc); %firing rates

%cov parts
covInst_F=zeros(Nc,Nc); %all cross covs F(Xj)
covInst_X=zeros(Nc,Nc); %all cross covs Xj

%in main loop, since N_relz x Nc, use (L*xi')'=xi*L', where L'=R_cor
% similarly for (Gm*FxMat')'=FxMat*Gm'
for j=2:Lt
    xi=randn(N_relz,Nc); %generate N_relz standard normal (uncorr in time)
    
    FxMat = 0.5*(1+tanh((xV-rv_Mat)./sp_Mat));
    
    %eqns for cells
    xV=xV+dt*tauInv_Mat.*( -xV+mu_Mat+sq_dt*(xi*R_cor)+FxMat*Gm' );
    
    %keep running sum of stats/dens
    mu_rsV=mu_rsV+xV;
    mu_Fv=mu_Fv+FxMat;
    
    covInst_X=(covInst_X*(j-2) + xV'*xV)./(j-1);     %running sum; entire cov matrix
    covInst_F=(covInst_F*(j-2) + FxMat'*FxMat)./(j-1); %running sum; entire cov matrix
    
end
%normalize and store stats/dens
mn_X=sum(mu_rsV)'./(N_relz*Lt);
mn_F=sum(mu_Fv)'./(N_relz*Lt);

covInst_F=covInst_F./N_relz;
covInst_X=covInst_X./N_relz;
cov_F=covInst_F-mn_F*mn_F';
cov_X=covInst_X-mn_X*mn_X';


