function [var_F,cov_F,mn_F,var_X,cov_X,mn_X]=mc_WC_large(Nc,mu_vec,sig_vec,tau_vec,rv_vec,sp_vec,Gm,CinMat)
%very much like mc_WC.m but exploiting sparsity structure of Gm,CinMat
% ALSO, outputs split into var & cov, so diag is of cov is NOT var
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
covM_xi=sparse(CinMat.*(sig_vec*sig_vec')); %pure cov matrix
R_cor=chol(covM_xi); %R_cor'*randn(Nc,1) gives correct noise term (RHS)

%indices of cov/corr save, assuming same structure as CinMat
crMat_p=sparse(triu(CinMat-diag(diag(CinMat)))); %entries to save
ind_nnzCr=find(crMat_p);   %linear indices of non-zero correl
[rw_ind,cl_ind]=ind2sub([Nc Nc],ind_nnzCr); %row and column indices of non-zero correl
num_Cvs=length(ind_nnzCr);

%OUTPUTS
mn_X=zeros(Nc,1);    %mean of Xj
var_X=zeros(Nc,1);   %var of Xj
cov_X=zeros(num_Cvs,1);  %subset of cov matrix of Xj (no diag)
mn_F=zeros(Nc,1);    %mean of F(Xj)
var_X=zeros(Nc,1);   %var of F(Xj)
cov_F=zeros(num_Cvs,1);  %subset of cov matrix of Fj (no diag)

%reset vals
mu_rsV=zeros(N_relz,Nc);
mu_Fv=zeros(N_relz,Nc); %firing rates

%cov parts
varInst_F=zeros(Nc,1); %all cross covs F(Xj)
varInst_X=zeros(Nc,1); %all cross covs Xj
covInst_F=zeros(num_Cvs,1); %all cross covs F(Xj)
covInst_X=zeros(num_Cvs,1); %all cross covs Xj

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
    
    varInst_X=(varInst_X*(j-2) + sum(xV.^2)')./(j-1);
    varInst_F=(varInst_F*(j-2) + sum(FxMat.^2)')./(j-1);
    
    covInst_X=(covInst_X*(j-2) + sum(xV(:,rw_ind).*xV(:,cl_ind))')./(j-1);     %running sum; entire cov matrix
    covInst_F=(covInst_F*(j-2) + sum(FxMat(:,rw_ind).*FxMat(:,cl_ind))')./(j-1); %running sum; entire cov matrix
    
end
%normalize and store stats/dens
mn_X=sum(mu_rsV)'./(N_relz*Lt);
mn_F=sum(mu_Fv)'./(N_relz*Lt);

covInst_F=covInst_F./N_relz;
covInst_X=covInst_X./N_relz;
varInst_F=varInst_F./N_relz;
varInst_X=varInst_X./N_relz;

var_X=varInst_X-mn_X.^2;
var_F=varInst_F-mn_F.^2;
cov_F=covInst_F-mn_F(rw_ind).*mn_F(cl_ind);
cov_X=covInst_X-mn_X(rw_ind).*mn_X(cl_ind);


