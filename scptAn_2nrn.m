%script to run iter_method.m function on 2 cell network

load d2FF_b
Nc=2; %2 cell, Feedforward inputs
tau_vec=ones(Nc,1);
%the following variables are assumed to be saved in .mat loaded above:
%mu1,mu2,sig1,sig2,s_rv,s_sp,Gs,Corrs
mu_vec=[mu1;mu2];
sig_vec=[sig1;sig2];
rv_vec=ones(Nc,1)*s_rv;
sp_vec=ones(Nc,1)*s_sp;
Gm=zeros(Nc,Nc);
CinMat=diag(ones(Nc,1));

Corrs=(0:0.2:0.8)'; %overwrite d2FF_b value
len_cr=length(Corrs);
g_21=0.4;
Gs=(-2:0.2:2)';     %overwrite d2FF_b value
len_gs=length(Gs);

% -- outputs to save --
mnX_Ma=zeros(Nc,len_gs,len_cr);
covX_Ma=zeros(Nc,Nc,len_gs,len_cr);
mnF_Ma=zeros(Nc,len_gs,len_cr);
covF_Ma=zeros(Nc,Nc,len_gs,len_cr);
convg=zeros(len_gs,len_cr);
corrVld=zeros(len_gs,len_cr);

Gm(2,1)=g_21; %fixed throughout

tic
for ind_gs=1:len_gs
    Gm(1,2)=Gs(ind_gs);
    for ind_cr=1:len_cr
        CinMat(1,2)=Corrs(ind_cr);
        CinMat(2,1)=Corrs(ind_cr);
        
        
        [convged,Corr_valid,cov_Fa,mn_Fa,cov_Xa,mn_Xa]=iter_method(Nc,mu_vec,sig_vec,tau_vec,rv_vec,sp_vec,Gm,CinMat);
        
        convg(ind_gs,ind_cr)=convged;
        corrVld(ind_gs,ind_cr)=Corr_valid;
        mnX_Ma(:,ind_gs,ind_cr)=mn_Xa;
        covX_Ma(:,:,ind_gs,ind_cr)=cov_Xa;
        mnF_Ma(:,ind_gs,ind_cr)=mn_Fa;
        covF_Ma(:,:,ind_gs,ind_cr)=cov_Fa;
    end
end
toc
save dAn2_A covF_Ma covX_Ma mnX_Ma mnF_Ma convg corrVld Gs Corrs g_21