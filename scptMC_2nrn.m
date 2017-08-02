%script to test mc_WC.m function on 2 cell network

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
mnX_M=zeros(Nc,len_gs,len_cr);
covX_M=zeros(Nc,Nc,len_gs,len_cr);
mnF_M=zeros(Nc,len_gs,len_cr);
covF_M=zeros(Nc,Nc,len_gs,len_cr);

Gm(2,1)=g_21; %fixed throughout

tic
for ind_gs=1:len_gs
    Gm(1,2)=Gs(ind_gs);
    for ind_cr=1:len_cr
        %ind_cr=2;
        CinMat(1,2)=Corrs(ind_cr);
        CinMat(2,1)=Corrs(ind_cr);
        
        
        [cov_F,mn_F,cov_X,mn_X]=mc_WC(Nc,mu_vec,sig_vec,tau_vec,rv_vec,sp_vec,Gm,CinMat);
        
        mnX_M(:,ind_gs,ind_cr)=mn_X;
        covX_M(:,:,ind_gs,ind_cr)=cov_X;
        mnF_M(:,ind_gs,ind_cr)=mn_F;
        covF_M(:,:,ind_gs,ind_cr)=cov_F;
    end
end
toc

save dmc2_A ind_gs covF_M covX_M mnX_M mnF_M Nc mu_vec sig_vec rv_vec sp_vec tau_vec g_21 Gs Corrs