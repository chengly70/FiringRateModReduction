%script to run larger various tau's
% .. !!! HAVE to run scptAn_testTau.m first to get parms (at least), make sure regime ok

%command loads Nc, mu_vec, sig_vec, tau_vec, rv_vec, sp_vec, Gm, CinMat,corVld
%change file name below IF necessary
load('dAn_testTau_strCoup','Nc','mu_vec','sig_vec','tau_vec','rv_vec','sp_vec','Gm','CinMat','corrVld')

% -- outputs to save --
mnX_M=zeros(Nc,1);
covX_M=zeros(Nc,Nc);
mnF_M=zeros(Nc,1);
covF_M=zeros(Nc,Nc);


tic 
    
[cov_F,mn_F,cov_X,mn_X]=mc_WC(Nc,mu_vec,sig_vec,tau_vec,rv_vec,sp_vec,Gm,CinMat);

mnX_M(:,1)=mn_X;
covX_M=cov_X;
mnF_M=mn_F;
covF_M=cov_F;
    
toc

save dmc_testTau_strCoup covF_M covX_M mnX_M mnF_M