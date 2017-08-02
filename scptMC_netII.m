%script to run larger network, Network II
% .. !!! HAVE to run scptAn_netII.m first to get parms (at least), make
% sure regimes ok (convergence)

%command loads Nc, mu_vec, sig_vec, tau_vec, rv_vec, sp_vec, Gm, CinMat,corVld
load('dAn_netTw_aut','Nc','mu_vec','sig_vec','tau_vec','rv_vec','sp_vec','Gm','CinMat','corrVld')

len_vr=length(corrVld);

% -- outputs to save --
mnX_M=zeros(Nc,len_vr);
covX_M=zeros(Nc,Nc,len_vr);
mnF_M=zeros(Nc,len_vr);
covF_M=zeros(Nc,Nc,len_vr);


tic
for j=1:len_vr
    
    
    [cov_F,mn_F,cov_X,mn_X]=mc_WC(Nc,mu_vec,sig_vec,tau_vec,rv_vec,sp_vec,Gm(:,:,j),CinMat);
    
    mnX_M(:,j)=mn_X;
    covX_M(:,:,j)=cov_X;
    mnF_M(:,j)=mn_F;
    covF_M(:,:,j)=cov_F;
    
end
toc

save dmc_netTw_aut covF_M covX_M mnX_M mnF_M