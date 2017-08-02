%script to test mc_WC.m function for larger tau_vec values, produced
%dAn_testTau (Gm std=0.1) and dAn_testTau_strCoup (Gm std=0.25)

Nc=50; %larger network
tau_vec=(.5: 4.5/(Nc-1) : 5)';
mu_vec=0.7*ones(Nc,1);
sig_vec=1.3*ones(Nc,1);

rv_vec=0.1*ones(Nc,1);
sp_vec=0.35*ones(Nc,1);
CinMat=0.3*diag(ones(Nc-1,1),1)+0.3*diag(ones(Nc-1,1),-1) + diag(ones(Nc,1));

Gm=zeros(Nc,Nc); %vary coupling matrix, randomly chosen

% -- outputs to save --
mnX_Ma=zeros(Nc,1);
covX_Ma=zeros(Nc,Nc);
mnF_Ma=zeros(Nc,1);
covF_Ma=zeros(Nc,Nc);
convg=0;
corrVld=0;

tic

%change Gm or ..?
%Gm=0.1*randn(Nc,Nc); %weaker coupling
Gm=0.25*randn(Nc,Nc); %stronger coupling

[convged,Corr_valid,cov_Fa,mn_Fa,cov_Xa,mn_Xa]=iter_method(Nc,mu_vec,sig_vec,tau_vec,rv_vec,sp_vec,Gm,CinMat);

convg=convged;
corrVld=Corr_valid;
mnX_Ma(:,1)=mn_Xa;
covX_Ma=cov_Xa;
mnF_Ma=mn_Fa;
covF_Ma=cov_Fa;
 
toc

save dAn_testTau_strCoup covF_Ma covX_Ma mnX_Ma mnF_Ma convg corrVld Nc mu_vec sig_vec tau_vec rv_vec sp_vec Gm CinMat