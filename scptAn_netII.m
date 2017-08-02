%script to run iter_method.m function on large random 50 cell network
%(dense)

Nc=50; %larger network
tau_vec=1+.05*randn(Nc,1);
mu_vec=2*rand(Nc,1) - 1; %[-1,1]
sig_vec=rand(Nc,1) + 1;

rv_vec=0.1*randn(Nc,1);
sp_vec=0.35*rand(Nc,1) + 0.05;
Cin_sub=0.8*randn(Nc,Nc);
Cin_sub=Cin_sub'*Cin_sub; %make symme
scl_dia=diag(1./sqrt(diag(Cin_sub))); %diag matrix to put 1's diag CinMat
CinMat=scl_dia*Cin_sub*scl_dia;  %should be positive definite, check:
R=chol(CinMat);

len_vr=4; %# times vary params
Gm=zeros(Nc,Nc,len_vr); %vary coupling matrix, randomly chosen
gm_tmp=zeros(Nc,Nc);

% -- outputs to save --
mnX_Ma=zeros(Nc,len_vr);
covX_Ma=zeros(Nc,Nc,len_vr);
mnF_Ma=zeros(Nc,len_vr);
covF_Ma=zeros(Nc,Nc,len_vr);
convg=zeros(len_vr,1);
corrVld=zeros(len_vr,1);

tic
for j=1:len_vr
    %change Gm or ..?
    gm_tmp=j*.1*randn(Nc,Nc);
    %Gm(:,:,j)=gm_tmp-diag(diag(gm_tmp));
    Gm(:,:,j)=gm_tmp;
    
    [convged,Corr_valid,cov_Fa,mn_Fa,cov_Xa,mn_Xa]=iter_method(Nc,mu_vec,sig_vec,tau_vec,rv_vec,sp_vec,Gm(:,:,j),CinMat);
    
    convg(j,1)=convged;
    corrVld(j,1)=Corr_valid;
    mnX_Ma(:,j)=mn_Xa;
    covX_Ma(:,:,j)=cov_Xa;
    mnF_Ma(:,j)=mn_Fa;
    covF_Ma(:,:,j)=cov_Fa;
 
end
toc

save dAn_netTw_aut covF_Ma covX_Ma mnX_Ma mnF_Ma convg corrVld Nc mu_vec sig_vec tau_vec rv_vec sp_vec Gm CinMat