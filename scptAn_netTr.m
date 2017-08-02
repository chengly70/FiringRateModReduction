%script for theory Network III, clustered coupling in E; 50E's/50I's

Nc=100; %larger network
N_e=50;
N_i=50; %N_e+N_i=Nc
tau_vec=1+.075*randn(Nc,1); %diff than scptAn_netII.m
mu_vec=2*rand(Nc,1) - 1; %[-1,1]
sig_vec=rand(Nc,1) + 1;

rv_vec=0.1*randn(Nc,1);
sp_vec=0.4*rand(Nc,1) + 0.05; %diff than scptAn_netII.m
%CinMat is diff than scptAn_netII.m
Cin_ee=diag(.5*ones(N_e,1))+diag(0.1+.1*randn(N_e-1,1),1);
Cin_ii=diag(.5*ones(N_e,1))+diag(0.12+.1*randn(N_i-1,1),1);
Cin_ei=zeros(N_e,N_i);
for j=1:N_e
    Cin_ei(N_e-(j-1),j)=0.3+.1*randn;
end
%Cin_sub1=[Cin_ee zeros(N_i,N_e);Cin_ei Cin_ii];
Cin_sub1=[Cin_ee Cin_ei; zeros(N_e,N_i) Cin_ii];
CinMat=Cin_sub1+Cin_sub1';
R=chol(CinMat);

len_vr=4; %# times vary params
Gm=zeros(Nc,Nc,len_vr); %vary coupling matrix, randomly chosen
gm_tmp=zeros(Nc,Nc);
clust_blk=ones(10,10)-diag(ones(10,1)); %create clustered connect in 5 segments
for j=1:4
    clust_blk=blkdiag(clust_blk,ones(10,10)-diag(ones(10,1)));
end%remove autaptic coupling
pct_ItoE=0.35;
pct_EtoI=0.35;
pct_ItoI=0.35;

% -- outputs to save --
mnX_Ma=zeros(Nc,len_vr);
covX_Ma=zeros(Nc,Nc,len_vr);
mnF_Ma=zeros(Nc,len_vr);
covF_Ma=zeros(Nc,Nc,len_vr);
convg=zeros(len_vr,1);
corrVld=zeros(len_vr,1);
g_ItoE=zeros(len_vr,1);
g_EtoI=zeros(len_vr,1);
g_ItoI=zeros(len_vr,1);
g_EtoE=zeros(len_vr,1);
tic
for j=1:len_vr
    g_EtoE(j,1)=(rand)/(N_e*1/5);
    g_ItoE(j,1)=(-6*rand-2)/(N_i*pct_ItoE);
    g_EtoI(j,1)=(6*rand+2)/(N_e*pct_EtoI);
    g_ItoI(j,1)=(-6*rand-2)/(N_i*pct_ItoI);
    
    gm_tmp=[clust_blk*g_EtoE(j,1) (rand(N_e,N_i)<pct_ItoE)*g_ItoE(j,1) ;...
        (rand(N_i,N_e)<pct_EtoI)*g_EtoI(j,1) (rand(N_i,N_i)<pct_ItoI)*g_ItoI(j,1)];
    %remove autaptic coupling
    gm_tmp(N_e+1:end,N_e+1:end)=gm_tmp(N_e+1:end,N_e+1:end)-diag(diag(gm_tmp(N_e+1:end,N_e+1:end)));
    
    Gm(:,:,j)=sparse(gm_tmp);
    
    [convged,Corr_valid,cov_Fa,mn_Fa,cov_Xa,mn_Xa,mnAll]=iter_method(Nc,mu_vec,sig_vec,tau_vec,rv_vec,sp_vec,Gm(:,:,j),CinMat);
    
    convg(j,1)=convged;
    corrVld(j,1)=Corr_valid;
    mnX_Ma(:,j)=mn_Xa;
    covX_Ma(:,:,j)=cov_Xa;
    mnF_Ma(:,j)=mn_Fa;
    covF_Ma(:,:,j)=cov_Fa;
    norm(mnAll(:,end-1)-mnAll(:,end))
end
toc

save dAn_netEI covF_Ma covX_Ma mnX_Ma mnF_Ma convg corrVld Nc mu_vec sig_vec tau_vec rv_vec sp_vec Gm CinMat g_* pct_* N_*