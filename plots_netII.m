
%compare the MC sims with thery with Nc=50 (netTw)

load dAn_netTw_aut %with autaptic connections
load dmc_netTw_aut

% load dAn_netTw %without autaptic connections
% load dmc_netTw

len_vr=length(convg);
cc=jet(len_vr);

figure
hold on
% -- start off values of mn/mx_diag--
mn_diag=mnF_M(1,1);
mx_diag=mn_diag;
for j=len_vr:-1:1
    plot(mnF_M(:,j),mnF_Ma(:,j),'.','MarkerSize',18,'color',cc(j,:))
    % update extreme values for diagonal
    mn_diag=min([mnF_M(:,j);mnF_Ma(:,j);mn_diag]);
    mx_diag=max([mnF_M(:,j);mnF_Ma(:,j);mx_diag]);
end
xl=(mn_diag: (mx_diag-mn_diag)/99: mx_diag)';
plot(xl,xl,'k','LineWidth',1) %diagonal line
set(gca,'FontSize',18)
xlabel('Monte Carlo')
ylabel('Analytic Theory')
title('Mean Firing Rate')

figure
hold on
% -- start off values of mn/mx_diag--
mn_diag=covF_M(1,1,1);
mx_diag=mn_diag;
for j=len_vr:-1:1
    vrM_tmp=diag(covF_M(:,:,j));
    vrAn_tmp=diag(covF_Ma(:,:,j));
    
    plot(vrM_tmp,vrAn_tmp,'.','MarkerSize',18,'color',cc(j,:))
    mn_diag=min([vrM_tmp;vrAn_tmp;mn_diag]);
    mx_diag=max([vrM_tmp;vrAn_tmp;mx_diag]);
end
xl=(mn_diag: (mx_diag-mn_diag)/99: mx_diag)';
plot(xl,xl,'k','LineWidth',1) %diagonal line
set(gca,'FontSize',18)
xlabel('Monte Carlo')
ylabel('Analytic Theory')
title('Variance of Firing Rate')

figure
hold on
% -- start off values of mn/mx_diag--
mn_diag=covF_M(1,2,1);
mx_diag=mn_diag;
for j=len_vr:-1:1
    cvM_tmp=nonzeros(triu(covF_M(:,:,j)-diag(diag(covF_M(:,:,j)))));
    cvAn_tmp=nonzeros(triu(covF_Ma(:,:,j)-diag(diag(covF_Ma(:,:,j)))));
    
    plot(cvM_tmp,cvAn_tmp,'.','MarkerSize',18,'color',cc(j,:))
    % update extreme values for diagonal
    mn_diag=min([cvM_tmp;cvAn_tmp;mn_diag]);
    mx_diag=max([cvM_tmp;cvAn_tmp;mx_diag]);
end
xl=(mn_diag: (mx_diag-mn_diag)/99: mx_diag)';
plot(xl,xl,'k','LineWidth',1) %diagonal line
set(gca,'FontSize',18)
xlabel('Monte Carlo')
ylabel('Analytic Theory')
title('Covariance of Firing Rate')





% ------ repeat but do it for the activity variable X -----
figure
hold on
% -- start off values of mn/mx_diag--
mn_diag=mnX_M(1,1);
mx_diag=mn_diag;
for j=len_vr:-1:1
    plot(mnX_M(:,j),mnX_Ma(:,j),'.','MarkerSize',18,'color',cc(j,:))
    % update extreme values for diagonal
    mn_diag=min([mnX_M(:,j);mnX_Ma(:,j);mn_diag]);
    mx_diag=max([mnX_M(:,j);mnX_Ma(:,j);mx_diag]);
end
xl=(mn_diag: (mx_diag-mn_diag)/99: mx_diag)';
plot(xl,xl,'k','LineWidth',1) %diagonal line
set(gca,'FontSize',18)
xlabel('Monte Carlo')
ylabel('Analytic Theory')
title('Mean Activity')


figure
hold on
% -- start off values of mn/mx_diag--
mn_diag=covX_M(1,1,1);
mx_diag=mn_diag;
for j=len_vr:-1:1
    vrM_tmp=diag(covX_M(:,:,j));
    vrAn_tmp=diag(covX_Ma(:,:,j));
    
    plot(vrM_tmp,vrAn_tmp,'.','MarkerSize',18,'color',cc(j,:))
    mn_diag=min([vrM_tmp;vrAn_tmp;mn_diag]);
    mx_diag=max([vrM_tmp;vrAn_tmp;mx_diag]);
    
end
xl=(mn_diag: (mx_diag-mn_diag)/99: mx_diag)';
plot(xl,xl,'k','LineWidth',1) %diagonal line
set(gca,'FontSize',18)
xlabel('Monte Carlo')
ylabel('Analytic Theory')
title('Variance of Activity')

figure
hold on
% -- start off values of mn/mx_diag--
mn_diag=covX_M(1,2,1);
mx_diag=mn_diag;
for j=len_vr:-1:1
    cvM_tmp=nonzeros(triu(covX_M(:,:,j)-diag(diag(covX_M(:,:,j)))));
    cvAn_tmp=nonzeros(triu(covX_Ma(:,:,j)-diag(diag(covX_Ma(:,:,j)))));
    
    plot(cvM_tmp,cvAn_tmp,'.','MarkerSize',18,'color',cc(j,:))
    % update extreme values for diagonal
    mn_diag=min([cvM_tmp;cvAn_tmp;mn_diag]);
    mx_diag=max([cvM_tmp;cvAn_tmp;mx_diag]);
end
xl=(mn_diag: (mx_diag-mn_diag)/99: mx_diag)';
plot(xl,xl,'k','LineWidth',1) %diagonal line
set(gca,'FontSize',18)
xlabel('Monte Carlo')
ylabel('Analytic Theory')
title('Covariance of Activity')

if(0) %plotting the firing rate correlation
   figure
   hold on
   mn_diag=0;
   mx_diag=mn_diag;
   for j=len_vr:-1:1
       vrM_tmp=diag(covF_M(:,:,j));
       varM_matr=vrM_tmp*vrM_tmp'; %used to divide to get correl matrix
       %calc the corr
       rho_M=covF_M(:,:,j)./sqrt(varM_matr);
    
       vrAn_tmp=diag(covF_Ma(:,:,j));
       varA_matr=vrAn_tmp*vrAn_tmp';
       %calc the corr
       rho_A=covF_Ma(:,:,j)./sqrt(varA_matr);
       
       rhM_tmp=nonzeros(triu(rho_M-diag(diag(rho_M))));
       rhAn_tmp=nonzeros(triu(rho_A-diag(diag(rho_A))));
       
       plot(rhM_tmp,rhAn_tmp,'.','MarkerSize',18,'color',cc(j,:))
       % update extreme values for diagonal
       mn_diag=min([rhM_tmp;rhAn_tmp;mn_diag]);
       mx_diag=max([rhM_tmp;rhAn_tmp;mx_diag]);
   end
   xl=(mn_diag: (mx_diag-mn_diag)/99: mx_diag)';
   plot(xl,xl,'k','LineWidth',1) %diagonal line
   set(gca,'FontSize',18)
   xlabel('Monte Carlo')
   ylabel('Analytic Theory')
   title('Firing Rate Correlation')
end