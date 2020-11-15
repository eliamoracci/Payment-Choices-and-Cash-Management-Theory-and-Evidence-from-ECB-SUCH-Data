%%% This code should solve a model of dynamic payment choice
%%% Elia Moracci
%%% 31.05.2020

% Housekeeping
clear all
close all
clc
% Setting LaTeX default interpreter for text in plots
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'DefaultLineLineWidth', 2);

newcolors = [0.65 0.0 0.0
             0.40 0.4 0.4
             0.75 0.75 0.75
             0.0 0.0 0.65];
tic
% Plug in already-estimated values and moments from STATA
% data (ECB Survey on the Use of Cash by Households).
phi         = 0.7;              % Average European value
% mu_unc    = 1.948;            % Calibrated using STATA on all payments<250€
% sigma_unc = 1.277;            % Calibrated using STATA on all payments<250€

%% Solve model with estimated parameters

beta      = 0.941;
eta_h     = 2.047;
eta_l     = 2.047;
R         = 0.009;
sigma     = 0;
kappa     = 0.461;
mu_unc    = 2.3915;
sigma_unc = 0.8789;

N         = 1;
T         = 1;
[V,pol,money,payment,size,grid_m,grid_s]=...
    SolveModel(beta,eta_h,eta_l,R,sigma,kappa,phi,mu_unc,sigma_unc,N,T);
%% Produce Figures

% figure(1)
% subplot(2,1,1)
% contourf(grid_m,grid_s,pol.adj_yn',[0 1]);
% colormap(gray)
% grid on
% colorbar
% xlabel('$m$')
% ylabel('$s$')
% title('Adjustment decision')
% 
% subplot(2,1,2)
% contourf(grid_m,grid_s,pol.adj_size');
% colorbar
% grid on
% xlabel('$m$')
% ylabel('$s$')
% title('Adjustment size')
% sgtitle('Optimal adjustment policy')
% 
% figure(2)
% subplot(2,1,1)
% contourf(grid_m,grid_s(2:end),pol.cash(:,2:end)',[0 1]);
% colormap([0.55 0.55 0.55; 0.85 0.85 0.85])
% hold on
% plot(grid_m(2:end),grid_m(2:end),'r')
% xlabel('$m$')
% ylabel('$s$')
% grid on
% title('Cash vs card decision')
% subplot(2,1,2)
% contourf(grid_m,grid_s(2:end),pol.cash(:,2:end)',[0 1]);
% colormap([0.55 0.55 0.55; 0.85 0.85 0.85])
% hold on
% plot(grid_m(2:end),grid_m(2:end),'r')
% xlabel('$m$')
% ylabel('$s$')
% ylim([2 150])
% xlim([2 150])
% grid on
% title('Cash vs card decision (detail)')

figure('Renderer', 'painters', 'Position', [400 300 900 400])
ax(1)=subplot(1,2,1);
contourf(grid_m,grid_s,pol.adj_size');
colormap(ax(1),gray)
colorbar('Location','southoutside')
grid on
xlabel('$\tilde{m}$')
ylabel('$\tilde{s}$','Rotation',0,'Position',[-40 125 2])
title('Voluntary adjustment policy function $a(\tilde{m},\tilde{s})$')
sgtitle('Policy functions for estimated parameter values $\Theta^{\ast}$')
annotation('textbox',[.09 .1 .5 .01],'String','Optimal post-forced-adjustment cash holding $m^{\ast}=43.96$.','interpreter', 'latex','EdgeColor','none','FontSize',11)
annotation('textbox', [0.61, 0.5, 0.2, 0.2], 'String', "Cashless",'EdgeColor','none','interpreter', 'latex','FontSize',15)
annotation('textbox', [0.78, 0.17, 0.2, 0.2], 'String', "Cash",'EdgeColor','none','interpreter', 'latex','FontSize',15)
ax(2)=subplot(1,2,2);
contourf(grid_m,grid_s(2:end),pol.cash(:,2:end)',[0 1]);
colormap(ax(2),[0.55 0.55 0.55; 0.85 0.85 0.85])
hold on
plot(grid_m(2:end),grid_m(2:end),'r')
xlabel('$m$')
ylabel('$s$','Rotation',0,'Position',[-20 75 2])
ylim([2 150])
xlim([2 150])
grid on
title('Payment choice policy function $p(m,s)$')



% figure(5)
% subplot(2,1,1)
% histogram(money.afternoon(:,22:end))
% subplot(2,1,2)
% histogram(nonzeros(money.adjust(:,22:end)))
% title('Histogram of money holdings')

% [dens_cash,cash_vals]=ksdensity(nonzeros(payment.cash(:)),'Support',[0 251]);
% [dens_cashless,cashless_vals]=ksdensity(nonzeros(payment.cashless(:)),'Support', [0 251]);
% [dens_all,all_vals]=ksdensity(nonzeros(payment.all(:)),'Support',[0 251]);
% 
% 
% figure('Renderer', 'painters', 'Position', [400 300 700 400])
% ax=gca;
% h1=plot(cash_vals,dens_cash,'Color',[0 0.7 0.1]);
% h1.Color(4)=0.6;
% hold on
% h2=plot(cashless_vals,dens_cashless,'Color',[0 0.1 0.7]);
% h2.Color(4)=0.6;
% hold on
% h3=plot(all_vals,dens_all,'Color',[1 0 0],'Linestyle','--','LineWidth',0.5);
% ylim([-0.01 0.17])
% xlim([-5    250])
% xlabel('Payment size')
% ylabel('Density')
% grid on
% ax.GridLineStyle=':';
% ax.GridAlpha=0.35;
% title('Kernel density estimate of simulated payment size distribution')
% legend('Cash','Cashless','All','Location','southoutside')
% legend boxoff
% box(ax,'off')