%%
%*****plot Figure 3 **********
%
%
clc;clear all
a = load('256个样本点\GMC.mat');    %GMC
b = load('256个样本点\CCMF.mat');  %CCMF
c = load('256个样本点\WCSD.mat');  %WCSD
d = load('256个样本点\CMF.mat');    %CMF
e = load('256个样本点\NAC.mat');    %NAC

figure(1);
plot( a.pf(:,2),a.pd(:,2),'o-b',...
        b.pf_2(:,2),b.pd_2(:,2),'*-b',...
        c.pf_2(:,3),c.pd_2(:,3),'s-b',...
        d.pf_2(:,2),d.pd_2(:,2),'+-b',...
        e.pf,e.pd,'v-b','LineWidth',1.5);
xlabel('False alarm probability,{\it P_f} , GSNR = -10dB'); 
ylabel('Detection probability,{\it P_d }');
leg = legend('GMC','CCMF','WCSD','CMF','NAC');
set(leg,'Location','SouthEast');
set(gcf, 'position', [400 400 400 300]);
grid on
hold on