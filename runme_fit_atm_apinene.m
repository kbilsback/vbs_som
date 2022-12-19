clc; clear; close all;

% Setup
sp = 'napthalene';
noxcase = 'hi';
wlcase = 'wl0';
vbs.cstar = 10.^[-2, -1, 0, 1];

voc.initial = 1e-6; % initial VOC concentration, ppm (confirm with Charles)

%a-pinene
voc.mw = 136.23; % molecular weight, g/mole
voc.koh = 5.23E-11; % reaction rate constant with OH, cm3/molecules/s

% limonene
%voc.mw = 136.24; % molecular weight, g/mole
%voc.koh = 1.64E-10; % reaction rate constant with OH, cm3/molecules/s

%b-caryophyllene
%voc.mw = 204.36; % molecular weight, g/mole
%voc.koh = 1.97E-10; % reaction rate constant with OH, cm3/molecules/s

% benzene
%voc.mw = 78.11; % molecular weight, g/mole
%voc.koh = 1.22E-12; % reaction rate constant with OH, cm3/molecules/s

% toluene
%voc.mw = 92.14; % molecular weight, g/mole
%voc.koh = 5.62E-12; % reaction rate constant with OH, cm3/molecules/s

% m-xylene
%voc.mw = 106.16; % molecular weight, g/mole
%voc.koh = 2.31E-11; % reaction rate constant with OH, cm3/molecules/s

% naphthalene
%voc.mw = 128.2; % molecular weight, g/mole
%voc.koh = 2.44E-11; % reaction rate constant with OH, cm3/molecules/s


% Other
env.OHconc = 1.5e6; % OH concentration, molecules/cm3
env.COA = [0.1 1 10]; % OA mass concentration, µg/m3
env.Pres = 101325; % pressure, Pa
env.Temp = 298; % temperature, K
env.Rgas = 8.314; % universal gas constant, gm2/s2/K/mole
dt = 100; % time step, s
maxtime = 7*24*3600; % total time, s
%env.time = 601:dt:maxtime;
env.time = 0*24*3600:dt:maxtime; %only fit 48 hours on

% VBS
%vbs.cstar = 10.^[-2, -1, 0, 1]; % saturation concentration, µg/m3
vbs.nbins = length(vbs.cstar); % number of VBS bins
%
% Measured SOA mass yield data
load(strcat(sp,'_pt1ug_',noxcase,'_',wlcase,'.mat'));
ts1 = timeseries(meas.yield,meas.time*3600); % set data as a timeseries for 48 hours
%ts1 = timeseries(ts1.Data(find(ts1.Time >= 5*24*3600)),ts1.Time(find(ts1.Time >= 5*24*3600))- 5*24*3600);
ts2 = resample(ts1,env.time); % resample at a higher time resolution
meas2(1).time = env.time;
meas2(1).yield = ts2.data;

load(strcat(sp,'_1ug_',noxcase,'_',wlcase,'.mat'));
ts1 = timeseries(meas.yield,meas.time*3600); % set data as a timeseries
%ts1 = timeseries(ts1.Data(find(ts1.Time >= 5*24*3600)),ts1.Time(find(ts1.Time >= 5*24*3600))- 5*24*3600);
ts2 = resample(ts1,env.time); % resample at a higher time resolution
meas2(2).time = env.time;
meas2(2).yield = ts2.data;

load(strcat(sp,'_10ug_',noxcase,'_',wlcase,'.mat'));
ts1 = timeseries(meas.yield,meas.time*3600); % set data as a timeseries
%ts1 = timeseries(ts1.Data(find(ts1.Time >= 5*24*3600)),ts1.Time(find(ts1.Time >= 5*24*3600))- 5*24*3600);
ts2 = resample(ts1,env.time); % resample at a higher time resolution
meas2(3).time = env.time;
meas2(3).yield = ts2.data;

% save full time series to .csv (comment out for fitting)
time = transpose(env.time/3600);
st_pt1 = meas2(1).yield;
st_1 = meas2(2).yield;
st_10 = meas2(3).yield;

m_save = [time, st_pt1, st_1, st_10];
csvwrite(strcat(sp,'_',noxcase,'_',wlcase,'_no_fit.csv'), m_save);
'Done'
%%
% Initial guesses
x0 = [0.001+zeros(1,vbs.nbins) log10(voc.koh) log10(1)];

%% Perform fits
f = @(x)fityield(x,meas2,voc,vbs,env);
options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',20000,'FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',4000);
%[x, resnorm] = lsqnonlin(f,x0,[zeros(1,vbs.nbins) log10(1e-13) log10(0.1)],[ones(1,vbs.nbins) log10(1e-9) log10(10)],options);
[x, resnorm] = lsqnonlin(f,x0,[zeros(1,vbs.nbins) log10(1e-100) log10(1e6)],[ones(1,vbs.nbins)*2 log10(1e-100) log10(1e6)],options); % no aging

%% Check fits
% env.COA = [0.1 1 10 100]; % OA mass concentration, µg/m3
[yield, Cgas, Cpar] = calcyield(x,voc,vbs,env);

%% Plots
% Plot SOM and VBS yields
figure(1); title('a-pinene');
subplot(1,2,1); hold on;
plot(env.time/3600,meas2(1).yield,'-r','LineWidth',3);
plot(env.time/3600,meas2(2).yield,'-b','LineWidth',3);
plot(env.time/3600,meas2(3).yield,'-g','LineWidth',3);
plot(env.time/3600,yield,'-.','LineWidth',3);
xlabel('Time (hours)'); ylabel('SOA Mass Yield');
legend('SOM-TOMAS (0.1 ug/m^3)','SOM-TOMAS (1 ug/m^3)','SOM-TOMAS (10 ug/m^3)',...
    'VBS_{SOM} (0.1 ug/m^3)','VBS_{SOM} (1 ug/m^3)','VBS_{SOM} (10 ug/m^3)','location','northwest');
set(gca,'tickdir','out'); set(gca,'LineWidth',1.5,'TickLength',[0.05 0.05]); set(gcf,'color','w'); set(gca,'fontsize', 16);
box off; axis square;

% Plot final fit values
subplot(1,2,2); hold on;
bar(log10(vbs.cstar),x(1:vbs.nbins));
% bar(log10(vbs.cstar),x(1:vbs.nbins)/sum(x(1:vbs.nbins)));
xlabel('log_{10}C^{*}'); ylabel('Mass Yield');
title1 = sprintf('k_{aging}=%7.2e, p_{frag}=%5.2f',10^(x(vbs.nbins+1)),10^(x(vbs.nbins+2)));
title(title1);
set(gca,'tickdir','out'); set(gca,'LineWidth',1.5,'TickLength',[0.05 0.05]); set(gcf,'color','w'); set(gca,'fontsize', 16);
box off; axis square;

x(1:vbs.nbins)

%% save to .csv
time = transpose(env.time/3600);
st_pt1 = meas2(1).yield;
st_1 = meas2(2).yield;
st_10 = meas2(3).yield;
vbs_pt1 = transpose(yield(1,:));
vbs_1 = transpose(yield(2,:));
vbs_10 = transpose(yield(3,:));

m_save = [time, st_pt1, st_1, st_10, vbs_pt1, vbs_1, vbs_10];
csvwrite(strcat(sp,'_',noxcase,'_',wlcase,'.csv'), m_save);

% % Plot gas and particle concentrations
% figure(3);
% subplot(1,2,1); hold on;
% plot(meas.time,meas.tgas,'--k','LineWidth',2);
% plot(env.time/3600,sum(Cgas,2),'-b','LineWidth',3);
% xlim([0 maxtime(end)/3600]);
% xlabel('Time (hours)'); ylabel('Concentrations (\mug m^{-3})');
% legend('SOM-TOMAS','VBS_{SOM}','location','southwest'); title('Gas');
% set(gca,'tickdir','out'); set(gca,'LineWidth',1.5,'TickLength',[0.05 0.05]); set(gcf,'color','w'); set(gca,'fontsize', 16);
% box off; axis square;
% subplot(1,2,2); hold on;
% plot(meas.time,meas.tpar,'--k','LineWidth',2);
% plot(env.time/3600,sum(Cpar,2),'-r','LineWidth',3);
% xlim([0 maxtime(end)/3600]);
% xlabel('Time (hours)'); ylabel('Concentrations (\mug m^{-3})');
% legend('SOM-TOMAS','VBS_{SOM}','location','southwest'); title('Particle');
% set(gca,'tickdir','out'); set(gca,'LineWidth',1.5,'TickLength',[0.05 0.05]); set(gcf,'color','w'); set(gca,'fontsize', 16);
% box off; axis square;
% 
% % Plot gas/particle ratio
% figure(4); hold on;
% plot(meas.time,meas.tpar./meas.tgas,'--k','LineWidth',2);
% plot(env.time/3600,sum(Cpar,2)./sum(Cgas,2),'-r','LineWidth',3);
% xlim([0 maxtime(end)/3600]);
% xlabel('Time (hours)'); ylabel('Particle:Gas');
% legend('SOM-TOMAS','VBS_{SOM}');
% set(gca,'tickdir','out'); set(gca,'LineWidth',1.5,'TickLength',[0.05 0.05]); set(gcf,'color','w'); set(gca,'fontsize', 16);
% box off; axis square;

% % Plot volatility distributions
% [ntime nsp] = size(meas.pargrid);
% vbsforplot = 10.^[-1:1:10];
% nvbsforplot = length(vbsforplot);
% Cvolgp = zeros(ntime,nvbsforplot,2);
% Cvol = zeros(ntime,nvbsforplot);
% for j=1:nsp
%     k = round(log10(meas.cstar(j))) - log10(vbsforplot(1)) + 1;
%     if k>nvbsforplot
%         k = nvbsforplot;
%     elseif k<1
%         k = 1;
%     end
%     Cvolgp(:,k,1) = Cvolgp(:,k,1) + meas.pargrid(:,j);
%     Cvolgp(:,k,2) = Cvolgp(:,k,2) + meas.gasgrid(:,j);
%     Cvol(:,k) = Cvol(:,k) + meas.pargrid(:,j) + meas.gasgrid(:,j);
% end
% figure(5); hold on;
% i1 = min(find(meas.time>maxtime/3600)); 
% bar(log10(vbsforplot),squeeze(Cvol(i1,:)),'k')
% plot(log10(vbs.cstar),Cgas(end,:)+Cpar(end,:),'-r','LineWidth',3);
% xlabel('log_{10}C^*'); ylabel('G+P (\mug m^{-3})');
% legend('SOM-TOMAS','VBS_{SOM}','location','northwest');
% set(gca,'tickdir','out'); set(gca,'LineWidth',1.5,'TickLength',[0.05 0.05]); set(gcf,'color','w'); set(gca,'fontsize', 16);
% box off; axis square;
% 
% % Plot absolute and normalized contribution of volatility bins
% figure(6); 
% subplot(1,2,1);
% area(meas.time,Cvol);
% xlabel('Time (hours)'); ylabel('G+P (\mug m^{-3})');
% set(gca,'tickdir','out'); set(gca,'LineWidth',1.5,'TickLength',[0.05 0.05]); 
% set(gcf,'color','w'); set(gca,'fontsize', 16);
% box off; axis square;
% subplot(1,2,2);
% for i=1:ntime
%     Cvolnorm(i,:) = Cvol(i,:)/sum(Cvol(i,:));
% end
% area(meas.time,Cvolnorm);
% xlabel('Time (hours)'); ylabel('G+P (norm.)');
% set(gca,'tickdir','out'); set(gca,'LineWidth',1.5,'TickLength',[0.05 0.05]); 
% set(gcf,'color','w'); set(gca,'fontsize', 16);
% box off; axis square;
% 
% %% Simulate dilution from Coa=10 µg/m3 to 1 µg/m3
% dr = 1:1:100; % dilution ratio
% for i=1:length(dr)
%     coa(i) = 10/dr(i);
%     som_cpardr(i) = 0;
%     vbs_cpardr(i) = 0;
%     cpar_dr(i) = 0;
%     for j=1:nsp
%         zeta = (1 + meas.cstar(j)/coa(i))^(-1);
%         som_cpardr(i) = som_cpardr(i) + zeta * (meas.pargrid(i1,j)+meas.gasgrid(i1,j));
%     end
%     for j=1:vbs.nbins
%         zeta = (1 + vbs.cstar(j)/coa(i))^(-1);
%         vbs_cpardr(i) = vbs_cpardr(i) + zeta * (Cpar(end,j)+Cgas(end,j));
%     end
% end
% figure(7); hold on;
% plot(coa,som_cpardr,'--k','LineWidth',3);
% plot(coa,vbs_cpardr,'-r','LineWidth',3);
% % xlim([0 maxtime(end)/3600]);
% xlabel('C_{OA} (\mug m^{-3})'); ylabel('Particle concentration (\mug m^{-3})');
% legend('SOM-TOMAS','VBS_{SOM}');
% set(gca,'tickdir','out'); set(gca,'LineWidth',1.5,'TickLength',[0.05 0.05]); set(gcf,'color','w'); set(gca,'fontsize', 16);
% box off; axis square;