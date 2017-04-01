function plotModel(RTpred_cdf,traj_mat,SDF_mat,params,save_plots) 
%PlotModel
%
%   Input:
%       RTpred_cdf, struct, two fields
%           'easy' cumulative predicted RT distribution for easy condition
%           'hard' cumulative predicted RT distribution for hard condition
%       traj_mat, struct, two fields
%           'T' predicted model output (simulated trials x time) for target in RF accumulator
%           'D' predicted model output (simulated trials x time) for distractor in RF accumulator
%       SDF_mat, struct, two fields
%           'Tin' model input for target in RF
%           'Din' model input for distractor in RF
%       params, struct, with the following parameter fields
%           'theta' threshold parameter
%           'g' gating parameter
%       save_plots, logical, if true, save a png copy of all figures.
%   

%observed quantiles for comparison to model (divided at 10th, 30th, 50th, 70th, 90th percentiles)
obs_quantiles.easy=[189.2 205.73 219.34 233 258.11];
obs_quantiles.hard=[220.64 252.28 276.81 305.23 370.34];

%plot predicted CDFs on top of data
figure(1)
axes('position',[.35 .45 .3 .3])
hold on
plot(-500:1000,RTpred_cdf.easy,'color',[.6 0 0],'linewidth',2)
plot(-500:1000,RTpred_cdf.hard,'color',[0 .6 0],'linewidth',2)
plot(obs_quantiles.easy,[.1 .3 .5 .7 .9],'linestyle','none','marker','o','color',[.6 0 0],'markersize',6,'linewidth',2)
plot(obs_quantiles.hard,[.1 .3 .5 .7 .9],'linestyle','none','marker','o','color',[0 .6 0],'markersize',6,'linewidth',2)
legh=legend('Easy - Predicted','Hard - Predicted','Easy - Observed','Hard - Observed');
legend('boxoff')
set(legh,'position',[.75 .55 .1 .1])
xlim([100 500])
ylabel('CDF')
xlabel('Time from array onset (ms)')
if save_plots
    print(gcf,'RT_CDFs.png','-dpng')
end

%plot model input
figure(2)
x_limits=[-200 500];
y_limits=[0 1]; 
subplot('position',[.35 .55 .3 .3])
hold on
plot(-500:1000,mean(SDF_mat.Tin),'color','k','linewidth',3)
plot(-500:1000,mean(SDF_mat.Din),'color','k','linewidth',1)
plot([-500 1000],[params.g params.g],'r','linestyle',':')
xlim([x_limits])
ylim([0 1])
set(gca,'xticklabel',[])
ylabel('Normalized Activity')
title('Average Model Input')
legh=legend('Target in RF','Distractor in RF','g');
legend('boxoff')
set(legh,'position',[.8,.8,.1,.1])
%plot example output (model activation)
subplot('position',[.35 .15 .3 .3])
hold on
nSims = size(traj_mat.T,1);
n=round(nSims/2);
plot(-500:1000,mean(traj_mat.T(n:n+10,:)),'color','k','linewidth',3)
plot(-500:1000,mean(traj_mat.D(n:n+10,:)),'color','k','linewidth',1)
plot([x_limits],[params.theta params.theta],'--','color','k','linewidth',2)
ylim([0 params.theta])
xlim(x_limits)
ylabel('Model Activation')
xlabel('Time from array onset (ms)')
title('Sample Model Output')
if save_plots
    print(gcf,'model_input_output.png','-dpng')
end
    
    
    
    
   