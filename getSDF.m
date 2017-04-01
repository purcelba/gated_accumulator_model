function SDF = getSDF(Spike, Align_Time, Plot_Time, triallist, TrialStart_, maxSDFactivity_recip,normalize)
% creates spike-density function from spike times

%Spike = input all
%align = input all
%trialstart = input only selected trials
%max activity = input only selected trials
%rt vector = input all


% In:
% (1) the spike times [Spike], 
% (2) the event, on which the histogram
% will be aligned [Align_event], 
% (3) the time before and after the align
% event, that should be plotted [Plot_Time], 
% (4) the list of used trials
% [triallist], 
% (5) beginning of each trial [TrialStart_], 
% (6) the reciprocal
% of the maximum activity for the cell in each trial - normally will be the
% same for each trial unless different cells are being used, 
% (7) the onset of
% the event after/before which spikes are truncated [TruncEvent],
%
% Out:
% SDF: the spike density function


%initialize
window_end=0;
time_before_saccade=0-Plot_Time(1);
EmptyTrials=0;%Trials in which there are no spikes in the interval
ListTimes = []; times = []; Kernel=[]; poslist=[]; temp=[]; pl=[];
%align time event
if(isequal(Align_Time,TrialStart_))
   Align_Time(1:length(Align_Time)) = 0;
end
%Pre-Time & Post-Time
Pre_Time = Plot_Time(1)-100; Post_Time = Plot_Time(2);
BinCenters = [round(Pre_Time):round(Post_Time)]';
%Create raw histogram
hist_times_matrix = zeros(length(triallist),length(BinCenters));
%Add up spikes
for pl=1:length(triallist)
   trl = triallist(pl);
   times = nonzeros(Spike(trl,:));
   if(isempty(nonzeros(times)))
       EmptyTrials = EmptyTrials+1;
   else
       times = nonzeros(times)-Align_Time(triallist(pl),1);
       times = times(times >= Pre_Time & times <= Post_Time);
       %ListTimes = [[ListTimes];[times]]
       hist_times = hist(times,BinCenters);
       if normalize==1
            hist_times_matrix(pl,:) = maxSDFactivity_recip(pl)*hist_times;
       elseif normalize==0
            hist_times_matrix(pl,:) = hist_times;
       end
   end
end
% 3. Convolve with exponential Kernel
Growth=1; Decay=20;
Half_BW=round(Decay*8);
BinSize=(Half_BW*2)+1;
Kernel=[0:Half_BW];
Half_Kernel=(1-(exp(-(Kernel./Growth)))).*(exp(-(Kernel./Decay)));
Half_Kernel=Half_Kernel./sum(Half_Kernel);
Kernel(1:Half_BW)=0;
Kernel(Half_BW+1:BinSize)=Half_Kernel;
Kernel=Kernel.*1000;
Kernel=Kernel';
if ~isempty(hist_times_matrix)
    trialcount=zeros(1,length(hist_times_matrix(1,:)));
    for pl = 1:length(hist_times_matrix(1,:))
       trialcount(pl)=sum(~isnan(hist_times_matrix(:,pl)));
    end
    temp_hist = nansum(hist_times_matrix,1)./trialcount;
    Hist_raw=temp_hist(:);
    SDF=convn(Hist_raw,Kernel,'same');
    SDF(1:100,:,:)=[];
else
   SDF=nan(length(Plot_Time(1):Plot_Time(2)),1);
end
