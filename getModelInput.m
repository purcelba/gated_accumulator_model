function data = getModelInput(data_set,data_directory)
%   Given a requested dataset, load the relevant behavioral and visual neuron spiking data 
%   and concatenate across sessions.
%
%   Input:
%       data_set, str, must be one of the following to indicate which set of data to use
%           - 'F'       monkey F
%           - 'L'       monkey L
%           - 'MM'      monkey M, motion discrimination task
%           - 'MC'      monkey M, color discrimination task
%           - 'pooled'  all data pooled
%       data_directory, str, path where the data files are saved
%   Ouput:
%       data, struct, with the following fields concatenated for requested data set
%           - Spikes, array, raster plot (1s = spikes) trials x time (ms) x neurons/sessions
%                            time is in units of ms starting 500ms before stimulus onset
%           - RT vector, array, response times (in ms) trials x neuron/sessionss
%           - Target, array, stimulus onset (in ms) trials x neurons/sessions
%           - Saccade, array, saccade latency (in ms) trials x neurons/sessions
%

%get list of neurons for the requested data set
if strcmp(data_set,'F')
    mletter_list = {'F'};
    neuron_list = {[1,6,7,10,11,12,13,14,22,24,29]};
elseif strcmp(data_set,'L')
    mletter_list = {'L'};
    neuron_list = {[1,2,4]};
elseif strcmp(data_set,'MM')
    mletter_list = {'MM'};
    neuron_list = {[1,2,3,5,6,7]};
elseif strcmp(data_set,'MC')
    mletter_list = {'MC'};
    neuron_list = {[1,2,3,4]};
elseif strcmp(data_set,'pooled')
    mletter_list = {'F','L','MM','MC'};
    neuron_list = {[1,6,7,10,11,12,13,14,22,24,29],[1,2,4],[1,2,3,5,6,7],[1,2,3,4]};
end

%initialize
data = struct('Spikes',[],'Target',[],'Saccade',[],'TrialStart',[],'T_in_EC',[],'T_in_HC',[],'D_in_EC',[],'D_in_HC',[],'maxSDFactivity',[]);
for m = 1:length(mletter_list)
    for n=neuron_list{m}
        %load data
        load(sprintf('%s/CELL_%s%d.mat',data_directory,mletter_list{m},n))
        %sort
        data.Spikes=cat(3,data.Spikes,Spikes);
        data.Target=cat(2,data.Target,Target);
        data.Saccade=cat(2,data.Saccade,Saccade);
        data.TrialStart=cat(2,data.TrialStart,TrialStart);
        data.T_in_EC=cat(2,data.T_in_EC,T_in_EC);
        data.T_in_HC=cat(2,data.T_in_HC,T_in_HC);
        data.D_in_EC=cat(2,data.D_in_EC,D_in_EC);
        data.D_in_HC=cat(2,data.D_in_HC,D_in_HC);
        data.maxSDFactivity=cat(2,data.maxSDFactivity,maxSDFactivity);
    end
end


