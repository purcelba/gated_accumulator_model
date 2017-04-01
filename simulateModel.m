function simulateModel(params, options)
% function simulateModel(params, options)
%
% Code to run versions of the models presented in Purcell, Heitz, Cohen, 
% Logan, Schall, & Palmeri (2010).  Neurally Constrained Modeling of
% Perceptual Decision Making.  Psychological Review.  
% 
% By default, the code uses the best fit parameters of the Gated Race Model 
% to the pooled data set (Table 2) and compares it to the pooled observed behavior.
% Plots (1) Observed and Predicted response times (RTs) (2) Average model inputs 
% and (3) Average model activations.
% 
% Input:
%   params, struct with fields indicating model parameters
%       N, integer, indicating number of input neurons to sample
%       theta, float, indicating response threshold
%       tballistic, float, indicating motor response time
%       k, float, leakage parameter 
%       g, float, gating parameter
%       beta, lateral inhibition
%       u, feed-forward inhibition
%   options, struct with fields indicating simulation options
%       nSims, integer, number of simulations
%       diff_start, integer, when to start the integration process. 
%       data_directory str, directory in which the data files are located.
%       save_plots, logical, if true, save a png copy of all plots.
% 
% Code written by: Braden Purcell (braden@nyu.edu; Vanderbilt University, New York University)
% Original Date: April 5, 2010
% Revised: March 30, 2017
% 
clear all; close all;

%Default parameters
if nargin<1
    %best-fitting values for Gated Race Model (Table 2), pooled data.
    params.N=22;               %sample size
    params.theta=17.18;        %threshold
    params.tballistic=15.03;   %post-decision time
    params.k=0.0003;           %leakage
    params.g=0.5782;           %gate
    params.beta=0;             %lateral inhibition (0 for race, 0 for diffusion, >0 for competitive)
    params.u=0;                %feed-forward inhibition (0 for race, 1 for diffusion, 0 for competitive)
end
if nargin<2
    %default simulation options
    options.nSims = 1000;       %5000 used for all simulations in paper
    options.diff_start = 200;   %200
    options.data_set = 'pooled';    
    options.data_directory = 'Data'; %by default, the data are stored in the Data subdirectory.
    options.save_plots = true;
end  
    
%default seed
rand('twister',5489)
randn('state',1)
    
%load dataset
fprintf('Loading data...')
data = getModelInput(options.data_set,options.data_directory);
num_cells = size(data.Spikes,3);
fprintf('Loaded.\n')

%Loop through easy and hard conditions (2=easy, 3=hard)
for easy_hard=[2 3]
    %display
    if easy_hard == 2
        fprintf('Simulating easy condition.\n')
    elseif easy_hard == 3
        fprintf('Simulating hard condition.\n')
    end
     %random sample cells for each simulation
    random_cell_samples = nan(params.N,options.nSims);
    for i = 1:params.N
        random_cell_samples(i,:) = randsample(num_cells,options.nSims,true);
    end
    %replace each cell number with a T-in and D-in trial for that cell.
    random_Tin_samples = zeros(size(random_cell_samples));
    random_Din_samples = zeros(size(random_cell_samples));
    for ce = 1:num_cells
        ki = find(random_cell_samples==ce);
        if easy_hard==2
            random_Tin_samples(ki) = randsample(nonzeros(data.T_in_EC(:,ce)),length(ki),true);
            random_Din_samples(ki) = randsample(nonzeros(data.D_in_EC(:,ce)),length(ki),true);
        elseif easy_hard==3
            random_Tin_samples(ki) = randsample(nonzeros(data.T_in_HC(:,ce)),length(ki),true);
            random_Din_samples(ki) = randsample(nonzeros(data.D_in_HC(:,ce)),length(ki),true);
        end
    end
    
    %Convert spike times to spike rate for each trial (spike density functions, SDFs)
    fprintf('\tGenerating SDFs...')
        %specify some options
    time_index = -500:1000;
    Plot_Time = [-500 1000];
    normalize = 1;
        %initialize variables to convert spikes to spike density functions.
    SDF_mat.Tin=zeros(options.nSims,length(time_index));
    SDF_mat.Din=zeros(options.nSims,length(time_index));
    Spikes_Tin = zeros(params.N,size(data.Spikes,2));
    Spikes_Din = zeros(params.N,size(data.Spikes,2));
    Target_Tin = zeros(params.N,1);
    Target_Din = zeros(params.N,1);
    TrialStart_Tin = zeros(params.N,1);
    TrialStart_Din = zeros(params.N,1);
    maxSDFactivity_reciprocal = zeros(params.N,1);
        %begin generating SDFs for model input 
    for n = 1:options.nSims
        %initialize
        cell_vector = random_cell_samples(:,n);
        Tin_trials = random_Tin_samples(:,n);
        Din_trials = random_Din_samples(:,n);
        %get data from each simulation
        for i = 1:params.N
            %save spike data for current simulation in trial x spikes matrix
            Spikes_Tin(i,:) = data.Spikes(Tin_trials(i),:,cell_vector(i));
            Spikes_Din(i,:) = data.Spikes(Din_trials(i),:,cell_vector(i));
            %save Target_
            Target_Tin(i) = data.Target(Tin_trials(i),cell_vector(i));%align on target onset
            Target_Din(i) = data.Target(Din_trials(i),cell_vector(i));%align on target onset
            %save TrialStart_
            TrialStart_Tin(i) = data.TrialStart(Tin_trials(i),cell_vector(i));
            TrialStart_Din(i) = data.TrialStart(Din_trials(i),cell_vector(i));
            %save max activity for each cell
            maxSDFactivity_reciprocal(i) = 1/data.maxSDFactivity(cell_vector(i));
        end
        %generate SDFs.  save in matrix
        SDF_mat.Tin(n,:)=getSDF(Spikes_Tin, Target_Tin, Plot_Time, 1:params.N, TrialStart_Tin,maxSDFactivity_reciprocal,normalize);
        SDF_mat.Din(n,:)=getSDF(Spikes_Din, Target_Din, Plot_Time, 1:params.N, TrialStart_Din,maxSDFactivity_reciprocal,normalize);
    end
    fprintf('DONE.\n')
    
    %Simulate model
    fprintf('\tRunning model...')
        %generate Gaussian noise (SD=0.2)
    noiseT=randn([options.nSims,length(Plot_Time(1):Plot_Time(2))])*0.2;
    noiseD=randn([options.nSims,length(Plot_Time(1):Plot_Time(2))])*0.2;  
        %initialize for model runs
    traj_matrix_T = zeros(options.nSims,length(time_index));
    traj_matrix_D = zeros(options.nSims,length(time_index));
    RTpred=zeros(1,options.nSims);
    correct_vector=zeros(1,options.nSims);
        %run model. save predicted behavior and model trajectories.
    Tinput = zeros(1,size(SDF_mat.Tin,2));
    Dinput = zeros(1,size(SDF_mat.Din,2));
    for n = 1:options.nSims
        Tinput(options.diff_start:end) = SDF_mat.Tin(n,options.diff_start:end);
        Dinput(options.diff_start:end) = SDF_mat.Din(n,options.diff_start:end);
        [RTpred(n),correct_vector(n),Tunit,Dunit]=model(Tinput,Dinput,options.diff_start,time_index,params.theta,params.tballistic,params.k,params.g,params.beta,noiseT(n,:),noiseD(n,:),params.u);
        traj_matrix_T(n,1:length(Tunit)) = Tunit;
        traj_matrix_D(n,1:length(Dunit)) = Dunit;
    end%simulation loop
	fprintf('DONE\n')
       
    %Analyze output
    fprintf('\tAnalyzing results...')
        %predicted correct RT distributions
    RTpred_correct=RTpred(correct_vector==1);
    traj_mat.T=traj_matrix_T(correct_vector==1,:);
    traj_mat.D=traj_matrix_D(correct_vector==1,:);
        %sort
    [~,index]=sort(RTpred_correct);
    traj_mat.T=traj_mat.T(index,:);
    traj_mat.D=traj_mat.D(index,:);    
        %Get predicted CDFs
    RTpred_hist=hist(RTpred_correct,-500:1000);
    if easy_hard==2
        RTpred_cdf.easy=cumsum(RTpred_hist)/length(RTpred_correct);
    elseif easy_hard==3
        RTpred_cdf.hard=cumsum(RTpred_hist)/length(RTpred_correct);
    end
   fprintf('DONE.\n')
    
end%easy_hard loop

%Plot
%   (1) Observed/Predicted CDFs for easy and hard condition
%   (2) Average model input
%   (3) Average sample of 10 trajectories (sampled at median)
fprintf('Generating plots...')
plotModel(RTpred_cdf,traj_mat,SDF_mat,params,options.save_plots)
fprintf('Plots generated.\n')
fprintf('Finished.\n')

             