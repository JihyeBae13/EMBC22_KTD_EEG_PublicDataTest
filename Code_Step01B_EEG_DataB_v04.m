%% Clearing Workspace, Command Window and Closing All the Figures
clear; clc; close all

%% Setting Up Gloabl Variables (including directories)- Loading of the Dataset 
%dataPath: It is the path for folder where all the raw dataset provided for
%all 9 subjects.
%resultPath: This is the folder where you want your results (feature data and images) to be stored.
%subjectName = name of the Subject that you want to visualize.

dataPath = 'C:\Users\jbae2\Desktop\Jihye Bae\20211202 467C 5820\UK\Lab\2020_Bhoj02\Papers\EMBC01_KTDEEG\Codes\BCIcompetitionIV-2a';
resultPath = 'C:\Users\jbae2\Desktop\Jihye Bae\20211202 467C 5820\UK\Lab\2020_Bhoj02\Papers\EMBC01_KTDEEG\Codes\Results';
dataName = 'A09';
%'A01' or 'A02' or 'A03' or 'A04' or 'A05'
%'A06' or 'A07' or 'A08' or 'A09'

cd(resultPath)
if exist(dataName,'dir')~= 7
    mkdir(dataName)
end

resultPath = [resultPath '\' dataName];
cd(dataPath)
ChannelNames = ["Fz";...
    "FC3";"FC1";"FCz";"FC2";"FC4";...
    "C5";"C3";"C1";"Cz";"C2";"C4";"C6";...
    "CP3";"CP1";"CPz";"CP2";"CP4"; ...
    "P1";"Pz";"P2"; ...
    "POz"];%Channels Name.

%% Loading Selected Subject Dataset

AE = load([dataName 'E.mat']); AT = load([dataName 'T.mat']); % Loading Both Training and Evaluation Data. 
RawEEG = [AT.data(end-5:end),AE.data(end-5:end)]; % EEG data. adding two sessions of data for respective subject. Some subject does not have EOG blocks, so load from (end-5:end)

nrun = length(RawEEG);
itr = 1;
for irun = 1: nrun
    if irun == 1
        [tlength, nch] = size(RawEEG{1,irun}.X); % tlength=time lenghth, nch=number of channels.
        fs = RawEEG{1,1}.fs; % Sampling frequency.
        classID = RawEEG{1,irun}.y;
        artifacts = RawEEG{1,irun}.artifacts;
    else
        [tlength_ts, nch_ts] = size(RawEEG{1,irun}.X);
        fs_ts = RawEEG{1,irun}.fs;
        if tlength ~= tlength_ts
            fprintf('==========Error============\n')
            fprintf('Inconsistant Time Length in %d Session \n', irun)
        end
        
        if nch ~= nch_ts
            fprintf('==========Error============\n')
            fprintf('Inconsistant Channel Number in %d Session \n', irun)
        end
        
        if fs ~= fs_ts
            fprintf('==========Error============\n')
            fprintf('Inconsistant Sampling Frequency in %d Session \n', irun)
        end
        classID = [classID; RawEEG{1,irun}.y];
        artifacts = [artifacts; RawEEG{1,irun}.artifacts];
    end
end

%% Visualization of the Trials Distribution and sequence of the Training and Evaluation Datasets

% This Section shows how the sessions were recorded in terms of motor
% imagery task distribution

h = figure;
subplot(2,1,1)
hold on
plot(classID,'k')
plot(classID == 1,'r*') %Class 1: Left hand
plot((classID == 2)*2,'g*')%Class 2: Right hand
plot((classID == 3)*3,'b*')%Class 3: Both Feet
plot((classID == 4)*4,'m*')%Class 4: Tongue
yticks([1 2 3 4]);
yticklabels({'Left Hand','Right Hand' ,'Both Feet' ,'Tongue'});
ylim([1 4])
xlabel('Trial Index')
ylabel('Class')
title(['Class Information ' dataName])
set(gca,'fontsize', 18);

subplot(2,1,2)
hold on
plot(classID,'k')
plot(classID == 1,'r*')
plot((classID == 2)*2,'g*')
plot((classID == 3)*3,'b*')
plot((classID == 4)*4,'m*')
ylim([1 4])
xlim([100 200])
yticks([1 2 3 4]);
yticklabels({'Left Hand','Right Hand' ,'Both Feet' ,'Tongue'});
xlabel('Trial Index')
ylabel('Class')
title('(Zoomed In)')
set(gca,'fontsize', 18);

%%
%%%%%%%%%%%%%%%%%%%%% DATA SEGMENTATION %%%%%%%%%%%%%%%%%%%%%%%%%
tstart = 2; % Define start time to extract EEG. Relative to the strat of a trial.
tend = 2.85; % Define end time to extract EEG. Relative to the strat of trial.
ntpoints = round(tend*fs)-round(tstart*fs);%Number of samples in defined time.
ntrialall = length(classID);
trialEEG = nan(ntpoints,nch,ntrialall);
for irun = 1: nrun
    nrunTrial = length(RawEEG{1,irun}.y);
    classTrialIndex = RawEEG{1,irun}.trial; % start trial index
    for irunTrial = 1: nrunTrial
        trialEEG(:,:,itr) = RawEEG{1,irun}.X(classTrialIndex(irunTrial)+round(tstart*fs):classTrialIndex(irunTrial)+round(tend*fs)-1,:);
        itr = itr + 1;
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Feature1: RAW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classID = classID(not(artifacts),1);
trialEEG = trialEEG(:,:,not(artifacts));
nclass = max(classID); %Number of classes.
ntrial = length(classID); % Number of trials
numEOGch = 3; % The last 3 channels are for EOG, known from visual inspection.
RAWfeature = nan(ntrial,ntpoints*(nch-numEOGch));
for itr = 1:ntrial
    for ich = 1: nch-numEOGch
        RAWfeature(itr,ntpoints*(ich-1)+1:ntpoints*ich) = trialEEG(:,ich,itr);
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Feature2: Cartesian FTA, REAL & IMAGINARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trialEEGfft = fft(trialEEG); % FFT of the segmented EEG
lFFT = size(trialEEG,1); % length of FFTf

trialEEGfftImaginary = imag(trialEEGfft); % 2-sided spectrum
trialEEGfftReal = real(trialEEGfft); % 2-sided spectrum
trialEEGfftImaginary1side = 2*trialEEGfftImaginary(1:lFFT/2+1,:,:);
trialEEGfftReal1side = trialEEGfftReal(1:lFFT/2+1,:,:);
trialEEGfftReal1side(2:end,:,:) = 2*trialEEGfftReal1side(2:end,:,:);

nFreqComp = 13;
nFTAfeature = 2*nFreqComp-1;
FTAfeature = nan(ntrial,(nch-numEOGch)*nFTAfeature);
for itr = 1:ntrial     
    for ich = 1: nch-numEOGch
        FTAfeature(itr,nFTAfeature*(ich-1)+1:nFTAfeature*ich) = [trialEEGfftReal1side(1:nFreqComp,ich,itr); ...
            trialEEGfftImaginary1side(2:nFreqComp,ich,itr)]'; 
    end
end

cd(resultPath)
save('Data_EEGfeatureRAW.mat','RAWfeature','classID')
save('Data_EEGfeatureFTA.mat','FTAfeature','classID')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting Features: 1) ERP and 2) Cartesian FTA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ChToPlot = "C3"; % Channel name to plot

ChToPlotIndex = find(ChToPlot == ChannelNames);
if isempty(ChToPlotIndex)
    fprintf('=====================ERROR========================\n')
    fprintf('The entered channel to plot does not exist or is wrong\n')
else
    % Plotting RAW EEG
    trialEEGToPlot = squeeze(trialEEG(:,ChToPlotIndex,:));

    h = figure;
    tplot = tstart:1/fs:tend;
    tplot = tplot';
    hold on;
    plot(tplot,mean(trialEEGToPlot(:,classID==1),2),'r','LineWidth',2);
    plot(tplot,mean(trialEEGToPlot(:,classID==2),2),'b','LineWidth',2);
    plot(tplot,mean(trialEEGToPlot(:,classID==3),2),'g','LineWidth',2);
    plot(tplot,mean(trialEEGToPlot(:,classID==4),2),'m','LineWidth',2);
    plot(tplot,mean(trialEEGToPlot(:,classID==1),2)+std(trialEEGToPlot(:,classID==1),[],2),':r','LineWidth',1.5);
    plot(tplot,mean(trialEEGToPlot(:,classID==1),2)-std(trialEEGToPlot(:,classID==1),[],2),':r','LineWidth',1.5);
    plot(tplot,mean(trialEEGToPlot(:,classID==2),2)+std(trialEEGToPlot(:,classID==2),[],2),':b','LineWidth',1.5);
    plot(tplot,mean(trialEEGToPlot(:,classID==2),2)-std(trialEEGToPlot(:,classID==2),[],2),':b','LineWidth',1.5);
    plot(tplot,mean(trialEEGToPlot(:,classID==3),2)+std(trialEEGToPlot(:,classID==3),[],2),':g','LineWidth',1.5);
    plot(tplot,mean(trialEEGToPlot(:,classID==3),2)-std(trialEEGToPlot(:,classID==3),[],2),':g','LineWidth',1.5);
    plot(tplot,mean(trialEEGToPlot(:,classID==4),2)+std(trialEEGToPlot(:,classID==4),[],2),':m','LineWidth',1.5);
    plot(tplot,mean(trialEEGToPlot(:,classID==4),2)-std(trialEEGToPlot(:,classID==4),[],2),':m','LineWidth',1.5);
    grid on;
    xlim([tstart tend]);
    ylim([-20 20])
    xlabel('Relative time from a trial start [sec]')
    ylabel('Voltage [uV]')
    legend({'Left Hand','Right Hand','Both Feet','Tongue'},'Location','southeast')
    title([dataName     ', ' convertStringsToChars(ChToPlot)]);
    set(gca,'fontsize', 18);
    
    fileName= ['Fig_', dataName ,'_', convertStringsToChars(ChToPlot) , '_ERP'];
    cd(resultPath)
    saveas(h,[fileName,'.fig']);
    saveas(h,[fileName,'.tif']);

    
    % Plotting FTA
    FTrealToPlot = squeeze(trialEEGfftReal1side(:,ChToPlotIndex,:));
    FTimaginaryToPlot = squeeze(trialEEGfftImaginary1side(:,ChToPlotIndex,:));
    
    h = figure;
    fplot = fs*(0:(lFFT/2))/lFFT;
    xlimRange = [fplot(1) fplot(nFreqComp)+1];
    ylimRange = [-900 900];
    subplot(2,1,1)
    hold on
    plot(fplot,mean(FTrealToPlot(:,classID==1),2),'r','LineWidth',2);
    plot(fplot,mean(FTrealToPlot(:,classID==2),2),'b','LineWidth',2);
    plot(fplot,mean(FTrealToPlot(:,classID==3),2),'g','LineWidth',2);
    plot(fplot,mean(FTrealToPlot(:,classID==4),2),'m','LineWidth',2);
    plot(fplot,mean(FTrealToPlot(:,classID==1),2)+std(FTrealToPlot(:,classID==1),[],2),':r','LineWidth',1.5);
    plot(fplot,mean(FTrealToPlot(:,classID==1),2)-std(FTrealToPlot(:,classID==1),[],2),':r','LineWidth',1.5);
    plot(fplot,mean(FTrealToPlot(:,classID==2),2)+std(FTrealToPlot(:,classID==2),[],2),':b','LineWidth',1.5);
    plot(fplot,mean(FTrealToPlot(:,classID==2),2)-std(FTrealToPlot(:,classID==2),[],2),':b','LineWidth',1.5);
    plot(fplot,mean(FTrealToPlot(:,classID==3),2)+std(FTrealToPlot(:,classID==3),[],2),':g','LineWidth',1.5);
    plot(fplot,mean(FTrealToPlot(:,classID==3),2)-std(FTrealToPlot(:,classID==3),[],2),':g','LineWidth',1.5);
    plot(fplot,mean(FTrealToPlot(:,classID==4),2)+std(FTrealToPlot(:,classID==4),[],2),':m','LineWidth',1.5);
    plot(fplot,mean(FTrealToPlot(:,classID==4),2)-std(FTrealToPlot(:,classID==4),[],2),':m','LineWidth',1.5);
    grid on;
    xlim(xlimRange);
    ylim(ylimRange)
    title([dataName     ', ' convertStringsToChars(ChToPlot)]);
    xlabel('Frequency [Hz]')
    ylabel('Real')
    set(gca,'fontsize', 18);
    
    subplot(2,1,2)
    hold on
    plot(fplot,mean(FTimaginaryToPlot(:,classID==1),2),'r','LineWidth',2);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==2),2),'b','LineWidth',2);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==3),2),'g','LineWidth',2);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==4),2),'m','LineWidth',2);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==1),2)+std(FTimaginaryToPlot(:,classID==1),[],2),':r','LineWidth',1.5);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==1),2)-std(FTimaginaryToPlot(:,classID==1),[],2),':r','LineWidth',1.5);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==2),2)+std(FTimaginaryToPlot(:,classID==2),[],2),':b','LineWidth',1.5);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==2),2)-std(FTimaginaryToPlot(:,classID==2),[],2),':b','LineWidth',1.5);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==3),2)+std(FTimaginaryToPlot(:,classID==3),[],2),':g','LineWidth',1.5);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==3),2)-std(FTimaginaryToPlot(:,classID==3),[],2),':g','LineWidth',1.5);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==4),2)+std(FTimaginaryToPlot(:,classID==4),[],2),':m','LineWidth',1.5);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==4),2)-std(FTimaginaryToPlot(:,classID==4),[],2),':m','LineWidth',1.5);
    grid on;
    xlim(xlimRange);
    ylim(ylimRange)
    xlabel('Frequency [Hz]')
    ylabel('Imaginary')
    set(gca,'fontsize', 18);
        
    fileName= ['Fig_', dataName ,'_', convertStringsToChars(ChToPlot) , '_FTA'];
    cd(resultPath)
    saveas(h,[fileName,'.fig']);
    saveas(h,[fileName,'.tif']);
end
