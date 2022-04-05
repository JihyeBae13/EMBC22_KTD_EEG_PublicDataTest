clear; clc; close all

dataPath = 'C:\Users\jbae2\Desktop\Jihye Bae\20211202 467C 5820\UK\Lab\2020_Bhoj02\Papers\EMBC01_KTDEEG\Codes\NatureData-CLA';
resultPath = 'C:\Users\jbae2\Desktop\Jihye Bae\20211202 467C 5820\UK\Lab\2020_Bhoj02\Papers\EMBC01_KTDEEG\Codes\Results';
dataName = 'CLASubjectC1511263StLRHand';

cd(resultPath)
if exist(dataName,'dir')~= 7
    mkdir(dataName)
end
resultPath = [resultPath '\' dataName];

cd(dataPath)
load([dataName '.mat']);
marker = o.marker; % Marker array
RawEEG = o.data; % EEG data (22 Channels)
ChannelNames = o.chnames; % Channel Names
fs = o.sampFreq; % Sampling frequency
[tlength, nch] = size(RawEEG); % Number of channels
nclass = 3; %Number of classes

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% CHANNEL_DECODIFICATION_3_ CLASSES %%%%%%%%%%%%%%%%%%%%%%%%%%%
j = 1;
k = 1;
l = 1;
m = 1;
for it = 1:tlength-1
    if (marker(it) == 0) && (marker(it+1) == 1) % Class1 = Left Hand
        indexC1(j,1) = it+1;
        j=j+1;
    elseif (marker(it) == 0) && (marker(it+1) == 2) % Class2 = Right Hand
        indexC2(k,1) = it+1;
        k=k+1;
    elseif (marker(it) == 0) && (marker(it+1) == 3) % Class3 = Neutral
        indexC3(l,1) = it+1;
        l=l+1;
    else
        indexND(m,1) = it; % Nothing displayed to validate total trial numbers
        m=m+1;
    end
end

% Validating the trial extraction based on the marker signal
h = figure;
subplot(2,1,1)
hold on
plot(marker,'k')
plot(indexC1,1,'r*')
plot(indexC2,2,'g*')
plot(indexC3,3,'b*')
xlabel('Time Index')
ylabel('Marker')
title('Single Entire Session')
set(gca,'fontsize', 18);
subplot(2,1,2)
hold on
plot(marker,'k')
plot(indexC1,1,'r*')
plot(indexC2,2,'g*')
plot(indexC3,3,'b*')
xlim([1*10^5 1.15*10^5])
xlabel('Time Index')
ylabel('Marker')
title('Selected Time Interval (Zoomed In)')
set(gca,'fontsize', 18);

% cd(resultPath)
% saveas(h,'Fig_TrialExtraction.fig');
% saveas(h,'Fig_TrialExtraction.tif');

%%
%%%%%%%%%%%%%%%%%%%%% DATA SEGMENTATION %%%%%%%%%%%%%%%%%%%%%%%%%

tstart = 0; % Define start time to extract EEG. Relative to the strat of a trial.
tend = 0.85; % Define end time to extract EEG. Relative to the strat of trial.
ntrial = length(indexC1)+length(indexC2)+length(indexC3); % Number of trials
classID = [ones(size(indexC1)).*1; ones(size(indexC2)).*2; ones(size(indexC3)).*3;];
classTrialIndex = [indexC1; indexC2; indexC3];
trialEEG = zeros(tend*fs,nch,ntrial);

for itr = 1: ntrial
    trialEEG(:,:,itr) = RawEEG(classTrialIndex(itr)+tstart*fs:classTrialIndex(itr)+tend*fs-1,:);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Feature1: RAW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ntpoints = size(trialEEG,1);
RAWfeature = nan(ntrial,ntpoints*(nch-1));
for itr = 1:ntrial
    for ich = 1: nch-1
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

nFreqComp = 5;
nFTAfeature = 2*nFreqComp-1;
FTAfeature = nan(ntrial,(nch-1)*nFTAfeature);
for itr = 1:ntrial     
    for ich = 1: nch-1
        FTAfeature(itr,nFTAfeature*(ich-1)+1:nFTAfeature*ich) = [trialEEGfftReal1side(1:nFreqComp,ich,itr); ...
            trialEEGfftImaginary1side(2:nFreqComp,ich,itr)]'; 
    end
end

% cd(resultPath)
% save('Data_EEGfeatureRAW.mat','RAWfeature','classID')
% save('Data_EEGfeatureFTA.mat','FTAfeature','classID')


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
    tplot = tplot(1:end-1)';
    hold on;
    plot(tplot,mean(trialEEGToPlot(:,classID==1),2),'r','LineWidth',2);
    plot(tplot,mean(trialEEGToPlot(:,classID==2),2),'b','LineWidth',2);
    plot(tplot,mean(trialEEGToPlot(:,classID==3),2),'g','LineWidth',2);
    plot(tplot,mean(trialEEGToPlot(:,classID==1),2)+std(trialEEGToPlot(:,classID==1),[],2),':r','LineWidth',1.5);
    plot(tplot,mean(trialEEGToPlot(:,classID==1),2)-std(trialEEGToPlot(:,classID==1),[],2),':r','LineWidth',1.5);
    plot(tplot,mean(trialEEGToPlot(:,classID==2),2)+std(trialEEGToPlot(:,classID==2),[],2),':b','LineWidth',1.5);
    plot(tplot,mean(trialEEGToPlot(:,classID==2),2)-std(trialEEGToPlot(:,classID==2),[],2),':b','LineWidth',1.5);
    plot(tplot,mean(trialEEGToPlot(:,classID==3),2)+std(trialEEGToPlot(:,classID==3),[],2),':g','LineWidth',1.5);
    plot(tplot,mean(trialEEGToPlot(:,classID==3),2)-std(trialEEGToPlot(:,classID==3),[],2),':g','LineWidth',1.5);
    grid on;
    xlim([tstart tend]);
    ylim([-7 7])
    xlabel('Relative time from a trial start [sec]')
    ylabel('Voltage [uV]')
    legend({'Left Hand','Right Hand','Neutral'},'Location','southeast')
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
    ylimRange = [-300 300];
    subplot(2,1,1)
    hold on
    plot(fplot,mean(FTrealToPlot(:,classID==1),2),'r','LineWidth',2);
    plot(fplot,mean(FTrealToPlot(:,classID==2),2),'b','LineWidth',2);
    plot(fplot,mean(FTrealToPlot(:,classID==3),2),'g','LineWidth',2);
    plot(fplot,mean(FTrealToPlot(:,classID==1),2)+std(FTrealToPlot(:,classID==1),[],2),':r','LineWidth',1.5);
    plot(fplot,mean(FTrealToPlot(:,classID==1),2)-std(FTrealToPlot(:,classID==1),[],2),':r','LineWidth',1.5);
    plot(fplot,mean(FTrealToPlot(:,classID==2),2)+std(FTrealToPlot(:,classID==2),[],2),':b','LineWidth',1.5);
    plot(fplot,mean(FTrealToPlot(:,classID==2),2)-std(FTrealToPlot(:,classID==2),[],2),':b','LineWidth',1.5);
    plot(fplot,mean(FTrealToPlot(:,classID==3),2)+std(FTrealToPlot(:,classID==3),[],2),':g','LineWidth',1.5);
    plot(fplot,mean(FTrealToPlot(:,classID==3),2)-std(FTrealToPlot(:,classID==3),[],2),':g','LineWidth',1.5);
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
    plot(fplot,mean(FTimaginaryToPlot(:,classID==1),2)+std(FTimaginaryToPlot(:,classID==1),[],2),':r','LineWidth',1.5);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==1),2)-std(FTimaginaryToPlot(:,classID==1),[],2),':r','LineWidth',1.5);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==2),2)+std(FTimaginaryToPlot(:,classID==2),[],2),':b','LineWidth',1.5);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==2),2)-std(FTimaginaryToPlot(:,classID==2),[],2),':b','LineWidth',1.5);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==3),2)+std(FTimaginaryToPlot(:,classID==3),[],2),':g','LineWidth',1.5);
    plot(fplot,mean(FTimaginaryToPlot(:,classID==3),2)-std(FTimaginaryToPlot(:,classID==3),[],2),':g','LineWidth',1.5);
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