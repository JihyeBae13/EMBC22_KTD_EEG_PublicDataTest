clear; clc; close all

codePath = pwd;
dataPathIni = 'C:\EMBC22_KTDEEG\Codes\Results';
resultPath = 'C:\EMBC22_KTDEEG\Manuscript\Figures';
featureName = 'FTA';
classifierName = 'KTD';
dataName = {'CLASubjectA1601083StLRHand', ...
    'CLASubjectB1510193StLRHand','CLASubjectB1510203StLRHand','CLASubjectB1512153StLRHand',...
    'CLASubjectC1511263StLRHand', 'CLASubjectC1512163StLRHand','CLASubjectC1512233StLRHand', ...
    'CLASubjectD1511253StLRHand', ...
    'CLASubjectE1512253StLRHand', 'CLASubjectE1601193StLRHand','CLASubjectE1601223StLRHand',...
    'CLASubjectF1509173StLRHand', 'CLASubjectF1509283StLRHand'};

nepoch = 100;
ndata = length(dataName);
results = cell(1,ndata);
successRateMean = nan(ndata,nepoch);
successRateStd = nan(ndata,nepoch);
for idata = 1:ndata
    dataPath = [dataPathIni '\' dataName{1,idata} '\' classifierName '_' featureName];
    cd(dataPath)
    results{1,idata} = load('ResultsWparms_stepsize05_10runs.mat');
    successRateMean(idata,:) = mean(results{1,idata}.successRate,1);
    successRateStd(idata,:) = std(results{1,idata}.successRate,1);
end

h=figure;
cd(codePath)
colorSet=varycolor(ndata);
xepoch = 1:nepoch;
hold on
for idata = 1:ndata
    plot(xepoch,successRateMean(idata,:),'LineWidth',2,'Color',colorSet(idata,:));
end
xlim([0 xepoch(50)])
xlabel('Epochs')
ylabel('Average Success Rate')
grid on
set(gca,'fontsize', 18);

cd(resultPath)
saveas(h,['Fig_DataA_' classifierName '_' featureName '.fig']);
saveas(h,['Fig_DataA_' classifierName '_' featureName '.tif']);
