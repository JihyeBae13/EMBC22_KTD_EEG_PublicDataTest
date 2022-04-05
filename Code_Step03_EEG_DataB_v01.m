clear; clc; close all

codePath = pwd;
dataPathIni = 'C:\Users\jbae2\Desktop\Jihye Bae\20211202 467C 5820\UK\Lab\2020_Bhoj02\Papers\EMBC01_KTDEEG\Codes\Results';
resultPath = 'C:\Users\jbae2\Desktop\Jihye Bae\20211202 467C 5820\UK\Lab\2020_Bhoj02\Papers\EMBC01_KTDEEG\Manuscript\Figures';
featureName = 'FTA';
classifierName = 'KTD';
dataName = {'A01', 'A02','A03','A04','A05',...
    'A06', 'A07','A08','A09'};

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
saveas(h,['Fig_DataB_' classifierName '_' featureName '.fig']);
saveas(h,['Fig_DataB_' classifierName '_' featureName '.tif']);
