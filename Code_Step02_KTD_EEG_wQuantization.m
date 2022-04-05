%% This code implements offline center out reaching task using KTD algorithm based on EEG data
% The code is framed to conduct a reaching task.
% The reward is assigned when the cursor is closed to the taget: i.e.
% distance to the target is less then a determined threshold.
% This code is updated by Jihye Bae on 12/21/2021

clear
close all
clc

codePath = pwd;
dataPath = 'C:\EMBC22_KTDEEG\Codes\Results';
dataName = 'A09';
featureName = 'RAW';
classifierName = 'KTD';
folderName = [classifierName '_' featureName ];
testRunName = 'stepsize05_10runs';

cd(dataPath)
if exist(dataName,'dir')~= 7
    mkdir(dataName)
end
dataPath = [dataPath '\' dataName]; 

cd(dataPath)
if exist(folderName,'dir')~= 7
    mkdir(folderName)
end
resultPath = [dataPath '\' folderName]; 

cd(dataPath)
dataFileName = ['Data_EEGfeature' featureName '.mat'];
load(dataFileName);
NSinput = eval([featureName 'feature']);
TargetIndex = classID;

%% setting up parameters
[parms.ntrial, parms.nch] = size(NSinput);
parms.ntarget = max(classID);
parms.npossibleAction = parms.ntarget;

% experiment parameters
parms.nepoch = 100; % number of epoches
parms.nMCrun = 10; % numger of Monte Carlos runs

% normalization parameters
parms.NormalizationUpperBound= 1; % Data Normalization factors
parms.NormalizationLowerBound = -1;

% KTD parameters
parms.KTDstepsize = 0.5; % learning rate
parms.KTDkernelBWini = 1; % initial kernel size, one dummy value is tentatively assigned for the first time
parms.KTDkernelBWfactor = 0.3; % multipication facor for the kernel size
parms.KTDquantizationThr = 0.1; % quantization threshold

% RL parameters
parms.RLepsilon = 0.01; % exploration rate for the epeilson greedy method
parms.RLepsilonDecayEpoch= 50; % every 50th epoch, apply epsilonEpochDecayRate to parms.RLepsilon
parms.RLepsilonDecayRate = 0.5;
parms.RLgamma = 0.9; % discounting factor in reward
parms.RLpreward= parms.ntarget-1; % positive reward
parms.RLnreward = -1; % negative reward

% reching task
parms.ReachingRadius = 1;
parms.ReachingDisanceThr = 0.1;
parms.ReachingCenterXY = [0,0];
parms.ReachingTargetXY = func_GenerateTargetXY(parms.ntarget, parms.ReachingRadius, parms.ReachingCenterXY);

% plotting generated target locations
h = figure;
plotThr = 0.2;
plotRange = [parms.ReachingRadius*(-1)-plotThr parms.ReachingRadius+plotThr];
plotTheta = 0:pi/50:2*pi;
plot(parms.ReachingRadius * cos(plotTheta), parms.ReachingRadius * sin(plotTheta),'-.','LineWidth',2)
hold on
plot(parms.ReachingTargetXY(:,1),parms.ReachingTargetXY(:,2),'ro','LineWidth',4,'MarkerFaceColor','r')
plot(0,0,'gs','LineWidth',6,'MarkerFaceColor','g')
grid on
xlim(plotRange)
ylim(plotRange)
axis equal
xlabel('X')
ylabel('Y')
set(gca,'fontsize', 18);

cd(resultPath)
saveas(h,'Fig_TargetXY.tif')
saveas(h,'Fig_TargetXY.fig')

% KTD implementation
ktdReachingIndexAll = cell(1,parms.nMCrun);
sucessAll = cell(1,parms.nMCrun);
QAll = cell(1,parms.nMCrun);
kernelBWtraceAll = cell(1,parms.nMCrun);
successRate = nan(parms.nMCrun,parms.nepoch);
trialIndex = nan(parms.nMCrun,parms.ntrial);
for imcrun = 1:parms.nMCrun
    ikernel = 1;   
    ktdReachingIndex = nan(parms.nepoch,parms.ntrial);
    Q = nan(parms.ntrial*parms.nepoch,parms.npossibleAction);
    kernelBWtrace = nan(parms.ntrial*parms.nepoch,1);
    success = zeros(parms.nepoch,parms.ntarget); % need to assign 0 values in sucess 
    % so that it starts adding +1 whenever a successful trial is observed.

    trialIndex(imcrun,:) = randperm(parms.ntrial);
    
    for iepoch = 1: parms.nepoch
        
        fprintf('Currnet Implementation: "%d" Monte Carlo run and "%d" epoch\n', imcrun, iepoch)
        
        if mod(iepoch,parms.RLepsilonDecayEpoch) == 0
            parms.RLepsilon = parms.RLepsilon*parms.RLepsilonDecayRate;
        end
        
        for itrial = 1: parms.ntrial
            cursorPreXY = parms.ReachingCenterXY;
            % normalization NSinput
            ns_max = max(NSinput(trialIndex(imcrun,itrial),:));
            ns_min = min(NSinput(trialIndex(imcrun,itrial),:));
            if ns_max == ns_min
                fprintf('=====================ERROR========================\n')
                fprintf('NSinput has max and min the same value, at %d th epoch and %d th trial', iepoch, itrial)
                pause
            end
            ns_amp = (parms.NormalizationUpperBound - parms.NormalizationLowerBound)/(ns_max-ns_min);
            ns_off = parms.NormalizationUpperBound - ns_amp*ns_max;
            NormalizedNS = ns_amp*NSinput(trialIndex(imcrun,itrial),:) + ns_off;
            
            % assigning the first unit values in KTD
            if iepoch == 1 && itrial == 1
                Unit(ikernel,:) = NormalizedNS;
                Weight(ikernel,1:parms.npossibleAction) = zeros(ikernel,parms.npossibleAction);
                % assgining dummy kernel size
                kernelBW = parms.KTDkernelBWini;
            else
                                
                % computing the current Q(t)
                indif = bsxfun(@minus,Unit,NormalizedNS);
                inputDifference = sum(indif.*indif,2);
                kernelBW = sqrt(mean(inputDifference))*parms.KTDkernelBWfactor; % computing kernel size
                kernelf = exp(inputDifference/(-2*kernelBW.^2));
                out1 = parms.KTDstepsize*kernelf'*Weight;
                
                maxout1 = find(out1 == max(out1));
                action1 = func_ActionSelection(maxout1,parms.npossibleAction,parms.RLepsilon);
                Q1 = out1(action1);
                Q((iepoch-1)*parms.ntrial-1+itrial,:) = out1;
                % tracking selected action
                ktdReachingIndex(iepoch,itrial) = action1;
                
                % quantization application
                quantizationIndex = find(inputDifference < parms.KTDquantizationThr);
                if isempty(quantizationIndex)
                    ikernel = ikernel+1;
                    Unit(ikernel,:) = NormalizedNS;
                    Weight(ikernel,:) = Weight(ikernel-1,:);
                    centerIndex = ikernel;
                else
                    if length(quantizationIndex)>1
                        fprintf('=====================ERROR========================\n')
                        fprintf('Quantizatio Weight! Please modify your quantization threshold.\n')
                        pause
                    else
                        centerIndex = quantizationIndex;
                    end
                end
                
                % computing the future Q(t+1)
                if itrial < parms.ntrial
                    indif = bsxfun(@minus, Unit, NormalizedNS);
                    inputDifference = sum(indif.*indif,2);
                    out2 = parms.KTDstepsize*exp(inputDifference/(-2*kernelBW.^2))'*Weight;
                    maxout2 = find(out2 == max(out2));
                    action2 = func_ActionSelection(maxout2,parms.npossibleAction,parms.RLepsilon);
                    Q2 = out2(action2);
                end

                % updating cursor xy location
                cursorNextXY = func_CursorUpdateXY(cursorPreXY, action1, parms.npossibleAction, parms.ReachingCenterXY, parms.ReachingRadius);
                cursorDistance = sqrt(sum((cursorNextXY-parms.ReachingTargetXY(classID(trialIndex(imcrun,itrial)),:)).^2));
                
                % evaluating cursor to target distance
                if cursorDistance < (parms.ReachingRadius-parms.ReachingDisanceThr)
                    reward = parms.RLpreward;
                else
                    reward = parms.RLnreward;
                end
                
                % computing TD value
                if itrial == parms.ntrial
                    TD = reward - Q1;
                else
                    TD = reward + parms.RLgamma*Q2 - Q1;
                end
                
                % counting success trials
                if reward == parms.RLpreward
                    success(iepoch,classID(trialIndex(imcrun,itrial))) = ...
                        success(iepoch,classID(trialIndex(imcrun,itrial))) + 1;
                end
                
                % assigning weights
                if isempty(quantizationIndex)
                        Weight(centerIndex,:) = zeros(1,parms.ntarget);
                        Weight(centerIndex,action1) = TD;
                else
                    Weight(centerIndex,action1) = Weight(centerIndex,action1)+TD;
                end
            end
            
            % success rate calculation
            if itrial == parms.ntrial
                successRate(imcrun,iepoch) = (sum(success(iepoch,:)))/parms.ntrial;
                fprintf('Current success rates: %d.\n',successRate(imcrun,iepoch))
            end
            
            kernelBWtrace((iepoch-1)*parms.ntrial+itrial,1) = kernelBW;
        end % end of trial
        
    end % end of epoch
    
    ktdReachingIndexAll{1,imcrun} = ktdReachingIndex;
    sucessAll{1,imcrun} = success;
    QAll{1,imcrun} = Q;
    kernelBWtraceAll{1,imcrun} = kernelBWtrace;
    
    clear Unit
    clear Weight
end % end of MC run

% plotting success rate
h = figure;
plot(successRate','LineWidth',2)
grid on
xlabel('Epoch')
ylabel('Success Rate')
set(gca,'fontsize', 18);

cd(resultPath)
saveas(h,['Fig_SucessRatesAll_' testRunName '.tif'])
saveas(h,['Fig_SucessRatesAll_' testRunName '.fig'])

% saving parameters and performances
cd(resultPath)
save(['ResultsWparms_' testRunName '.mat'], ...
    'parms', ...
    'ktdReachingIndexAll','sucessAll','trialIndex', ...
    'QAll','kernelBWtraceAll', ...
    'successRate')

function ReachingTargetXY = func_GenerateTargetXY(ntarget, ReachingRadius, ReachingCenterXY)
% Generate a x-y target location in a center-out reaching task
% based on the the number of targets.
% The center is located at the origin, (0,0)
%
% ReachingTargetXY = func_GenerateTarget(ntarget)
%
% N(ntarget) regularly spaced center out format
% ntarget: number of targets
% ReachingRadius: ReachingRadius from the center (origin) to the target

theta = 0:2*pi/ntarget:2*pi-2*pi/ntarget;
ReachingTargetXY = round([cos(theta')+ReachingCenterXY(1,1) sin(theta')+ReachingCenterXY(1,2)].*ReachingRadius,3);

end

function action = func_ActionSelection(maxout,npossibleAction,RLepsilon)
% Select one action based on RLepsilon-greedy method
%
% action = func_ActionSelection(maxout,ntarget,RLepsilon)
% 
% maxout: maximum Q value
% ntarget: number of targets
% RLepsilon: RLepsilon in the RLepsilon-greedy method

if rand > 1-RLepsilon
    rnum = rand(1,npossibleAction);
    action = find(rnum == max(rnum));
else
    if isempty(maxout) == 1
        rnum = rand(1,npossibleAction);
        action = find(rnum == max(rnum));
    elseif length(maxout) > 1
        rnum = rand(1,npossibleAction);
        action = find(rnum == max(rnum));
    else
        action = maxout;
    end
end

end

function CursorNextXY = func_CursorUpdateXY(CursorPreXY, selectedAction, npossibleAction, ReachingCenterXY, ReachingRadius)
% Update a x-y cursor location in a center-out reaching task 
% based on a selected action from the current cursor position.
%
% CursorNextXY = func_CursorUpdateXY(CursorPreXY,npossibleAction, ReachingRadius)
%
% CursorPreXY = current X,Y position
% ntarget = number of targets
% ReachingRadius = ReachingRadius from the center (origin) to the target

theta = 0:2*pi/npossibleAction:2*pi-2*pi/npossibleAction;
possibleXY = round([cos(theta')+ReachingCenterXY(1,1) sin(theta')+ReachingCenterXY(1,2)].*ReachingRadius,3);
CursorNextXY = CursorPreXY + possibleXY(selectedAction,:);

end
