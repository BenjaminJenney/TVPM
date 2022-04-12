function [wsa] = OFAnal_Simulated_Real_Head(dirname)

%   OFAnal_Simulated_Head.m
%      usage: OFAnal_Simulated_Real_Head(dirname)
%         by: Charlie Burlingham
%
%       date: 07/25/19
%    purpose: Within-subjects averaging of run data to produce 9
%             psychometric curves, one per condition, for either simulated
%             or real head rotation conditions
%
%             Trial Types:
%             3 rotation velocities (-5,0,5 deg/s) x 3 headings (-10,0,10
%             deg)


%close all; %rng('shuffle')

contents = dir(dirname);

matches = strfind({contents.name},'.mat');

contentscell = struct2cell(contents);
filenames = contentscell(1,:);

cd(dirname);
numRuns = 0;

for ii = 1:length(filenames)
    if filenames{ii}(1) ~= '.'
        numRuns = numRuns+1;
        
        filenames2{numRuns} = extractBefore(filenames{ii},'.mat');
    end
end

for ii = 1:numRuns
    data(ii) = load(sprintf('%s',filenames2{ii}));
end



% 1 : -5 rotation and -10 heading
% 2 : -5  rotation and 0 deg heading
% 3 : -5  rotation and 10 deg heading
% 4 : 0 rotation and -10 heading
% 5 : 0 rotation and 0 heading
% 6 : 0 rotation and 10 heading
% 7 : 5 rotation and -10 heading
% 8 : 5 rotation and 0 heading
% 9 : 5 rotation and 10 heading

trialVec = 1:data(1).pa.nTrials;

for kk = 1:numRuns % run
    for ii = 1:data(1).pa.numRotationVelocities % rotations
        for jj = 1:data(1).pa.numHeadingAngles % headings
            %preallocate
            responsesPerCondition{ii,jj}  =[];
            groundTruthHeadingAnswer{ii,jj} = [];
            post{ii,jj} = [];
        end
    end
end


% here would be the place to cut the runs at the beginning! those trials are where the staircase is converging
% and you don't want to include them in the analysis- they make the
% psychometric function weird. 

% also maybe should use two 3-up-1-down staircases in opposite
% directions instead. as hormet suggested. makes it faster
firstRunToInclude = 1; % ATTN: start at run N because we don't want trials where staircase is converging

for kk = firstRunToInclude:numRuns % run
    for ii = 1:data(1).pa.numRotationVelocities % rotations
        for jj = 1:data(1).pa.numHeadingAngles % headings
            
            indexesPerCondition{ii,jj} = trialVec(data(kk).pa.headingAngleVecRandThisRun == jj & data(kk).pa.rotSpeedVecRandThisRun' == ii);
            
            responsesPerCondition{ii,jj} = [ responsesPerCondition{ii,jj};  data(kk).pa.leftRightResponse(indexesPerCondition{ii,jj}) ];
            
            postPosPerConditionPerRun{ii,jj} = data(kk).pa.postPosDeg(indexesPerCondition{ii,jj});
            
            post{ii,jj}  = [post{ii,jj} postPosPerConditionPerRun{ii,jj}'];
            
            
            groundTruthHeadingAnswer{ii,jj} = [groundTruthHeadingAnswer{ii,jj}; data(kk).pa.heading(indexesPerCondition{ii,jj})' <  postPosPerConditionPerRun{ii,jj}  ]; % is actual heading left or right of post?
            
            
        end
    end
end

for ii = 1:data(1).pa.numRotationVelocities % rotations
    for jj = 1:data(1).pa.numHeadingAngles % headings
        % load in data in psignifit format: input variable, binary response, number of trials at that that input level
        dataMatrix{ii,jj} = [post{ii,jj}'  responsesPerCondition{ii,jj} ones(length(responsesPerCondition{ii,jj}),1)];
        
        % remove bad trials
        clear badTrialInds;
        badTrialInds =  isnan(dataMatrix{ii,jj}(:,2))  | dataMatrix{ii,jj}(:,2) == -1 ;  % if they didn't respond or pressed wrong button delete that trial's data
        dataMatrix{ii,jj}(badTrialInds,:) = [];
        
        % sort everything
        [~, indicesSorted] = sort(dataMatrix{ii,jj}(:,1),1); % sort by input variable aka post positions
        dataMatrix{ii,jj} = dataMatrix{ii,jj}(indicesSorted,:);
        
        dataMatrix{ii,jj}(:,2) = ~dataMatrix{ii,jj}(:,2); % flip responses! so psychometric function can be fit- it doesn't fit it when it's not flipped
        
        
        % ok now make bins!! 
        
        % first way: unique levels get own bin.
        % take unique target positions and sum over all
        %those new sub-matrices to get new data matrix with number of
        %trials at each unique target position and number of trials correct
        %at each target position:
        %{
        [g g2 g3] = unique(dataMatrix{ii,jj}(:,1)); % get unique trials at each target position
        
        newSubMatrix = [];  
        for kk = 1:length(g2)
            
            indsUnique = [];
            indsUnique = find(g3 == kk);
            
            if ~isempty(indsUnique)
                
                newSubMatrix = [newSubMatrix; [g(kk) sum(dataMatrix{ii,jj}(indsUnique,2:3),1)] ];
            end
        end
        
        dataMatrix{ii,jj} = newSubMatrix;
        %}
        
        
        % second way, choose number of bins a priori, then force tested levels
        % into newly created bins
        bounds = [min(dataMatrix{ii,jj}(:,1)) max(dataMatrix{ii,jj}(:,1))]; % smallest and largest x-values
        nBins = 12; %sum(abs(bounds))/4;  divides nicely to give integer number of trials per bin

        dataMatrix{ii,jj} = (BinData(nBins,bounds,dataMatrix{ii,jj}(:,1:2)'))';
       
        dataMatrix{ii,jj} = [dataMatrix{ii,jj}(:,1) dataMatrix{ii,jj}(:,3) dataMatrix{ii,jj}(:,2)]; % reorganize
        
        %{
        % to pool data across headings...
        if jj == 1
           dataMatrix{ii,jj}(:,1) = dataMatrix{ii,jj}(:,1) + 10; 
        elseif jj == 3
             dataMatrix{ii,jj}(:,1) = dataMatrix{ii,jj}(:,1) - 10; 
        end
        %}
            
    end
    
    newDataMatrix{ii} = [dataMatrix{ii,1}; dataMatrix{ii,2}; dataMatrix{ii,3}];
end

keyboard
%% fit with psignifit

headings = [-10 0 10];
rotations = [-5 0 5];

% fit and plot psychometric functions, depending on how much data there is.
% also get confidence intervals on all the parameter fits
PsignifitParams = []; ninetyfivepcConfintBias = []; ninetyfivepcConfintVar = []; devianceVec = [];marginalMu = {}; xVals = {};

subplotcounter = 0;
for ii = 1:data(1).pa.numRotationVelocities % rotations
    for jj = 1:data(1).pa.numHeadingAngles % headings
        
        inputForPsignifit = dataMatrix{ii,jj};
        options.sigmoidName  = 'norm';
        options.expType = 'equalAsymptote';
        
        %optional - if you want to fix lambda at some value
        %options.fixedPars = NaN(5,1);
        %options.fixedPars(3) = .1;
        
        % optional - If you want lambda to be in the range 0 to .1
        %{
        priorLambda = @(x) (x>=0).*(x<=.1);
        options.priors{3} = priorLambda; % set a prior on lambda so it is less than or equal to .1
        options.borders = nan(5,2); % adjust borders because of lambda prior
        options.borders(3,:)=[0,.1];
        %}
        
        fitOutput = psignifit(inputForPsignifit,options);
        %threshold(ii,jj) = getThreshold(fitOutput,0.75)-getThreshold(fitOutput,0.50); % should be very close to sigma, aka the second row of Psignifitparams / 2*1.96
        
        threshold50(ii,jj) = getThreshold(fitOutput,.5);
        
        plotOptions.xLabel         = 'Post Direction in Visual Degrees';     % xLabel
        plotOptions.yLabel         = '% Leftward Responses';    % yLabel
        plotOptions.CIthresh       = true;
        %{
        if ii < 3 % To plot the left, PT, and right curves in independent colours!
            plotOptions.dataColor      = [0 0 1];  % color of the data
        elseif ii == 3
            plotOptions.dataColor      = [0 1 0];  % color of the data
        elseif ii>3
            plotOptions.dataColor      = [1 0 0];  % color of the data
        end
        %}
        
        plotOptions.dataColor      = [0 0 1]; 
        
        plotOptions.lineColor = plotOptions.dataColor;
        
        
        %figure
        subplotcounter = subplotcounter+1;
        figure(1)
        subplot(3,3,subplotcounter)
        [hline1] = plotPsych(fitOutput,plotOptions);
        title(['Simulated, ', 'Heading = ', num2str(headings(jj)),' , Rotation = ', num2str(rotations(ii)), ' deg/s'], 'Fontsize',24)
        %title(['Real, ', 'Rotation = ', num2str(rotations(ii)), ' deg/s'], 'Fontsize',24)
        xlim([-30 30])

        
        %hold on
        
        
        PsignifitParams = [PsignifitParams; fitOutput.Fit(1) fitOutput.Fit(2)];
        
        marginalMu = [marginalMu fitOutput.marginals{1, 1}];
        
        xVals = [xVals fitOutput.X1D{1, 1}];
        
        ninetyfivepcConfintBias = [ninetyfivepcConfintBias; fitOutput.conf_Intervals(1,:,1)];
        ninetyfivepcConfintVar = [ninetyfivepcConfintVar; fitOutput.conf_Intervals(2,:,1)];
        % The order of reported parameters is
        % [threshold,width,lambda,gamma,eta]
        
        % optional- get the deviance of the model from the real data
        [devianceResid devianceVal] = getDeviance(fitOutput);
        devianceVec = [devianceVec devianceVal];
        
        % Optional if you want to plot deviance of model from real data
        %figure;plotsModelfit(fitOutput);
        
        %clear options; clear plotOptions; clear fitOutput;
   end
end


%{
save('PsignifitParams','PsignifitParams')
save('ninetyfivepcConfintBias','ninetyfivepcConfintBias')
save('ninetyfivepcConfintVar','ninetyfivepcConfintVar')
save('devianceVec','devianceVec')
save('marginalMu','marginalMu')
save('xVals','xVals')
save('threshold','threshold')
%}
keyboard
%%
figure(1)
%savefig('EMSEMcurvesPlot')

%{
figure
stem(PsignifitParams(:,1))
xlabel('Condition')
ylabel('Fitted Mu (visual degrees)')
title('Bias Plot Psignifit','Fontsize',24)
set(gca,'XTick',1:size(sortedHeading,1))
hline(0)
savefig('biasPlotPsignifit')
%}


headingsVec = [-10 0 10 -10 0 10 -10 0 10]';
rotationVectors = [-5 -5 -5 0 0 0 5 5 5]';
bias = PsignifitParams(:,1)-headingsVec;
biasCI = ninetyfivepcConfintBias-headingsVec ;
% Optionally plot bias with 95% confidence intervals
figure
subplot(3,1,1)
errorbar((1:3)',bias([1,4,7]),bias([1,4,7])-biasCI([1,4,7],1),biasCI([1,4,7],2)-bias([1,4,7],1),'-or');
set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'-5 º/s' , '0 º/s' '+5 º/s'},'Fontsize',12)
title('Heading Left -10º')
ylim([-15 15])
xlim([.7 3.2])
hline(0,'k')
set(gca,'FontSize',20)


subplot(3,1,2)
errorbar((1:3)',bias([2,5,8]),bias([2,5,8])-biasCI([2,5,8],1),biasCI([2,5,8],2)-bias([2,5,8],1),'-or');
set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'-5 º/s' , '0 º/s' '+5 º/s'},'Fontsize',12)
ylim([-15 15])
title('Heading Straight Ahead 0º')
xlim([.7 3.2])
hline(0,'k')
set(gca,'FontSize',20)


subplot(3,1,3)
errorbar((1:3)',bias([3,6,9]),bias([3,6,9])-biasCI([3,6,9],1),biasCI([3,6,9],2)-bias([3,6,9],1),'-or');
ylim([-15 15])
xlabel('Rotation Condition','Fontsize',24)
ylabel('Heading Error (deg)','Fontsize',24)
title('Heading Right 10º')
set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'-5 º/s' , '0 º/s' '+5 º/s'},'Fontsize',12)
xlim([.7 3.2])
hline(0,'k')
set(gca,'FontSize',20)


threshold = PsignifitParams(:,2)/(2*1.96);
thresholdCI = ninetyfivepcConfintVar./(2*1.96);

% Optionally plot variance with 95% confidence intervals
figure
errorbar(1:size(PsignifitParams,1),threshold,threshold-thresholdCI(:,1),thresholdCI(:,2)-threshold,'ko')
xlabel('Rotation Condition','Fontsize',24)
ylabel('Threshold (visual degrees)','Fontsize',24)
title('Threshold Plot Psignifit w/ 95% CIs','Fontsize',30)
%set(gca,'XTickLabel',{'Left Fast', 'Left Slow', 'Pure Translation', 'Right Slow', 'Right Fast'},'Fontsize',12)
%set(gca,'XTick',1:size(sortedHeading,1))
hline(0)
%savefig('varPlotPsignifitwCIs')






%% to plot bias just for 3 rotations (collapsed over headings)

%bias = PsignifitParams(:,1);
%biasCI = ninetyfivepcConfintBias ;

%figure
errorbar((1:3)',bias, bias-biasCI(1:3,1),biasCI(1:3,2)-bias(1:3,1),'-bo','markerSize',3)
set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'-5 º/s' , '0 º/s' '+5 º/s'},'Fontsize',12)
title('Heading Error','Fontsize',30)
title('Simulated vs. Real')
ylim([-10 10])
xlim([.5 3.5])
hline(0,'k-')
xlabel('Rotation Velocity (º/s)')
ylabel('Heading Bias (º)')
box off

end