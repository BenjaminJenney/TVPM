function [wsa] = OFAnal_Simulated_Real_Head(dirname)

%   OFAnal_Simulated_Head.m
%      usage: OFAnal_Simulated_Real_Head(dirname)
%         by: Charlie Burlingham
%
%       date: 07/29/21
%    purpose: Within-subjects averaging of run data to produce
%             psychometric curves and estimates of bias and threshold
%
%             Trial Types:
%             6 headings (corresponding to 3 rotation speeds)


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


trialVec = 1:data(1).pa.nTrials;


for kk = 1:numRuns % run
        for jj = 1:data(1).pa.numHeadingAngles % headings
            %preallocate
            responsesPerCondition{jj}  =[];
            groundTruthHeadingAnswer{jj} = [];
            post{jj} = [];
            heading{jj} = [];
        end
end

% here would be the place to cut the runs at the beginning! those trials are where the staircase is converging
% and you don't want to include them in the analysis- they make the
% psychometric function weird.

% also maybe should use two 3-up-1-down staircases in opposite
% directions instead. as hormet suggested. makes it faster
firstRunToInclude = 1; % ATTN: start at run N because we don't want trials where staircase is converging

for kk = firstRunToInclude:numRuns % run
    for jj = 1:data(1).pa.numHeadingAngles % headings
        
        indexesPerCondition{jj} = trialVec(data(kk).pa.headingAngleVecRandThisRun == jj);
        
        responsesPerCondition{jj} = [responsesPerCondition{jj} data(kk).pa.leftRightResponse(indexesPerCondition{jj})'];
        
        post{jj} = [post{jj} data(kk).pa.postPosDeg(indexesPerCondition{jj})'];
        
        heading{jj} = [heading{jj} data(kk).pa.heading(indexesPerCondition{jj})'];
        
    end
end


for jj = 1:data(1).pa.numHeadingAngles % headings
    groundTruthHeadingAnswer{jj} = heading{jj} > post{jj};
    xVar{jj} = heading{jj}-post{jj};
    
    dataMatrix{jj} = [xVar{jj}' responsesPerCondition{jj}' ones(length(responsesPerCondition{jj}'),1) ];
     
    dataMatrixUncut{jj} = dataMatrix{jj};
    
    [iind jind] =ind2sub(size(dataMatrix{jj}), find(isnan(dataMatrix{jj})) );
    
    dataMatrix{jj}(unique(iind),:) = [];
end





%% fit with psignifit

headings = [-15 -10 -5 5 10 15];

% fit and plot psychometric functions, depending on how much data there is.
% also get confidence intervals on all the parameter fits
PsignifitParams = []; ninetyfivepcConfintBias = []; ninetyfivepcConfintVar = []; devianceVec = [];marginalMu = {}; xVals = {};

subplotcounter = 0;
for jj = 1:5%data(1).pa.numHeadingAngles % headings
    
    inputForPsignifit = dataMatrix{jj};
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
    
    plotOptions.xLabel         = 'Heading-Post in Visual Degrees';     % xLabel
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
    subplot(3,2,subplotcounter)
    [hline1] = plotPsych(fitOutput,plotOptions);
    title(['Real, ', 'Heading = ', num2str(headings(jj))], 'Fontsize',24)
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
savefig('psychFuncs')

%{
save('PsignifitParams','PsignifitParams')
save('ninetyfivepcConfintBias','ninetyfivepcConfintBias')
save('ninetyfivepcConfintVar','ninetyfivepcConfintVar')
save('devianceVec','devianceVec')
save('marginalMu','marginalMu')
save('xVals','xVals')
save('threshold','threshold')
%}

%%
figure(2)
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


headingsVec = [-15 -10 -5 5 10 15];
bias = 0-PsignifitParams(:,1);
biasCI = 0-ninetyfivepcConfintBias;

biasPooled3RotSpeeds = ( -1.*flipud(bias(1:3)) + (bias(4:6)))/2; % % CSB change this to pool over the binary data (after reversing sign for positive heading), then fit model, and get bias
biasCIPooled3RotSpeeds = ( -1.*flipud(biasCI(1:3,:)) + (biasCI(4:6,:)))/2;

% Optionally plot bias with 95% confidence intervals
subplot(2,1,1)
errorbar((1:3)',biasPooled3RotSpeeds,biasPooled3RotSpeeds-biasCIPooled3RotSpeeds(:,1),biasCIPooled3RotSpeeds(:,1)-biasPooled3RotSpeeds,'-or');
set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'5' , '10' '15'},'Fontsize',12)
ylabel('Bias (deg)','Fontsize',24)
title('Bias')
%ylim([-15 15])
%xlim([.7 3.2])
%hline(0,'k')
set(gca,'FontSize',20)
savefig('biasPlotPsignifit')

threshold = PsignifitParams(:,2)/(2*1.96);
thresholdCI = ninetyfivepcConfintVar./(2*1.96);

threshPooled3RotSpeeds = (flipud(threshold(1:3)) + (threshold(4:6)))/2;
threshCIPooled3RotSpeeds = (flipud(thresholdCI(1:3,:)) + (thresholdCI(4:6,:)))/2;

% Optionally plot threshold with 95% confidence intervals
subplot(2,1,2)
errorbar(1:3,threshPooled3RotSpeeds,threshPooled3RotSpeeds-threshCIPooled3RotSpeeds(:,1),threshCIPooled3RotSpeeds(:,2)-threshPooled3RotSpeeds,'ko')
set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'5' , '10' '15'},'Fontsize',12)
xlabel('Aboslute heading (deg)','Fontsize',24)
ylabel('Threshold (deg)','Fontsize',24)
title('Threshold','Fontsize',30)
%set(gca,'XTickLabel',{'Left Fast', 'Left Slow', 'Pure Translation', 'Right Slow', 'Right Fast'},'Fontsize',12)
%set(gca,'XTick',1:size(sortedHeading,1))
%hline(0)
savefig('varPlotPsignifitwCIs')


keyboard

end