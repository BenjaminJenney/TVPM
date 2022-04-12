function [ds, pa, kb, oc] = SetupNewRun(ds, pa, kb, oc)
% This sets up a new run.

% Randomize trial order (rotation velocities and heading directions order)
numTrialsPerRotVel = pa.nTrials/pa.numRotationVelocities;
pa.rotSpeedVec = ones(numTrialsPerRotVel,pa.numRotationVelocities).* [1:length(pa.rotationVelocities) ];  pa.rotSpeedVec = pa.rotSpeedVec(:); % unroll
pa.trialPermOrder = randperm(pa.nTrials); % will use for both headings and rotation velocities
pa.rotSpeedVecRandThisRun = pa.rotSpeedVec(pa.trialPermOrder);

numTrialsPerHeadingAngles = pa.nTrials/pa.numHeadingAngles;
pa.headingAngleVec = repmat(pa.headingAnglesIndexes,[1,numTrialsPerHeadingAngles]); % this will equalize the number of different headings at each rot velocity
pa.headingAngleVecRandThisRun = pa.headingAngleVec(pa.trialPermOrder);  % randomly permute headings with same order as rotation velocities

% Output and carry over current staircase thresholds to next run
%, but re-initalize everything else (reset stepsize and staircase update 
% algorithm rule, so it starts over again as it would if we re-initialized)
clear thresholdCarryOver;
for ii = 1:pa.numRotationVelocities
    for jj = 1:pa.numStaircases
        for kk = 1:pa.numHeadingAngles
            thresholdCarryOver(ii,jj,kk) = pa.stairs(ii,jj,kk).threshold;
        end
    end
end

% We also have to clean up struct fields so that the new structs have
% the right fields.. just bc "upDownStaircase.m" is dumb
pa = rmfield(pa,'stairs');

for ii = 1:pa.numRotationVelocities
    for jj = 1:pa.numStaircases
        for kk = 1:pa.numHeadingAngles
            pa.stairs(ii,jj,kk) = upDownStaircase(1,1, thresholdCarryOver(ii,jj,kk),pa.stepSize,'levitt'); % init staircases with threshold from previous run
            pa.stairs(ii,jj,kk).minThreshold = pa.minStairThreshold; % deg. min 28 deg. Maximum post position must be less than rotation rate * trial length (sec) + monocular horizontal field of view (<35 deg conservatively). Also depends on viewing distance inside hmd. see http://doc-ok.org/?p=1414
            pa.stairs(ii,jj,kk).maxThreshold = pa.maxStairThreshold; % deg. max 28 deg.            
        end
    end
end

for ii = 1:pa.numRotationVelocities
    for jj = 1:pa.numStaircases
        for kk = 1:pa.numHeadingAngles
            pa.stairs(ii,jj,kk).strength = []; % must preallocate fields that will be added later, so struct fields are the same between new and old staircase structs
            pa.stairs(ii,jj,kk).direction = [];
            pa.stairs(ii,jj,kk).reversals = [];
        end
    end
end

% Reset head pos info fields bc they are HUGE and will make data files
% keep growing bigger after every run
oc.modelViewDataLeft = [];
oc.modelViewDataRight = [];
pa.practiceSquare = [];

% reset readyToBeginRun
pa.readyToBeginRun = 0;

% Increment run counter
pa.runNumber = pa.runNumber + 1; 

% Reset the trial counter to 1 for new run
pa.trialNumber = 1;

% Start run onset timer from current vertical blank time
pa.runOnset = ds.vbl; % the run starts....NOW


end