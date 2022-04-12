function [ds, pa, kb, oc] = SetupNewTrial(ds, pa, kb, oc)
% This sets up a new trial.

% Flip a coin to determine which staircase to run for each rotationXheading
% condition (which contributes to one psychometric fucntion
pa.coinFlip = rand(1);
pa.stairCaseCodeRanAlready = [0 0]; % reset these dummy variables indicating whether staircase code ran already

% Generate random dot positions within spherical hull
% by carving out a sphere from a cube, then carving out a smaller sphere
% from that sphere to make a hull
pa.practiceSquare = pa.cubeWidth.*(rand(3,round(pa.nDots))-0.5); % here is where i randomly place dots in a cube on each trial
normCubeVectors = sqrt(sum(pa.practiceSquare .^ 2, 1));
pa.sphereRadius = pa.cubeWidth/2;
pa.sphereExclusionRegionRadius = pa.percentageDotsExcludedFromSphere.*pa.sphereRadius; % exclude dots that are closer than 25% of sphere radius in meters in any direction. just at beginning of trial. exclusion region doesn't move with you bc then dots would be appearing and reappearing and it could be distracting.
pa.sphereExclusion = pa.practiceSquare(:, normCubeVectors < pa.sphereExclusionRegionRadius);
pa.sphereMaster = pa.practiceSquare(:, normCubeVectors < pa.sphereRadius);
normVecs = sqrt(sum(pa.sphereMaster.^ 2, 1));
pa.sphere = pa.sphereMaster(:,normVecs>pa.sphereExclusionRegionRadius);

% Set dummy variables saying participant hasn't responded yet.
kb.responseGiven = 0; % init dummy variable saying that response and feedback haven't been given yet
pa.feedbackGiven = 0;
kb.answer = NaN; % resetting kb.answer on each trial (before the participant responds)

% Increment trial counter
pa.trialNumber = pa.trialNumber + 1; % increment trial counter
pa.thisTrial = pa.trialNumber; % variable to record which trial we're currently on

% Start trial onset timer from current vertical blank time
pa.trialOnset = ds.vbl; % the trial starts....NOW

