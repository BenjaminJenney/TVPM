function [ds,pa] = SetupParameters(ds,pa)

%% Basic experimental specs 

% Sound parameters
[yMiss, pa.freqMiss] = psychwavread('Sounds/Swish.wav');
pa.waveDataMiss = yMiss';
pa.nrchannelsMiss = size(pa.waveDataMiss,1); % Number of rows == number of channels.

[yHit, pa.freqHit] = psychwavread('Sounds/cowbell.wav');
pa.waveDataHit = yHit';
pa.nrchannelsHit = size(pa.waveDataHit,1); % Number of rows == number of channels.

% Setup sounds
InitializePsychSound(1);
pa.handleHit = PsychPortAudio('Open', [], [], 0, pa.freqHit, pa.nrchannelsHit);
PsychPortAudio('FillBuffer', pa.handleHit, pa.waveDataHit);
pa.handleMiss = PsychPortAudio('Open', [], [], 0, pa.freqMiss, pa.nrchannelsMiss);
PsychPortAudio('FillBuffer', pa.handleMiss, pa.waveDataMiss);

% setup date and data directory
pa.date = datestr(now,30);
s = RandStream('mt19937ar','Seed','shuffle'); % this sets up a new state
RandStream.setGlobalStream(s); % sets the stream according to the state
defaultStream = RandStream.getGlobalStream; % for 64-bit; getCurrentStream for old 32-bit
pa.savedState = defaultStream.State; %% these two lines allow us to recover the subject's experience
pa.baseDir = pwd;
pa.Desktop =  extractBefore(pwd,'\Optic-Flow-Stereo');
pa.Desktop = [pa.Desktop '\'];


%% Parameters of all visual objects - stimulus, targets, etc.

pa.instructionDisparity = 200; % arbitrary and in pixels, for instruction text

pa.cubeWidth = 60; % m. the width of the initial cube = diameter of spherical dot hull...  the sphere radius is half of this 
pa.percentageDotsExcludedFromSphere = .25; % percentage of sphere exluded- determines size of exclusion region in meters (just multiply diameter by this percentage).
% this number times the sphere radius (cubeWidth/2) = the distance of the
% closest dot at the beginning of the trial
pa.dotDensity = 1; % ATTN: need to actually compute this. but for now not using this and just setting the number of dots by hand below...  %((2.*ds.viewingDistance).*ds.pixelsPerM).*tand(0.5); % in dots/deg. per plane, calculated using the formula PPD = 2dr.*tan(0.5 degrees),
%where d is the viewing distance and r is the display resolution in pixels per length
pa.nDots = 1e5;%(pa.dotDensity.*ds.degreesPerM).*pa.planeWidth;
% CSB: attn, need to compute dots/deg on monitor, by finding out the dots
% per m on the monitor then converting using ds.degreesPerM.  

%Colors:
pa.red = [1 0 0];
pa.black = [0 0 0 1];

% response target AKA post parameters
pa.postDiameter = 25;

% pursuit taget parameters
pa.pursuitTargetDist = 1; %pa.cubeWidth/4; % depth of target in m.  csb july  18 2021, put the post very close to observer, like ehrlich, banks 1998 paper to overcome path perception issue / curved path percept. 1 m is what they used too

% fixation dot parameters - must be smaller than response target
pa.fixationDiameter = 20;
pa.fixationDist = -pa.cubeWidth/4 +.5; % depth of fixation dot in m
pa.fixationVertexPos = [0 0 pa.fixationDist];

pa.stimulusDuration_sec = .8; % seconds
pa.nFrames = 90 * pa.stimulusDuration_sec; % NEED TO FIX: 80 fps is by no means exact.
pa.viewingDistance_m = .05; % CHECK WITH DAVID: this is almost certainly incorrect. Is traditional viewing distance meaningful in vr? If so are we interested in distance to OLED's or lens'?
pa.dty = 0; % dtx and dty are accumulators inside the render loop for the phase change of the animated plaids.
pa.dtx = 0;
%% experimental structure/design
pa.trialNumber = 0; % start at trial 0

% trial timing setup
pa.targetMotionDuration = pa.stimulusDuration_sec; % s
pa.trialNumber = 0; % start at 0
pa.itiLength = 2; % seconds
if strcmp(ds.experimentType,'real')
    pa.recenterScreenLength = 1; %0.25 seconds for simulated condition, 0.7 for real head rotation condition %recentering screen duration in seconds
elseif strcmp(ds.experimentType,'simulated') || strcmp(ds.experimentType, 'stabilized') || strcmp(ds.experimentType, 'tvpmfull') || strcmp(ds.experimentType, 'tvpmplanes')
    pa.recenterScreenLength = 1; %0.2 seconds for simulated condition, 0.7 for real head rotation condition %recentering screen duration in seconds
end

% number of trials
pa.nTrials = 60; % Must be divisible by 6 bc there are 3 rot velocities and 6 headings. % csb june 21 2021

% number of runs
if pa.TRAINING == 1
    pa.numRuns = 2; %2 for practice trials
elseif pa.TRAINING == 0
    pa.numRuns = 18;  % number of runs in session. 36 runs per condition, and that takes ~2 hours. so 18 runs per session
                      % ~5 minutes per run, when trials are 2.5 seconds long
end

% setup rotation velocity, translation velocity, and heading directions
pa.transSpeed = .5;%.5; % m/s
%transpeed and planedepth get multiplied by eachother.
if pa.TRAINING ==1
    pa.rotationVelocities = 0; % deg/s.  if training, use no rotation only
else
    pa.rotationVelocities = 0; % deg/s. csb june 21 2021
end

pa.numRotationVelocities = length(pa.rotationVelocities);
numTrialsPerRotVel = pa.nTrials/pa.numRotationVelocities;
pa.rotSpeedVec = ones(numTrialsPerRotVel,pa.numRotationVelocities).* [1:length(pa.rotationVelocities) ];  pa.rotSpeedVec = pa.rotSpeedVec(:); % unroll
pa.trialPermOrder = randperm(pa.nTrials); % will use for both headings and rotation velocities
pa.rotSpeedVecRandThisRun = pa.rotSpeedVec(pa.trialPermOrder);

pa.headingAngles = [-15 -10 -5 5 10 15]; % csb june 21 2021. deg. need to add a small number so staircase thresh never actually converges to exactly this heading and throws error or renders wrong feedback.
pa.numHeadingAngles = length(pa.headingAngles);
pa.headingAnglesIndexes = 1:pa.numHeadingAngles;
numTrialsPerHeadingAngles = pa.nTrials/pa.numHeadingAngles;
pa.headingAngleVec = repmat(pa.headingAnglesIndexes,[1,numTrialsPerHeadingAngles]); % this will equalize the number of different headings at each rot velocity
pa.headingAngleVecRandThisRun = pa.headingAngleVec(pa.trialPermOrder);  % randomly permute headings with same order as rotation velocities

% Make staircases for each velocity speed
pa.numStaircases = 2; % how many staircases to interleave
pa.stepSize = 4; %amount staircase is being shifted by initially
pa.minStairThreshold = -20; % deg   Maximum distance btw heading and post (heading + rotation rate * trial length (sec)) must be less than monocular horizontal field of view (<35 deg conservatively). Also depends on viewing distance inside hmd. see http://doc-ok.org/?p=1414
pa.maxStairThreshold = 20; % deg

for ii = 1:pa.numRotationVelocities
    for jj = 1:pa.numStaircases
        for kk = 1:pa.numHeadingAngles
        
            if jj == 1
                pa.initStairCaseHeading = 15;
            elseif jj == 2
                pa.initStairCaseHeading = -15;
            end
            
        pa.stairs(ii,jj,kk) = upDownStaircase(1,1,pa.initStairCaseHeading,pa.stepSize,'levitt'); % init staircases with just initial values!
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

% pre-allocate response, post position, and correctness vectors
pa.response = NaN(pa.nTrials,3); % Eventual response matrix. 3 is num columns: [pa.trialNumber, pa.responseTime, pa.currentTime]
pa.postPosDeg = NaN(pa.nTrials,1);
pa.correctness = NaN(pa.nTrials,1);
pa.whichStairCase = NaN(pa.nTrials,1);
pa.leftRightResponse = NaN(pa.nTrials,1);
pa.RotationVel = NaN(pa.nTrials,1);
pa.heading = NaN(pa.nTrials,1);

pa.quitFlag = 0; % don't give up

%% Params for TVPM-FULL and TVPM-PLANES+MASK
pa.sf = 2;
pa.apertureDia_px = 20;
end