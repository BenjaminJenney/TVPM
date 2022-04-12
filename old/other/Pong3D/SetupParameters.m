function [ds,pa] = SetupParametersfor(ds)

% Establishes the parameters for the Oculus 3D Pong experiment using the
% 2nd generation Oculus (DK2)

% At the start of the session, should only update lines 11 and 12


%% Basic experimental specs 

pa.subjectName = 'test';
pa.feedbackFlag = 2;  % 0 --> no feedback;  1 --> only sound-based feedback;  2 --> visual + sound-based feedback

pa.criterion = 8*(pi/180);  % this is the criterion for the distance between the mid-point of the paddle and the midpoint of the ball; distances smaller than this are hits - any overlap counts

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



pa.date = datestr(now,30);
s = RandStream('mt19937ar','Seed','shuffle'); % this sets up a new state
RandStream.setGlobalStream(s); % sets the stream according to the state
defaultStream = RandStream.getGlobalStream; % for 64-bit; getCurrentStream for old 32-bit
pa.savedState = defaultStream.State; %% these two lines allow us to recover the subject's experience
pa.baseDir = pwd;


%% parameters for the texture surround and virtual 'room'

pa.fixationDiskRadius = 2.5; % deg
pa.fixationHalfSquareSize = 1;  % deg - just the 'box' in which the fixation lines will be drawn
pa.fixationHalfSquare = 0.0157; % m 
pa.fixationLineLength = 0.0079; % m 
    
pa.apertureRadius = 7.5;  % deg
   
pa.floorHeight = -1.5; % m
pa.ceilingHeight = 1.5; % m 
    
%% parameters for the adjustable paddle

pa.paddleHalfWidth = 0.005;% m
pa.paddleHalfHeight = 0.005;% m
pa.paddleHalfDepth = 0.01;% m
pa.paddleHeightFactor = 0.0057;% m 
pa.paddleAngle = 0; % deg - start at the rightward position
pa.shiftPaddle = 0.25;
pa.rotationSpeed = pa.shiftPaddle; % this will update according to the observer's response
pa.acceleratePaddle = 0.008;
pa.paddleOrbitShift = 0.1185;% m
 

%% parameters for the target - will establish the speed distributions below

pa.targetMotionDuration = 1; % s
pa.targetContrast = [1 0.15 0.075]; % fully-visible target, 15% and 7.5% contrast
pa.targetRadius = 0.25; % deg 
pa.targetSize = 0.004;% m
  

%% experimental structure/design

pa.trialNumber = 0; % gotta start somewhere

% Set up a full factorial design 
pa.nRepeats = 75; % each target contrast condition gets pa.nRepeats trials - 75*3 = 225 per block 
pa.fullFactorial = fullfact(length(pa.targetContrast)); 
pa.fullFactorial = repmat(pa.fullFactorial,pa.nRepeats,1); % repeat the full factorial design nRepeats times
pa.nTrials = size(pa.fullFactorial,1);
 
% *Gaussian* vx and vz sampling each from a distribution with mean=0 and
% stdev=0.02m - the 0.061 m/s and -0.061 m/s would then be 3 stdevs from the mean,
% allowing us to cut off only a small percentage of speeds that come up
for rtr=1:pa.nTrials 
    % sample the vx and vz independently
        inRange = 0;
        while inRange==0
            pa.xSpeedValues(rtr) = 0.02.*randn(1,1); % mean=0, stdev=0.02
            if pa.xSpeedValues(rtr)<0.061 && pa.xSpeedValues(rtr)>-0.061 % as long as the speed falls between -0.061 m/s and 0.061 m/s, it's acceptable
                inRange = 1;
            end
        end
        inRange = 0;
        while inRange==0;
            pa.zSpeedValues(rtr) = 0.02.*randn(1,1); % mean=0, stdev=0.02
            if pa.zSpeedValues(rtr)<0.061 && pa.zSpeedValues(rtr)>-0.061
                inRange = 1;
            end
        end
end

pa.fullFactorial(:,1) = pa.targetContrast(pa.fullFactorial(:,1)); % instead of having this be indices of the contrast values it instead will have the actual contrast values - much nicer for data analysis
pa.fullFactorial(:,2) = pa.xSpeedValues';
pa.fullFactorial(:,3) = pa.zSpeedValues';
pa.fullFactorial = pa.fullFactorial(randperm(pa.nTrials),:); % Randomly permute order of trial presentation

pa.response = []; % Eventual response matrix
pa.quitFlag = 0; % don't give up
end




