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
elseif strcmp(ds.experimentType,'simulated') || strcmp(ds.experimentType, 'stabilized')
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




%% textures


InitializeMatlabOpenGL(1);
 %BJ: We are only interested in the API reported res of 1344 * 1600
 % The render resolution is higher than the OLED res to counteract the barrel distortion/pin-cushion distortion caused by the lenses.

Screen('BeginOpenGL', ds.w); % Setup the OpenGL rendering context
glEnable(GL.TEXTURE_2D); % Enable 2D texture mapping

 % meters % This is according to Bas
%% THREE PLANES

%Set up texture for the three planes
planeTextureWidth_px   = (ds.hFOV_perPersonAvg/3) * ds.px_per_deg;
planeTextureHeight_px  = (ds.vFOV_perPersonAvg) * ds.px_per_deg;
 
xHalfTextureSize = planeTextureWidth_px/2;
yHalfTextureSize = floor(planeTextureHeight_px/2);

[x,y] = meshgrid(-xHalfTextureSize+1:xHalfTextureSize,-yHalfTextureSize+1:yHalfTextureSize);

spatialFrequency = 1;% cycles/deg % sf of plaid
%textureData  = makeSomePlaid(ds, pa, x, y, u, v);

%planeTexture = Texture(GL, GL_TEXTURE_2D, xHalfTextureSize*2, yHalfTextureSize*2, textureData);
% End texture for the three planes

% Setup vertices for the three planes

planeDepths_m     = [-.5, -1, -2]; % near, mid, far.
numPlanes         = length(planeDepths_m);
ds.planeWidths_m  = zeros(1, numPlanes);
ds.planeHeights_m = zeros(1, numPlanes);
% 
 for i = 1:numPlanes
    ds.planeWidths_m(i)  = 2 * -planeDepths_m(i) * tand((ds.hFOV_perPersonAvg/numPlanes)/3);
    ds.planeHeights_m(i) = 2 * -planeDepths_m(i) * tand(ds.vFOV_perPersonAvg/2);
 end
% 
% numVertices = 4;
% corners = [0 0;
%     1 0;
%     1 1;
%     0 1];
% 
% ds.planes = zeros(1, numPlanes);
% 
% for i = 1:numPlanes
%     
%     ithPlaneVertices =[-ds.planeWidths_m(i)/2 -ds.planeHeights_m(i)/2 planeDepths_m(i) ;... 
%                         ds.planeWidths_m(i)/2 -ds.planeHeights_m(i)/2 planeDepths_m(i) ;...
%                         ds.planeWidths_m(i)/2  ds.planeHeights_m(i)/2 planeDepths_m(i) ;...
%                        -ds.planeWidths_m(i)/2  ds.planeHeights_m(i)/2 planeDepths_m(i) ]';
%                    
%     %ds.initPositions{i} = opticFlow(ds.planeWidths_m(i), ds.planeHeights_m(i), planeDepths_m(i), 50);
%     ds.planes(i) = glGenLists(1);
% 
%     glNewList(ds.planes(i), GL.COMPILE);
%         planeTexture.bind;
%         glBegin(GL.POLYGON);
%         for j = 1:numVertices
%             glTexCoord2dv(corners(j,:));
%             glVertex3dv(ithPlaneVertices(:,j));
%         end
%         glEnd;
%     glEndList(); 
%     
% end
%planeTexture.unbind;


%  ds.point = glGenLists(1);
%  glPointSize(20);
%  glNewList(ds.point, GL.COMPILE);
%          glBegin(GL.POINT);
%          glVertex3dv([0,0,0]);
%          glEnd;
%  glEndList();
% %% MID PLANE
% midPlaneDepth_m    = -1;
% midPlaneWidth_deg  = ds.hFOV_perPersonAvg/3; % 87
% midPlaneHeight_deg = ds.vFOV_perPersonAvg;   % 84
% midPlaneWidth_px   = midPlaneWidth_deg  * px_per_deg;
% midPlaneHeight_px  = midPlaneHeight_deg * px_per_deg; 
% 
% %ds.midPlane_id = glGenTextures(1);
% xHalfTextureSize = midPlaneWidth_px;%ds.xc;  % needs to be in pixels for meshgrid
% yHalfTextureSize = midPlaneHeight_px;
% [x,y] = meshgrid(-xHalfTextureSize+1:xHalfTextureSize,-yHalfTextureSize+1:yHalfTextureSize);
% 
% %textureData = makeSomeNoise(x,y);
% textureData = makeSomePlaid(x,y,spatialFrequency);
% 
% ds.midPlaneTexture = Texture(GL, GL_TEXTURE_2D, xHalfTextureSize*2, yHalfTextureSize*2, textureData);
% 
% corners = [0 0;
%     1 0;
%     1 1;
%     0 1];
% 
% midPlaneWidth_m  = (2 * -midPlaneDepth_m * tand(midPlaneWidth_deg/2)); %
% midPlaneHeight_m = (2 * -midPlaneDepth_m * tand(midPlaneHeight_deg/2)); %
% 
% midPlaneVertices =[-midPlaneWidth_m/2 -midPlaneHeight_m/2 midPlaneDepth_m ;... 
%                     midPlaneWidth_m/2 -midPlaneHeight_m/2 midPlaneDepth_m ;...
%                     midPlaneWidth_m/2  midPlaneHeight_m/2 midPlaneDepth_m ;...
%                    -midPlaneWidth_m/2  midPlaneHeight_m/2 midPlaneDepth_m ]';
% ds.midPlane = glGenLists(1);
% glNewList(ds.midPlane, GL.COMPILE);
% 
% glBegin(GL.POLYGON);
% for i = 1:4
%     glTexCoord2dv(corners(i,:));
%     glVertex3dv(midPlaneVertices(:,i));
% end
% glEnd;
% glEndList(); 
% 
% %% FAR PLANE
% farPlaneDepth_m    = -2;
% farPlaneWidth_deg  = ds.hFOV_perPersonAvg/3; %
% farPlaneHeight_deg = ds.vFOV_perPersonAvg;
% farPlaneWidth_px   = farPlaneWidth_deg * px_per_deg;
% farPlaneHeight_px  = farPlaneHeight_deg * px_per_deg;
% 
% xHalfTextureSize = farPlaneWidth_px;
% yHalfTextureSize = farPlaneHeight_px;
% [x,y] = meshgrid(-xHalfTextureSize+1:xHalfTextureSize,-yHalfTextureSize+1:yHalfTextureSize);
% 
% textureData = makeSomePlaid(x,y,spatialFrequency);
% ds.farPlaneTexture = Texture(GL, GL_TEXTURE_2D, 2*xHalfTextureSize, 2*yHalfTextureSize, textureData);
%  
% corners = [0 0;
%     1 0;
%     1 1;
%     0 1]; 
% 
% ds.farPlaneWidth_m  = (2 * -farPlaneDepth_m * tand(farPlaneWidth_deg/2));
% farPlaneWidth_m  = (2 * -farPlaneDepth_m * tand(farPlaneWidth_deg/2));
% farPlaneHeight_m = (2 * -farPlaneDepth_m * tand(farPlaneHeight_deg/2));
% 
% farPlaneVertices =[-farPlaneWidth_m/2 -farPlaneHeight_m/2 farPlaneDepth_m ;... 
%                     farPlaneWidth_m/2 -farPlaneHeight_m/2 farPlaneDepth_m ;...
%                     farPlaneWidth_m/2  farPlaneHeight_m/2 farPlaneDepth_m ;...
%                    -farPlaneWidth_m/2  farPlaneHeight_m/2 farPlaneDepth_m ]';
% 
% ds.farPlane = glGenLists(1);
% glNewList(ds.farPlane, GL.COMPILE);
% 
% glBegin(GL.POLYGON);
% for i = 1:4
%     glTexCoord2dv(corners(i,:));
%     glVertex3dv(farPlaneVertices(:,i));
% end
% glEnd;
% glEndList(); 

% Close the OpenGL rendering context
Screen('EndOpenGL', ds.w);
end