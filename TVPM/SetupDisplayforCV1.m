function ds = SetupDisplayforCV1(ds)

ds.viewingDistance = 0.45; %0.60; % in m

% Setup Psychtoolbox for OpenGL 3D rendering support and initialize the
% mogl OpenGL for Matlab/Octave wrapper:
InitializeMatlabOpenGL(1);
% PsychDefaultSetup(2); % the input of 2 means: execute the AssertOpenGL command, execute KbName('UnifyKeyNames') routine, AND unifies the color mode and switches from 0-255 to 0.0-1.0 - color part only impacts the current function or script, not ones that are called

ds.screenid = max(Screen('Screens')); % Find the screen to use for display:

ds.multiSample = 8; 
ds.doSeparateEyeRender = 1; % to render two eyes' views

PsychImaging('PrepareConfiguration');

% even in the 'fixed' viewing condition, we still want to track the head
% movement to save it out later - we just don't update the display in
% response to the movement

% Might want to try turning off DebugDisplay to improve latency
ds.hmd = PsychVRHMD('AutoSetupHMD', 'Tracked3DVR', 'LowPersistence TimeWarp FastResponse DebugDisplay', 0); % 'TimingSupport' % on my setup, I have the Oculus runtime installed. As a result, it detects no HMDs if the Oculus is unplugged from the system, and then opens a simulated HMD window (in the wrong orientation), rather than simply displaying to the external monitor
%     PsychVRHMD('SetHSWDisplayDismiss', ds.hmd, -1);
%     PsychVRHMD('Verbosity',4); %3 = info; 4 = debug
%     PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma'); % gamma correction

if ~isempty(ds.hmd) % VR headset found
    ds.oculusConnected = 1;
    ds.hmdinfo = PsychVRHMD('GetInfo', ds.hmd); % information saved to HMDInfo.m in the Working Demo Files folder
else
    % TODO combine loaded files so everything needed is in one place
    fprintf('No VR-HMD available, using default values.\n');
    ds.oculusConnected = 0;
    temp = load('dsHmdInfoCV1.mat');
    ds.hmdinfo = temp.ds.hmdinfo;
    clear temp
    Screen('Preference', 'SkipSyncTests', 1);
end

% Use a fixed head position, regardless of head tracking
% TODO figure out how to read IPD out of headset
% Might be provided as a 3rd argument in GetStaticRenderParameters
ds.defaultIPD = 0.064; % in m
defaultState.modelViewDataLeft = [1 0 0 -ds.defaultIPD/2; 0 1 0 0; 0 0 1 -ds.viewingDistance; 0 0 0 1]; %+/- IPD/2
defaultState.modelViewDataRight = [1 0 0 ds.defaultIPD/2; 0 1 0 0; 0 0 1 -ds.viewingDistance; 0 0 0 1];
ds.defaultState.modelView = {defaultState.modelViewDataLeft, defaultState.modelViewDataRight};

% Oculus CV1 resolution is 1080 x 1200 in each eye
% rect = [0 0 2*1080 1200]./2; % show oculus view windowed on desktop
rect = [0 0 ds.hmdinfo.panelXRes ds.hmdinfo.panelYRes]./2; % show oculus view windowed on desktop.
[ds.w, ds.windowRect] = PsychImaging('OpenWindow', ds.screenid, 0, rect);  % keeps automatically setting the stereo mode to 6 for the oculus - this is because, as indicated in MorphDemo.m: "% Fake some stereomode if HMD is used, to trigger stereo rendering"    
[ds.xc, ds.yc] = RectCenter(ds.windowRect);

% PsychColorCorrection('SetEncodingGamma', ds.w, 1);%1./ds.gammaVals); % set required gamma

% Calculate display properties
ds.hFOV = ds.hmdinfo.fovL(1) + ds.hmdinfo.fovL(2); % in deg - symmetric for fovL and fovR
ds.vFOV = ds.hmdinfo.fovL(3) + ds.hmdinfo.fovL(4); % in deg - vertical field of view
ds.dFOV = sqrt(ds.hFOV^2 + ds.vFOV^2);

ds.Height = 1.7614; %value used in DK1 version %ds.viewingDistance*2*tand(ds.vFOV/2); % virtual height of the surround texture in meters, based on viewing distance - we want this to relate to the shorter dimension of the display
ds.halfHeight = ds.Height/2;
ds.Width = 1.7614; %value used in DK1 version % ds.viewingDistance*2*tand(ds.hFOV/2); % virtual width of the surround texture in meters, based on viewing distance - we want this to relate to the longer dimension of the display
ds.halfWidth = ds.Width/2;
ds.xc = RectWidth(ds.windowRect)/2; % the horizontal center of the display in pixels
ds.yc = RectHeight(ds.windowRect)/2; % the vertical center of the display in pixels
ds.textCoords = [ds.yc ds.xc];

ds.mPerDegree =  .1; %ds.viewportWidthM/ds.hFOV; % Has to be at the fixation distance
ds.degreesPerM = 1./ds.mPerDegree;

% Careful, because this suggests ppd is uniform across the display
ds.pixelsPerDegree = sqrt(RectHeight(ds.windowRect)^2 + RectWidth(ds.windowRect/2)^2) / ds.dFOV; 
ds.frameRate = 1./ds.hmdinfo.videoRefreshDuration;
ds.ifi = Screen('GetFlipInterval', ds.w); % Get duration of a single frame

% Initialize for saving purposes
ds.eyePoseModelViewDataLeft = []; % may as well save the model view matrix data as well - the hope is that this covers all of the critical information to later go back and analyze/reconstruct the head motion
ds.eyePoseModelViewDataRight = []; % may as well save the model view matrix data as well - the hope is that this covers all of the critical information to later go back and analyze/reconstruct the head motion

Screen('TextSize', ds.w, 18); % Set default textsize
Screen('BlendFunction', ds.w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Setup the OpenGL rendering context of the onscreen window for use by
% OpenGL wrapper. After this command, all following OpenGL commands will
% draw into the onscreen window 'ds.w':
Screen('BeginOpenGL', ds.w);

glDisable(GL.LIGHTING);
glEnable(GL.DEPTH_TEST);
glDepthFunc(GL.LEQUAL); % From NeHe, specify what kind of depth testing

glEnable(GL.LINE_SMOOTH);
glEnable(GL.POINT_SMOOTH);
glEnable(GL.POLYGON_SMOOTH);

glHint(GL.POINT_SMOOTH_HINT, GL.NICEST);
glHint(GL.LINE_SMOOTH_HINT, GL.NICEST);
glHint(GL.POLYGON_SMOOTH_HINT, GL.NICEST);

glEnable(GL.FRAMEBUFFER_SRGB); % for proper gamma correction as per Mario Kleiner's advice

% Set viewport properly:
glViewport(0, 0, RectWidth(ds.windowRect), RectHeight(ds.windowRect));  % this is how the viewport is specified in all of the demos, but what it does it makes the horizontal dimension the shorter one and the vertical dimension the longer one

% Enable alpha-blending for smooth dot drawing:
glEnable(GL.BLEND);
glBlendFunc(GL.SRC_ALPHA, GL.ONE_MINUS_SRC_ALPHA);

% Set projection matrix: This defines a perspective projection,
% corresponding to the model of a pin-hole camera - which is a good
% approximation of the human eye and of standard real world cameras --
% well, the best aproximation one can do with 3 lines of code ;-)
% glMatrixMode(GL.PROJECTION);

% % Retrieve and set camera projection matrix for optimal rendering on the HMD:
if ~isempty(ds.hmd)
    [ds.projMatrix{1}, ds.projMatrix{2}] = PsychVRHMD('GetStaticRenderParameters', ds.hmd);
else 
    temp = load('dsHmdInfoCV1.mat');
    ds.projMatrix = temp.ds.projMatrix;
    clear temp    
end

% Finish OpenGL rendering into PTB window. This will switch back to the
% standard 2D drawing functions of Screen and will check for OpenGL errors.
Screen('EndOpenGL', ds.w);
