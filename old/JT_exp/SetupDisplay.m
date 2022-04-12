function [ds,oc] = SetupDisplay(ds)
% CSB: 08/01/2019. ATTN, desparately need to clean up this function and remove
% unnecessary things, not hardcode whether oculus is connected, etc.  


% Setup Psychtoolbox for OpenGL 3D rendering support and initialize the
% mogl OpenGL for Matlab/Octave wrapper:
InitializeMatlabOpenGL(1);
PsychDefaultSetup(2); % the input of 2 means: execute the AssertOpenGL command, execute KbName('UnifyKeyNames') routine, AND unifies the color mode and switches from 0-255 to 0.0-1.0 - color part only impacts the current function or script, not ones that are called
% CSB: 08/1/2019, didn't we already do this in runExperiment.m??

ds.oculusConnected = 1; %0 % Is the HMD connected
ds.screenId = max(Screen('Screens')); % Find the screen to use for display:
ds.multiSample = 8;
ds.doSeparateEyeRender = 1; % render two eyes' views
PsychImaging('PrepareConfiguration');

% even in the 'fixed' viewing condition, we still want to track the head
% movement to save it out later - we just don't update the display in
% response to the movement
if ds.oculusConnected==1
    ds.hmd = PsychVRHMD('AutoSetupHMD', 'Tracked3DVR', 'LowPersistence TimeWarp FastResponse DebugDisplay', 0);
    PsychVRHMD('SetHSWDisplayDismiss', ds.hmd, -1);
    load('DefaultHMDParameters.mat');
    oc.defaultState = defaultState;
    
    % as of the PTB release for the CV1, there is built-in gamma correction
    % in Psychtoolbox for the devices
    
    %     load OculusDK2Gamma.mat; % load the oculus gamma table
    % for gamma correction
    %     ds.gammaVals = [GammaValue GammaValue GammaValue];
    %     PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma'); % gamma correction
else % oculus not connected
    ds.hmd = [];
    fprintf('No VR-HMD available, using default values.\n');
    load('OculusDK2Gamma.mat'); % load the oculus gamma table
    % for gamma correction
    ds.gammaVals = [GammaValue GammaValue GammaValue];
    load('DefaultHMDParameters.mat');
    oc.defaultState = defaultState;
    
    %PsychImaging('AddTask', 'General', 'SideBySideCompressedStereo'); % CSB: for specifying which eye channel we're displaying to. This isn't working with Screen('SelectStereoDrawBuffer')...
end


if ~isempty(ds.hmd) % Oculus connected
    [ds.w, ds.windowRect] = PsychImaging('OpenWindow', ds.screenId, 0, [], [], [], [], ds.multiSample);  % keeps automatically setting the stereo mode to 6 for the oculus - this is because, as indicated in MorphDemo.m: "% Fake some stereomode if HMD is used, to trigger stereo rendering"
else % Oculus not connected
    if ds.screenId % more than one screen connected, run fullscreen
        ds.winRect = [];
    else
        ds.winRect = [0 0 800 400];
    end
    [ds.w, ds.windowRect] = PsychImaging('OpenWindow', ds.screenId, 0, ds.winRect, [], [], 4, ds.multiSample);  % CSB. split screen mode; "4" sets this
end
%
% if ~DEBUG_FLAG
%     PsychColorCorrection('SetEncodingGamma', ds.w, 1./ds.gammaVals); % set required gamma
% end

% % Query info about this HMD: - note that this information has been saved in
% % the file HMDInfo.m in the Working Demo Files folder
if ~isempty(ds.hmd)
    ds.hmdinfo = PsychVRHMD('GetInfo', ds.hmd); % verifies that all of the basic requirements have been set up.
end

switch ds.experimentType
    case {'real'}
        ds.trackingFlag = 1; % screen will update with head movements
        ds.real = 1;
        ds.simulated = 0;
        ds.stabilized = 0;
    case {'simulated'}
        ds.trackingFlag = 1; % screen will update with head movements
        ds.real = 0;
        ds.simulated = 1;
        ds.stabilized = 0;
    case ('stabilized')
        ds.trackingFlag = 0;
        ds.real = 0;
        ds.simulated = 0;
        ds.stabilized = 1;
end


[ds.xc, ds.yc] = RectCenter(ds.windowRect);

if ~isempty(ds.hmd) % CSB: if using hmd
    ds.Height = 1.7614; % virtual height of the surround texture in meters, based on viewing distance - we want this to relate to the shorter dimension of the display
    ds.halfHeight = ds.Height/2;
    ds.Width = 1.7614; % virtual width of the surround texture in meters, based on viewing distance - we want this to relate to the longer dimension of the display
    ds.halfWidth = ds.Width/2;
    ds.yc = RectHeight(ds.windowRect)/2; % the horizontal center of the display in pixels
    ds.xc = RectWidth(ds.windowRect)/2; % the vertical center of the display in pixels
    ds.textCoords = [ds.yc ds.xc];
    
    if strcmp(ds.hmdinfo.modelName, 'Oculus Rift CV1')
        ds.hFOV = 80; % in deg - this is what is spit back from the Oculus readings at the start - horizontal field of view
        ds.vFOV = 90;  % in deg - vertical field of view
        ds.dFOV = sqrt(ds.hFOV^2 + ds.vFOV^2);
        ds.viewportWidthM = .09;% height component of the display - corresponding to the longer side of the display
        ds.metersPerDegree =  ds.viewportWidthM/ds.hFOV;
        ds.degreesPerM = ds.hFOV/ds.viewportWidthM;
        ds.viewportWidthDeg = ds.hFOV;
        ds.pixelsPerDegree = sqrt(RectHeight(ds.windowRect)^2 + RectWidth(ds.windowRect)^2) / ds.viewportWidthDeg;
        ds.pixelsPerM = sqrt(RectHeight(ds.windowRect)^2 + RectWidth(ds.windowRect)^2) / ds.viewportWidthM;
        ds.frameRate = 90;
    elseif strcmp(ds.hmdinfo.modelName, 'Oculus Rift DK2')
        ds.hFOV = 74;  % in deg - this is what is spit back from the Oculus readings at the start - horizontal field of view - but the display seems to flip them so h = v and vice versa
        ds.vFOV = 54;  % in deg - vertical field of view
        ds.dFOV = sqrt(ds.hFOV^2 + ds.vFOV^2);
        ds.viewportWidthM = 0.1512; % based on spec of 151.2 mm - height component of the display - corresponding to the longer side of the display
        ds.metersPerDegree =  ds.viewportWidthM/ds.hFOV;
        ds.viewportWidthDeg = ds.hFOV;
        ds.pixelsPerDegree = sqrt(RectHeight(ds.windowRect)^2 + RectWidth(ds.windowRect)^2) / ds.viewportWidthDeg;
        ds.pixelsPerM = sqrt(RectHeight(ds.windowRect)^2 + RectWidth(ds.windowRect)^2) / ds.viewportWidthM;
        ds.frameRate = 75;
    end
else % No hmd
    % CSB: screen params for default monitor, "whistler"
    ds.viewingDistance = .80; % CSB: viewing dist in m. for Heeger lab computer "whistler."
    ds.xyM = [.5 .31]; % width x height of "whistler" display in m.
    ds.xyPix = [RectWidth(ds.windowRect) RectHeight(ds.windowRect)]; % width by height of "whistler" display in pixels
    ds.xyDva = 2*atand(ds.xyM./(2.*ds.viewingDistance)); %  width by height of "whistler" display in degrees of visual angle
    ds.hFOV = ds.xyDva(1);  % in deg - horizontal field of view. % CSB: redundant, but one of these variables is being used elsewhere so I had to reproduce :/
    ds.vFOV = ds.xyDva(2);  % in deg - vertical field of view. % CSB: redundant, but one of these variables is being used elsewhere so I had to reproduce :/
    ds.viewportWidthM = ds.xyM(1);  % CSB: redundant, but one of these variables is being used elsewhere so I had to reproduce :/
    
    ds.dFOV = sqrt(ds.xyDva(1).^2 + ds.xyDva(2).^2);
    
    ds.Height = 1.7614; % virtual height of the surround texture in meters, based on viewing distance - we want this to relate to the shorter dimension of the display
    ds.halfHeight = ds.Height/2;
    ds.Width = 1.7614; % virtual width of the surround texture in meters, based on viewing distance - we want this to relate to the longer dimension of the display
    ds.halfWidth = ds.Width/2;
    
    ds.xyc = ds.xyPix./2; % the horizontal and vertical centers of the display in pixels
    ds.xc = ds.xyc(1); % the horizontal center of the display in pixels % CSB: redundant, but one of these variables is being used elsewhere so I had to reproduce :/
    ds.yc = ds.xyc(2); % the vertical center of the display in pixels % CSB: redundant, but one of these variables is being used elsewhere so I had to reproduce :/
    ds.textCoords = [ds.xc ds.yc];
    
    ds.metersPerDegree =  ds.viewportWidthM/ds.hFOV;
    ds.viewportWidthDeg = ds.hFOV;
    
    ds.pixelsPerDegree = sqrt(ds.xyPix(1).^2 + ds.xyPix(2).^2) ./ ds.xyDva(1);
    ds.pixelsPerM = sqrt(ds.xyPix(1).^2 + ds.xyPix(2).^2) ./ ds.xyM(1) ;
    
end


Screen('TextSize', ds.w, 18); % Size of text
Screen('TextStyle',ds.w,1); % 1=bold,2=italic
Screen('TextColor',ds.w,[255 255 255]);

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

% for proper gamma correction as per Mario Kleiner's advice
glEnable(GL.FRAMEBUFFER_SRGB);

% Set viewport properly:
glViewport(0, 0, RectWidth(ds.windowRect), RectHeight(ds.windowRect));  % this is how the viewport is specified in all of the demos, but what it does it makes the horizontal dimension the shorter one and the vertical dimension the longer one

% Enable alpha-blending for smooth dot drawing:
glEnable(GL.BLEND);
glBlendFunc(GL.SRC_ALPHA, GL.ONE_MINUS_SRC_ALPHA);

% Set projection matrix: This defines a perspective projection,
% corresponding to the model of a pin-hole camera - which is a good
% approximation of the human eye and of standard real world cameras --
% well, the best aproximation one can do with 3 lines of code ;-)
glMatrixMode(GL.PROJECTION);

% Retrieve and set camera projection matrix for optimal rendering on the HMD:
if ~isempty(ds.hmd)
    [ds.projMatrix{1} ds.projMatrix{2}] = PsychVRHMD('GetStaticRenderParameters', ds.hmd);%, 0.01, 5);  % add here the clipping plane distances; they are [clipNear=0.01],[clipFar=10000] by default
else
    
    ds.projMatrix{1} = [1.1903         0   -0.1486         0
        0    0.9998   -0.1107         0
        0         0   -1.0000   -0.0200
        0         0   -1.0000         0];
    
    ds.projMatrix{2} = [1.1903         0   0.1486         0
        0    0.9998   -0.1107         0
        0         0   -1.0000   -0.0200
        0         0   -1.0000         0];
    
end

% Initialize oculus modelview for head motion tracking
oc.modelViewDataLeft = []; % may as well save the model view matrix data as well - the hope is that this covers all of the critical information to later go back and analyze/reconstruct the head motion
oc.modelViewDataRight = []; % may as well save the model view matrix data as well - the hope is that this covers all of the critical information to later go back and analyze/reconstruct the head motion

% glLoadMatrixd(projMatrix);

% Setup modelview matrix: This defines the position, orientation and
% looking direction of the virtual camera:
glMatrixMode(GL.MODELVIEW);
glLoadIdentity;

% Set background clear color
glClearColor(0.5,0.5,0.5,1); % mid-gray
% glClearColor(0,0,0,1); % black

% Clear out the backbuffer: This also cleans the depth-buffer for
% proper occlusion handling: You need to glClear the depth buffer whenever
% you redraw your scene, e.g., in an animation loop. Otherwise occlusion
% handling will screw up in funny ways...
glClear;

% Finish OpenGL rendering into PTB window. This will switch back to the
% standard 2D drawing functions of Screen and will check for OpenGL errors.
Screen('EndOpenGL', ds.w);

ds.ifi = Screen('GetFlipInterval', ds.w); % Get duration of a single frame