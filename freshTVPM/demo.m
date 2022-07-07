function VRHMDDemo(stereoscopic, checkerboard, deviceindex)
% 
% VRHMDDemo([stereoscopic=1][, checkerboard=0][, deviceindex=0])
%
% A very basic demo for the most basic setup of
% VR HMDs, e.g., the Oculus VR Rift DK2. It shows the
% absolute minimum of steps needed - one line of code - to
% use the first connected HMD as mono or stereo display.
%
% 'stereoscopic' if set to 1 (which is the default), configures the
% HMD as a stereoscopic display. A setting of 0 configures it as a
% monoscopic display.
%
% 'deviceindex' if provided, selects the HMD with given index. Otherwise
% the first HMD (deviceindex 0) is chosen.
%
% The demo just renders one static simple 2D image, or image
% pair in stereoscopic mode, then displays it in a loop until a
% key is pressed.

% History:
% 05-Sep-2015  mk Written.
% 30-Mar-2017  mk Adapt to the reality of new VR runtimes.

% Setup unified keymapping and unit color range:
PsychDefaultSetup(2);

if nargin < 1 || isempty(stereoscopic)
  stereoscopic = 1;
end

if nargin < 2 || isempty(checkerboard)
  checkerboard = 0;
end

if nargin < 3
  deviceindex = [];
end

% Select screen with highest id as Oculus output display:
screenid = max(Screen('Screens'));

% Open our fullscreen onscreen window with black background clear color:
PsychImaging('PrepareConfiguration');
if ~stereoscopic
  % Setup the HMD to act as a regular "monoscopic" display monitor
  % by displaying the same image to both eyes:
  PsychVRHMD('AutoSetupHMD', 'Monoscopic', 'LowPersistence FastResponse DebugDisplay', [], [], deviceindex);
else
  % Setup for stereoscopic presentation:
  PsychVRHMD('AutoSetupHMD', 'Stereoscopic', 'LowPersistence FastResponse', [], [], deviceindex);
end

[win, rect] = PsychImaging('OpenWindow', screenid);

Screen('BeginOpenGL', win);
  % Set viewport properly:
    glViewport(0, 0, RectWidth(winRect), RectHeight(winRect));
    
    % Retrieve and set camera projection matrix for optimal rendering on the HMD:
    %{
        ‘projL’ is the 4x4 OpenGL projection matrix for the left eye rendering.
        ‘projR’ is the 4x4 OpenGL projection matrix for the right eye rendering.
        Please note that projL and projR are usually identical for typical rendering
        scenarios.
    %}
    %projL          %projR (I guess)
    [projMatrix{1}, projMatrix{2}] = PsychVRHMD('GetStaticRenderParameters', hmd);
    
    % Setup modelview matrix: This defines the position, orientation and
    % looking direction of the virtual camera:
    glMatrixMode(GL.MODELVIEW);
    glLoadIdentity;
    
    % Set background clear color to 'black' (R,G,B,A)=(0,0,0,0):
    glClearColor(0,0,0,0);
    
    % Clear out the backbuffer: This also cleans the depth-buffer for
    % proper occlusion handling: You need to glClear the depth buffer whenever
    % you redraw your scene, e.g., in an animation loop. Otherwise occlusion
    % handling will screw up in funny ways...
    glClear;
    
centerDotxyz   = [0;0;0];
centerDotSize  = 10;
centerDotColor = [1, 1, 1]; 

drawCenterDot = glGenLists(1);
glNewList(drawCenterDot, GL.COMPILE);
moglDrawDots3D(win, centerDotxyz, centerDotSize, centerDotColor, [], 1);
glEndList;


% Render one view for each eye in stereoscopic mode:
vbl = [];
while ~KbCheck
  state = PsychVRHMD('PrepareRender', hmd);
  for renderPass = 0:1
    modelView = state.modelView{renderPass + 1};
            
    % Set per-eye projection matrix: This defines a perspective projection,
    % corresponding to the model of a pin-hole camera - which is a good
    % approximation of the human eye and of standard real world cameras --
    % well, the best aproximation one can do with 2 lines of code ;-)
    glMatrixMode(GL.PROJECTION);
    glLoadMatrixd(projMatrix{renderPass+1});

    % Setup camera position and orientation for this eyes view:
    glMatrixMode(GL.MODELVIEW);
    glLoadMatrixd(modelView);
    glClear;

    glCallList(drawCenterDot);
    
    Screen('EndOpenGL', win);
    vbl(end+1) = Screen('Flip', win);
  end
end

KbStrokeWait;
sca;

close all;
plot(1000 * diff(vbl));

end