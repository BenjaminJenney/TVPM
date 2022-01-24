function RunTVPMExperiment(subjectInitials,condition,monocularFlag,trainingFlag,mask,planes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors:     Charlie Burlingham and Jonathan Trattner, based on code by
%              Jackie Fulvio and Bas Rokers
% Last Update: 08/01/2019
% Purpose:     Measure precision and accuracy of heading perception in the
%              presence of real, simulated, or stabilized head rotation.
% Usage:       Give function participant initials, condition, and whether
%              this will be a training session (with auditory feedback)
%
%% IMPORTANT!  The Oculus must be plugged in and turned on *before* starting MATLAB
% The individual must be wearing the device prior to starting the
% experimental program to achieve proper frame rate
% Also the Oculus VR runtime version '0.5.0.1' *must* be installed for PTB
% to properly interact with/recognize the Oculus DK2, otherwise, use the
% latest version of the runtime for the CV1 (note: cannot have 2 versions
% on one machine)

%% Clean up workspace
close all;
sca;

%% Init input variables

pa.subjectName = subjectInitials; % initials of participant
ds.experimentType = condition; % 'real', 'simulated', 'stabilized'
ds.monocularFlag = monocularFlag; % 1: monocular viewing, 2: binocular viewing
pa.TRAINING = trainingFlag; % 1 for training runs (will run two runs without rotation), or 0 for experiment (18 runs, all rotation velocities and headings)

global DEBUG_FLAG MONOCULAR AUDIO
DEBUG_FLAG = 1; % 1: if you're debugging and want to print info to command line, and to show display on screen at half opacity ,  0: turn off those things
MONOCULAR = ds.monocularFlag ;

if pa.TRAINING == 1
    AUDIO = 1;
else
    AUDIO = 0;
end
AUDIO = 1; % to override and have feedback even when you're not training

if DEBUG_FLAG == 1
    %Screen('Preference', 'SkipSyncTests', 1); % For debugging, old
    PsychDebugWindowConfiguration([],.4); % display on monitor at half opacity in addition to in the HMD
end

%% Setup Psychtoolbox for OpenGL 3D rendering support
% and initialize the mogl OpenGL for Matlab/Octave wrapper:
global GL; % GL data structure needed for all OpenGL programs
InitializeMatlabOpenGL(1);
PsychDefaultSetup(2); % the input of 2 means: execute the AssertOpenGL command, execute KbName('UnifyKeyNames') routine, AND unifies the color mode and switches from 0-255 to 0.0-1.0 - color part only impacts the current function or script, not ones that are called

%% Initialize screen, experimental parameters, and keyboard
% [ds,oc] = SetupDisplay(ds); %(oc); % set up the display, based on the DK2
% [ds,pa] = SetupParameters(ds,pa); % set up the experimental parameters for this session
% shape   = SetupShapes(ds,pa);
% kb = SetupKeyboard(ds,DEBUG_FLAG); % get the keyboard info for the participant's responses
% pa.runNumber = 1;


% [hGratings, yPhaseChange, xPhaseChange] = Preprocess(ds, pa, kb, oc, shape, GL);

[ds,oc] = SetupDisplay(ds); %(oc); % set up the display, based on the DK2
[ds,pa] = SetupParameters(ds,pa); % set up the experimental parameters for this session
shape   = SetupShapes(ds,pa);
kb = SetupKeyboard(ds,DEBUG_FLAG); % get the keyboard info for the participant's responses
ds.vbl = Screen('Flip', ds.w);
%% Start the experiment - opening screen, getting the participant set

pa.runNumber = 1; % initialize the run counter
pa.readyToBeginRun = 0; % initialize the ready-to-begin-a-new-run variable
readyToBegin = 0;
Screen('BeginOpenGL', ds.w);
                            sf =.1; % cycles per deg
    sf = (ds.deg_per_px).*sf; %1/px_per_deg %convert from deg to 
    af = 2*pi*sf;
    xGrating = (.5 + .5 * sin(af * shape.disk.texture.x - 0)); % vertically oriented grating
    yGrating = (.5 + .5 * cos(af * shape.disk.texture.y - 0)); % horiz grating
    plaidData = 255.*((xGrating + yGrating)/2);
    plaidData = repmat(plaidData,[ 1 1 3 ]);
    %x = single(plaidData);
    plaidData = permute(plaidData,[ 3 2 1 ]);

type = GL_TEXTURE_2D;
s1 = 512; s2 = 512;
sigma = 360;
[x,y] = meshgrid(-1*(s1/2):(s1/2)-1, -1*(s2/2): (s2/2)-1);
pa.rmin_bg = 45.6874;% pixels 
pa.rmax_bg = 137.7631;% pixels  
pa.rstrip = 11.6268;% this cuts a strip into the fixation disk that has a height the size of the paddle height
aperture = zeros(s1,s2).*255;
opaque = ones(size(x'));
opaque = min(opaque, ((sqrt((x').^2+(y').^2) > pa.rmax_bg))); % | ((abs(y') > pa.rstrip) & sqrt((x').^2+(y').^2) < pa.rmin_bg)));

gauss = ((1-exp(-((x/sigma).^2)-((y/sigma).^2)))).*255;
%keyboard
 aperture = repmat(aperture,[ 1 1 4 ]);
 aperture = permute(aperture,[ 3 2 1 ]);
 aperture(4,:,:) = shiftdim(255 .* opaque, -1);
  
     % keyboard
textureid = glGenTextures(2);
glBindTexture(type, textureid(1));
glTexParameterfv(type, GL.TEXTURE_WRAP_S, GL.REPEAT);
glTexParameterfv(type, GL.TEXTURE_WRAP_T, GL.REPEAT);
glTexParameterfv(type, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
glTexParameterfv(type, GL.TEXTURE_MIN_FILTER, GL.LINEAR);
glTexImage2D(type, 0, GL.RGB, 256, 256, 0, GL.RGB, GL.UNSIGNED_BYTE, uint8(plaidData));
glTexEnvi(GL.POINT_SPRITE, GL.COORD_REPLACE, GL.TRUE);
glBindTexture(type,0);

glBindTexture(type, textureid(2));
glTexParameterfv(type, GL.TEXTURE_WRAP_S, GL.REPEAT);
glTexParameterfv(type, GL.TEXTURE_WRAP_T, GL.REPEAT);
glTexParameterfv(type, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
glTexParameterfv(type, GL.TEXTURE_MIN_FILTER, GL.LINEAR); 
glTexImage2D(type, 0, GL.RGBA, s1, s2, 0, GL.RGBA, GL.UNSIGNED_BYTE, uint8(aperture));
glTexEnvi(GL.POINT_SPRITE, GL.COORD_REPLACE, GL.TRUE);
glBindTexture(type,0)
    %keyboard
 Screen('EndOpenGL', ds.w);
while ~readyToBegin % confirm everything's ready to go
    
    % Camera position when using head tracking + HMD: (according to SuperShapeDemo.m)
    globalPos = [0, 0, 0]; % x,y,z  % in meters - this is the starting position of head... should be at origin
    heading = 0; % yaw.
    ds.globalHeadPose = PsychGetPositionYawMatrix(globalPos, heading); % just an identity matrix to start with % initialize observer's start position to the default camera position specified above
    
    if isempty(ds.hmd) % Oculus not connected
        load DefaultHMDParameters.mat;
        oc.defaultState = defaultState;
        oc.initialState = defaultState.initialState;
    else % Oculus connected
        oc.initialState = PsychVRHMD('PrepareRender', ds.hmd, ds.globalHeadPose);  % get the state of the hmd now
    end
    
    instructionsCenterxL = ds.textCoords(2)-200+(pa.instructionDisparity/2);
    instructionsCenterxR = ds.textCoords(2)-200-(pa.instructionDisparity/2);
    
    for renderPass = 0:1 % loop over eyes
        if MONOCULAR && renderPass == 0 || ~MONOCULAR
            ds.renderPass = renderPass;
            Screen('SelectStereoDrawBuffer',ds.w,ds.renderPass);
            Screen('BeginOpenGL',ds.w);
            
            % Setup camera position and orientation for this eyes view:
            glMatrixMode(GL.PROJECTION)
            glLoadMatrixd(ds.projMatrix{renderPass + 1});
            
            modelView = oc.initialState.modelView{ds.renderPass + 1}; % Use per-eye modelView matrices
            
            % N/A, redundant compared to modelView above. pa.modelViewSaveOutForFixed{ds.renderPass + 1} =  oc.initialState.modelView{ds.renderPass + 1};
            
            glLoadMatrixd(modelView);
            
            Screen('EndOpenGL', ds.w);
            if renderPass == 0
                Screen('DrawText',ds.w,'Ready to start the experiment?',instructionsCenterxL,ds.textCoords(1)-200,[1 1 1]);
                Screen('DrawText',ds.w,'Press SPACE to confirm.',instructionsCenterxL,ds.textCoords(1)+50-200,[1 1 1]);
            elseif renderPass == 1
                Screen('DrawText',ds.w,'Ready to start the experiment?',instructionsCenterxR,ds.textCoords(1)-200,[1 1 1]);
                Screen('DrawText',ds.w,'Press SPACE to confirm.',instructionsCenterxR,ds.textCoords(1)+50-200,[1 1 1]);
            end
        end
    end
    frameiter = 0;
    Screen('DrawingFinished', ds.w);
    ds.vbl = Screen('Flip', ds.w);
    
    [kb.keyIsDown, kb.secs, kb.keyCode] = KbCheck(-1); % query the keyboard
    if kb.keyIsDown && kb.keyCode(kb.spacebarKey)
        readyToBegin=1;
    end
end


%% Participant is ready, so let's start first run of trials

ds.tElapsed = 0;
ds.fCount = 1;

[ds, pa, kb] = SetupNewTrial(ds, pa, kb);
ds.vbl = pa.trialOnset;
tStart = ds.vbl;
pa.experimentOnset = ds.vbl;

pa.block = 0;
kb.nextTrialKey = 0;

rotationMatrixYawHomo = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; % initialize for simulated condition
%[hGratings, hPhaseChange, vPhaseChange] = Preprocess(ds,pa,shape,GL);

while (pa.runNumber <= pa.numRuns) && ~kb.keyCode(kb.escapeKey)
    [vXw, vYw] = Preprocess(ds, pa, shape, GL);% m/f in world coords
    for i = 1:shape.plane.numPlanes
        posX{i} = shape.disk.X_m{i} + squeeze(vXw(1,1,1,:));
        posY{i} = shape.disk.Y_m{i} + squeeze(vYw(1,1,1,:));
        %posZ{i} = ones(size(posX{i})) .* shape.plane.depths_m(i); 
    end
    %k = 1;
    while (pa.trialNumber < pa.nTrials) && ~kb.keyCode(kb.escapeKey) % wait until all of the trials have been completed or the escape key is pressed to quit out
        
            
        
        % Get HMD state
        if isempty(ds.hmd) % Oculus is not connected - will display a poor imitation of the Oculus rift on your main computer screen
            if ~exist('state')
                state = oc.defaultState; % just set to a default, non-updating viewpoint
            end
        else   % Oculus is connected - uses PTB's code + openGL code to display in the HMD
            % Track and predict head position and orientation, retrieve modelview
            % camera matrices for rendering of each eye. Apply some global transformation
            % to returned camera matrices. In this case a translation + rotation, as defined
            % by the PsychGetPositionYawMatrix() helper function:
            state = PsychVRHMD('PrepareRender', ds.hmd, ds.globalHeadPose);  % Mark the start of the rendering cycle for a new 3D rendered stereoframe. Return a struct 'state' which contains various useful bits of information for 3D stereoscopic rendering of a scene, based on head tracking data
        end
        
        % Render the scene separately for each eye:
        
        
        for renderPass = 0:1 %0 left, 1 right eye
            
            if (MONOCULAR && renderPass == 0) || (~MONOCULAR)
                ds.renderPass = renderPass;
                if isempty(ds.hmd) % hmd not connected
                    % Get head position from keyboard input
                    
                    eye.eyeIndex = ds.renderPass; % We are switching eye index, but not the eye.modelview here
                    % CSB: this inits camera matrix if it is empty
                    % BR: This should happen outside of the loop init'd on line 131
                    if ~isfield(eye,'modelView') % checks if camera matrix exists. it shouldn't on first trial and it will be init'd.
                        if ds.renderPass==0 % drawing left eye
                            eye.modelView = oc.defaultState.modelViewDataLeft;
                        elseif ds.renderPass==1 % drawing right eye
                            eye.modelView =  oc.defaultState.modelViewDataRight;
                            % eye.modelView(1,4) =  eye.modelView(1,4)+100; % CSB: debug
                        end
                    end
                    
                    %[pa, kb, eye, ds, state]  = GetKeyboardHeadmotion(pa, ds, kb, eye, state);  % query the keyboard to allow the observer to rotate the paddle and eventually lock in his/her response to initiate a new trial
                    
                else % hmd connected
                    % Query which eye to render in this ds.renderPass, and query its
                    % eyePose vector for the predicted eye position to use for the virtual
                    % camera rendering that eyes view. The returned pose vector actually
                    % describes tracked head pose, ie. HMD position and orientation in space.
                    eye = PsychVRHMD('GetEyePose', ds.hmd, ds.renderPass, ds.globalHeadPose);
                    headRotationAloneInfo = eye.modelView; % this just tells you how head has rotated
                    
                    
                    % this is for saving purposes to recreate participants' head motion
                    if ds.renderPass % drawing right eye
                        oc.modelViewDataRight = [oc.modelViewDataRight; eye.modelView];
                    else % drawing left eye
                        oc.modelViewDataLeft = [oc.modelViewDataLeft; eye.modelView];
                    end
                    
                    if ~ds.trackingFlag % 'stabilized' condition loop - don't update the scene with tracked head motion, just use the default state
                        eye.modelView = modelView; %pa.modelViewSaveOutForFixed{ds.renderPass + 1}; % oc.defaultState.initialState.modelView{ds.renderPass + 1}; % comes back from the initial call...will not update the scene based on head tracking
                        state.tracked = 2;
                        eye.eyeIndex = ds.renderPass;
                    end
                    
                end
                
                Screen('SelectStereoDrawbuffer', ds.w, eye.eyeIndex); % Select 'eyeIndex' to render (left- or right-eye):
                
                if ds.simulated
                    modelView = rotationMatrixYawHomo*eye.modelView;
                else
                    modelView = eye.modelView; % Extract modelView matrix for this eye:
                end
                %Screen('DrawTextures', ds.fullWindowMask, ds.masktex, [], ds.dstRects); % see setupDisplay for rest of code doing this
                Screen('BeginOpenGL', ds.w); % Manually reenable 3D mode in preparation of eye draw cycle
                
                % Setup camera position and orientation for this eyes view:
                glMatrixMode(GL.PROJECTION)
                glLoadMatrixd(ds.projMatrix{renderPass + 1});
                
                glMatrixMode(GL.MODELVIEW);
                glLoadMatrixd(modelView); % updates openGL camera
                
                glClearColor(0, 0, 0, 1); % gray background
                glClear(); % clear the buffers - must be done for every frame
                glColor3f(1,1,1);
                
                glPushMatrix;
                %% Experiment Logic
                if ds.vbl <  pa.trialOnset + pa.targetMotionDuration % if current time < present until time: draw target, 1 s target motion
                    
                    %disp(['frame = ' num2str(pa.frame)])
                    %disp(['vbl = ' num2str(ds.vbl)])
                    
                    pa.flippy = 1; % binary variable for presenting fixation in simulated condition in a single spot
                    
                    pa.RotationVel(pa.trialNumber) = pa.rotSpeedVecRandThisRun(pa.trialNumber); %current rotation velocity to be used in the staircase determination and speed setting
                    
                    for cc = 1:pa.numRotationVelocities  %numHeadingAngles only when in practice trial mode so there are less staircases made when analyzing the data. When in real mode, use for cc = 1:pa.numRotationVelocities % loop through rotation speeds
                        if pa.RotationVel(pa.trialNumber) == cc
                            if ds.real == 1
                                pa.rotVelocityRadPursuit(pa.trialNumber) = deg2rad(pa.rotationVelocities(cc));
                            elseif ds.simulated == 1
                                pa.rotVelocityRad = deg2rad(pa.rotationVelocities(cc));
                            end
                        end
                    end
                    
                    for dd = 1:pa.numHeadingAngles
                        if pa.headingAngleVecRandThisRun(pa.trialNumber) == dd
                            pa.heading(pa.trialNumber) = pa.headingAngles(dd);
                        end
                    end
                    
                    if ds.real == 1
                        pa.rotVelocityRad = deg2rad(0); % don't rotate world in real head rotation condition
                    elseif ds.simulated == 1
                        pa.rotVelocityRadPursuit(pa.trialNumber) = deg2rad(0); % don't rotate pursuit target in simulated head rotation condition
                    end
                    
                    elapsedTime = (ds.vbl-pa.trialOnset); % elapsed time in seconds since trial start (pressing space bar)
                    
                    % Calculate instantaneous delta t
                    %dt = GetSecs() - ds.vbl; % delta t in secs
                    
                    xDotDisplacement = pa.transSpeed.*sind(-1.*pa.heading(pa.trialNumber)).*elapsedTime; % CSB: should be in meters. check. % ATTN: take negative of heading angle, so it displays correctly. otherwise displayed heading is inverse of ground truth
                    % pa.rotVelocityRadPursuit(pa.trialNumber) = deg2radtrialNumber)).*elapsedTime; % CSB: should be in meters. check. % ATTN: take negative of heading angle, so it displays correctly. otherwise displayed heading is inverse of ground truth
                    zDotDisplacement = pa.transSpeed.*cosd(-1.*pa.heading(pa.trialNumber)).*elapsedTime;  % JT: make sure this trig is right.
                    
                    if ds.simulated
                        % make simulated rotation matrix into homogenous cooridnates and
                        % invert it to transform back to original space before simulated rotation ... and
                        % additionally multiply by inverse modelView to undo head tracking,
                        % This will keep the fixation point in the line of sight.
                        %rotationMatrixYawHomo = [[rotationMatrixYaw; 0 0 0] [0 0 0 1]'];
                        %                         fixationPos =  inv(modelView) * [pa.fixationVertexPos, 1.0]'; % for fixation, must rotate back to center of screen according to BOTH head movement AND world rotation
                        %                         fixationPos = fixationPos(1:3);
                        %                         %Get ready to make the pursuit target fixed, too.
                        %                         % OLD % pursuitTargetCoords = inv(rotationMatrixYawHomo) * [xPosPursuit, 0, -zPosPursuit 1]'; % fox pursuit, we want it fixed in the world, despite simulated rotation, so just rotate back according to simulated rotation
                        %                         pursuitTargetCoords = inv(rotationMatrixYawHomo) * [xPosPursuit, 0, -zPosPursuit 1]'; % for pursuit, we want it fixed in the world, despite simulated rotation, so just rotate back according to simulated rotation
                        
                        %Translating the pursuit target
                        %moglDrawDots3D(ds.w, pursuitTargetCoords(1:3), 3*2, pa.red, [], 2); % fixation csb june 28 2021
                        %moglDrawDots3D(ds.w, pursuitTargetCoords(1:3), 36, pa.black, [], 2);
                        
                        %moglDrawDots3D(ds.w, fixationPos, pa.fixationDiameter,[1 1 1 1], [], 2); %drawing horizontal fixation dot  % CSB june 21 2021 uncomment
                        
                    end
                    
                    
                    if ds.simulated && planes == 0
                        
                        glPushMatrix;
                        
                        %if sign(pa.heading(pa.trialNumber)) ==1
                        
                        gluLookAt(-xDotDisplacement,0,-zDotDisplacement,0,0,-pa.cubeWidth/4,0,1,0);
                        
                        %elseif sign(pa.heading(pa.trialNumber)) ==-1
                        %    gluLookAt(-xDotDisplacement - (.064/2),0,-zDotDisplacement,0,0,-pa.cubeWidth/4,0,1,0);
                        %end
                        
                        moglDrawDots3D(ds.w, pa.sphere, 3, [1 1 1 1], [], 2); %drawing sphere
                        
                        glPopMatrix; % we do not want the fixation dot to move with the world
                        
                        if mask == 1
                            
                            drawTheHoleyBag(ds, 1);
                            
                        end
                        
                        moglDrawDots3D(ds.w, [0 0 -pa.cubeWidth/4 1]', 3*2, pa.red, [], 2); %drawing fixation dot in environ % CSB june 21 2021 uncomment
                    elseif ds.real && planes == 0
                        %{
                        %                         glPushMatrix;
                        %                         %glTranslatef(xDotDisplacement,0,zDotDisplacement); % shift the dot world to its position along its trajectory for this frame
                        %                         glEnable(GL.POINT_SPRITE);
                        %                         shape.disk.textureArray(mi).bind;
                        %                         %shape.disk.textureArray(mi).draw(randPoint_test)
                        %
                        %                         moglDrawDots3D(ds.w, shape.disk.coords, 20, [1 1 1 1], [], []);
                        %
                        %                         shape.disk.textureArray(mi).unbind;
                        %                         glDisable(GL.POINT_SPRITE);
                        %                         mi = mi+1;
                        %                         glPopMatrix;
                        %}
                        if mask == 1, drawTheHoleyBag(ds, 1), end % draws the ptb mask. Params: drawTheHoleyBag(window,[openglEnbled = 1]);
                        
                        glPushMatrix;
                        glTranslatef(xDotDisplacement,0,zDotDisplacement);
                        
                        moglDrawDots3D(ds.w, [0 0 -pa.cubeWidth/4 1]', 3*2, pa.red, [], 2); %drawing fixation dot in environ % CSB june 21 2021 uncomment
                        glPopMatrix;
                        
                     elseif planes == 1
                        %{
                        %                            glPushMatrix;
                        %
                        %                            glLoadIdentity;
                        %                            preMV = glGetFloatv(GL.MODELVIEW_MATRIX);
                        %                            preMV = reshape(preMV, [4,4]);
                        %                            gluLookAt(-xDotDisplacement,0,-zDotDisplacement,0,0,-pa.cubeWidth/4,0,1,0);
                        %                            %glTranslatef(-ds.planeWidths_m(2) + xDotDisplacement, 0.0, 0.0 + zDotDisplacement);
                        %                            curMV = glGetFloatv(GL.MODELVIEW_MATRIX);
                        %
                        %                            curMV = reshape(curMV, [4,4]);
                        %
                        %     t = t+1;
                        %                          for i = 1:ds.numFrames
                        
                        %                            if t == 2; keyboard;end
                        %                           for i = 1:shape.plane.numPlanes
                        %                             opticFlow(ds, pa, shape.disk.xpos_deg(:,:,i), shape.disk.ypos_deg(:,:,i), shape.plane.planes(i), 1050, curMV, t, omegaY, T);
                        %                           end
                        %                           glPopMatrix;
                        %keyboard;
                       % tic
                      % if pa.frame <= pa.nFrames % not necessary bc we are running under the frame rate of monitor, but whatever.
                       % glEnable(GL.POINT_SPRITE);
%                        for i = 1:shape.plane.numPlanes                
%                            for j = 1:shape.disk.numDisks
%                                 vCurPhase = vPhaseChange(pa.trialNumber,pa.frame,renderPass+1,i,j); if vCurPhase == 0; vCurPhase = 1; end
%                                 hCurPhase = hPhaseChange(pa.trialNumber,pa.frame,renderPass+1,i,j); if hCurPhase == 0; hCurPhase = 1; end
%                                 image = .5.*hGratings{vCurPhase} + .5.*hGratings{hCurPhase}';
%                                 image = 255.*(image./2);
%                                 image = repmat(image,[ 1 1 3 ]);
%                                 image = permute(image,[ 3 2 1 ]);
%                                glBindTexture(GL_TEXTURE_2D, diskTexture.id)
%                                glTexSubImage2D(GL_TEXTURE_2D, 0,0,0, diskTexture.width, diskTexture.height, GL_RGB, GL_UNSIGNED_BYTE, uint8(image));
                                 %keyboard
                                 %glEnable(GL_DEPTH_TEST);
                       %}
                         glEnable(GL.POINT_SPRITE);
                                 
                                 
                                 for i=1:shape.plane.numPlanes
                                    %diskTexture.id = glGenTextures(1);
                                    %glEnable(GL.TEXTURE_2D); % Enable 2D texture mapping
                                    %glActiveTexture(GL_TEXTURE0);
                                    %DRAW DISKS
                                    glBindTexture( type, textureid(1))

                                    moglDrawDots3D(ds.w, [posX{i},posY{i},shape.disk.Z_m{i}]', shape.disk.size_px, [], [], 0);
                                    
                                    glBindTexture(type, 0);
                                    %END DRAW DISKS
                                     %APERTURES
%                                     glBindTexture(type, textureid(2));
%                                     
%                                     moglDrawDots3D(ds.w, [posX{i},posY{i}, shape.disk.Z_m{i} + .001]', shape.disk.size_px + 10, [], [], []);
% 
%                                     glBindTexture(type, 0);
                                     %END APERTURES
                                 end
                                 glDisable(GL.POINT_SPRITE);
                                 %glBindTexture(GL.TEXTURE_2D, 0);
                                %
                                %glBindTexture(GL.TEXTURE_2D, 0)
 %                           end
 %                       end
                       % glDisable(GL.POINT_SPRITE);
                      % end
                        %glCallList(ds.point)
                        %
                       % toc
                        %glCallList(ds.planes(1));
                        %                           glPopMatrix;
                        %
                        %                           glPushMatrix;
                        %                           glTranslatef(0.0 + xDotDisplacement, 0.0, 0.0 + zDotDisplacement);
                        %                           glCallList(ds.planes(2));
                        %                           glPopMatrix;
                        %
                        %                           glPushMatrix;
                        %                           glTranslatef(ds.planeWidths_m(3) + xDotDisplacement, 0.0, 0.0 + zDotDisplacement);
                        %                           glCallList(ds.planes(3));
                        %                           glPopMatrix;
                        
                        %if mask == 1, drawTheHoleyBag(ds, 1), end
                        
                    end
                    %if ds.simulated % added csb july 3 2021. for plotting original fixation pt, for debugging. the two fixation pts should be in the same retinal location
                    %                         glPushMatrix;
                    %                         moglDrawDots3D(ds.w, inv(rotationMatrixYawHomo)*[0 0 -pa.cubeWidth/4 1]', 3*2, pa.red, [], 2); %drawing fixation dot in environ % CSB june 21 2021 uncomment
                    %                         glPopMatrix;
                    %                     end
                    
                elseif ~kb.responseGiven && ds.vbl >  pa.trialOnset + pa.targetMotionDuration && ds.vbl <=  pa.trialOnset + pa.targetMotionDuration + pa.itiLength  % show paddle and allow observers to adjust its position - no time constraints - they press the space bar to lock in their response and start a new trial
                    
                    if ds.real && planes == 0
                        glPushMatrix;
                        glTranslatef(xDotDisplacement,0,zDotDisplacement); % shift the target to its position along its trajectory for this frame
                        moglDrawDots3D(ds.w, pa.sphere, 3, [1 1 1 1], [], 2); %drawing sphere
                        moglDrawDots3D(ds.w, [0 0 -pa.cubeWidth/4], 3*2, pa.red, [], 2); %drawing fixation dot in environ % CSB june 21 2021 uncomment
                        glPopMatrix;
                    elseif ds.simulated && planes == 0
                        glPushMatrix;
                        gluLookAt(-xDotDisplacement,0,-zDotDisplacement,0,0,-15,0,1,0)
                        moglDrawDots3D(ds.w, pa.sphere, 3, [1 1 1 1], [], 2); %drawing sphere
                        glPopMatrix;
                        moglDrawDots3D(ds.w, [0 0 -pa.cubeWidth/4 1]', 3*2, pa.red, [], 2); %drawing fixation dot in environ % CSB june 21 2021 uncomment
                    elseif planes == 1
                        
                        %                         glAlphaFunc(GL_GREATER, 0.5); %diskards px fragments below .5 alpha value
                        %                         glEnable(GL_ALPHA_TEST);      %This way we can see the post without render ordering problems
                        %
                        %\                          glPushMatrix;
                        %                         glTranslatef(-ds.planeWidths_m(1), 0.0, 0.0);
                        %                         glCallList(ds.planes(1));
                        %                         glPopMatrix;
                        %
                        %                         glPushMatrix;
                        %                         glTranslatef(0.0, 0.0, 0.0);
                        %                         glCallList(ds.planes(2));
                        %                         glPopMatrix;
                        %
                        %                         glPushMatrix;
                        %                         glTranslatef(ds.planeWidths_m(3), 0.0, 0.0);
                        %                         glCallList(ds.planes(3));
                        %                         glPopMatrix;
                        %
                        %                         glDisable(GL_ALPHA_TEST);
                        
                    end
                    
                    
                    pa.responseOnset = ds.vbl; % start the timer on the response time
                    
                    if ~pa.stairCaseCodeRanAlready(1) == 1 % so just run staircase code ONCE per trial and then reset dummy vairable in setupNewTrial
                        
                        % flip a coin to determine which staircase you use
                        for cc = 1:pa.numRotationVelocities % loop through rotation speeds
                            for dd = 1:pa.numHeadingAngles
                                if pa.RotationVel(pa.trialNumber) == cc && pa.headingAngleVecRandThisRun(pa.trialNumber) == dd
                                    if pa.coinFlip >= .5
                                        pa.whichStairCase(pa.trialNumber) = 1;  %  use staircase 1 for this rotation velocity
                                    elseif pa.coinFlip < .5
                                        pa.whichStairCase(pa.trialNumber) = 2;  %  use staircase 2 for this rotation velocity.
                                    end
                                end
                            end
                        end
                        
                        currStaircaseIndex = pa.whichStairCase(pa.trialNumber);
                        currVelocityIndex = pa.RotationVel(pa.trialNumber);
                        currHeadingIndex = pa.headingAngleVecRandThisRun(pa.trialNumber);
                        
                        pa.postPosDeg(pa.trialNumber) = pa.heading(pa.trialNumber) + pa.stairs(currVelocityIndex,currStaircaseIndex,currHeadingIndex).threshold; % csb june 21 2021. added true heading to staircase distance to get post pos in deg
                        
                        % Compute whether observer is heading left or right of the post in deg:
                        if pa.postPosDeg(pa.trialNumber) < pa.heading(pa.trialNumber)
                            pa.groundTruthLeftOrRight = 1; % right % if the post is to the left of the true heading, the correct answer is "right"
                        elseif pa.postPosDeg(pa.trialNumber) > pa.heading(pa.trialNumber)
                            pa.groundTruthLeftOrRight = 0; % left % if the post is to the right of the true heading, the correct answer is "left"
                        elseif pa.postPosDeg(pa.trialNumber) == pa.heading(pa.trialNumber) % if they are equal, we want the staircase to move up or down randomly
                            pa.groundTruthLeftOrRight = 777;
                        end
                        
                        % Drawing the post based on staircase
                        postRadius = abs(pa.pursuitTargetDist); % The post distance from you = pursuit target distance from you = fixation distance from you in virtual environment
                        postPosX = postRadius*sind(pa.postPosDeg(pa.trialNumber));
                        postPosZ =  postRadius*cosd(pa.postPosDeg(pa.trialNumber));
                        
                        % uncomment if you want to plot a test target which
                        % should appear in the line of sight at the end of
                        % the trial in simulated AND real head rotation
                        % conditions (given some rotation). this will
                        % verify that the coordinate system is rotating
                        % appropriately in simulated and that the post is
                        % in the correct current coordinate system
                        %{
                        testEndLocation = pa.targetMotionDuration*pa.rotationVelocities(pa.RotationVel(pa.trialNumber));
                        testPosX = postRadius*sind(testEndLocation);
                        testPosZ =  postRadius*cosd(testEndLocation);
                        %}
                        
                        pa.stairCaseCodeRanAlready(1) = 1; % will only run this code in conditional once
                    end
                    glPushMatrix
                    %f = inv(eye.modelView);
                    if ds.simulated
                        postCoords = [postPosX , 0, -postPosZ]'; %  csb june 28 2021 - this needs to stay like this, bc the heading shifts on the retina over time, and so should the post accordingly!. if coordinate system rotates, causing post to rotate with world!!! this is essential so people don't do curved path task
                        %OLD if coordinate system doesn't rotate %postCoords = inv(rotationMatrixYaw)*[postPosX , 0, -postPosZ]'; % does this do the right thing??? i'm just updating according to the head rotation info from beginning of frame... important to test explicitly somehow if the dot actually appears in a fixed straight ahead in world coords
                    else
                        postCoords = [postPosX , 0, -postPosZ]'; % does this do the right thing??? i'm just updating according to the head rotation info from beginning of frame... important to test explicitly somehow if the dot actually appears in a fixed straight ahead in world coords
                    end
                    moglDrawDots3D(ds.w, postCoords, pa.postDiameter, [1 0 1 1], [], 2); % CSB june 21 2021 uncomment
                    
                    % moglDrawDots3D(ds.w, [testPosX , 0, -testPosZ]', pa.postDiameter, [.7 0 .3 1], [], 2); % uncomment if you want to plot a test target which should appear in the line of sight at the end of the trial in simulated AND real head rotation conditions (given some rotation). this will verify that the coordinate system is rotating appropriately in simulated and that the post is in the correct current coordinate system
                    glPopMatrix;
                    
                    
                    [pa, kb] = GetResponse(pa, ds, kb);  % query the keyboard to allow the observer to rotate the paddle and eventually lock in his/her response to initiate a new trial
                    
                    pa.feedbackOnset = ds.vbl;
                    
                elseif kb.responseGiven && ds.vbl >  pa.trialOnset + pa.targetMotionDuration && ds.vbl <=  pa.trialOnset + pa.targetMotionDuration + pa.itiLength
                    
                    pa.leftRightResponse(pa.trialNumber) = kb.answer; % save out key press response of observer
                    
                    if ~pa.stairCaseCodeRanAlready(2) == 1 % will only run this code in conditional once
                        
                        if pa.groundTruthLeftOrRight == kb.answer % correct
                            pa.correctness(pa.trialNumber) = 1;
                        elseif pa.groundTruthLeftOrRight ~= kb.answer && kb.answer ~= -1
                            pa.correctness(pa.trialNumber) = 0;
                        elseif kb.answer == -1
                            pa.correctness(pa.trialNumber) = NaN; % so NaN means no response or WRONG key (not left or right) do nothing to staircase in that case.. O
                        elseif pa.groundTruthLeftOrRight == 777
                            pa.correctness(pa.trialNumber) = NaN; % if the true heading = the post position
                        end
                        
                        if ~isnan(pa.correctness(pa.trialNumber)) || pa.groundTruthLeftOrRight == 777 % if participant responded and didn't press wrong key, or the heading = post position, then update staircase. otherwise, don't update staircase
                            if pa.stairs(currVelocityIndex,currStaircaseIndex,currHeadingIndex).threshold >= pa.heading(pa.trialNumber)  % if the current staircase threshold is greater than the true heading , correct -> subtract step.. if if less than true heading, then correct -> ADD step
                                pa.stairs(currVelocityIndex,currStaircaseIndex,currHeadingIndex) = upDownStaircase(pa.stairs(currVelocityIndex,currStaircaseIndex,currHeadingIndex),pa.correctness(pa.trialNumber));
                            elseif pa.stairs(currVelocityIndex,currStaircaseIndex,currHeadingIndex).threshold < pa.heading(pa.trialNumber)
                                pa.stairs(currVelocityIndex,currStaircaseIndex,currHeadingIndex) = upDownStaircase(pa.stairs(currVelocityIndex,currStaircaseIndex,currHeadingIndex),~pa.correctness(pa.trialNumber)); % if you are on staircase, YOU MUST ADD A STEP SIZE WHEN YOU ARE CORRECT, NOT SUBTRACT, BC STAIRCASE STARTS AT NEG VAL!!! IMPORTANT.
                            elseif pa.stairs(currVelocityIndex,currStaircaseIndex,currHeadingIndex).threshold == pa.heading(pa.trialNumber) % if they're equal, flip a coin to determine direction
                                % flip a coin to determine which direction to go in if staircase = true heading
                                coin = rand() >= .5;
                                if coin
                                    pa.stairs(currVelocityIndex,currStaircaseIndex,currHeadingIndex) = upDownStaircase(pa.stairs(currVelocityIndex,currStaircaseIndex,currHeadingIndex),pa.correctness(pa.trialNumber));
                                else
                                    pa.stairs(currVelocityIndex,currStaircaseIndex,currHeadingIndex) = upDownStaircase(pa.stairs(currVelocityIndex,currStaircaseIndex,currHeadingIndex),~pa.correctness(pa.trialNumber));
                                end
                            end
                        end
                        
                        pa.stairCaseCodeRanAlready(2) = 1;
                        
                        if DEBUG_FLAG % print response info and heading and target location
                            disp(['post =' num2str(pa.postPosDeg(pa.trialNumber))])
                            disp(['heading =' num2str(pa.heading(pa.trialNumber))])
                            disp(['response =' num2str(pa.leftRightResponse(pa.trialNumber))]);
                            disp(['correctness =' num2str(pa.correctness(pa.trialNumber))])
                        end
                    end
                    
                    if AUDIO
                        % play appropriate sound
                        if ~pa.feedbackGiven == 1
                            if pa.correctness(pa.trialNumber) == 1
                                % PsychPortAudio('Start', pahandle [, repetitions=1] [, when=0] [, waitForStart=0] [, stopTime=inf] [, resume=0])
                                PsychPortAudio('Start', pa.handleHit);
                            elseif pa.correctness(pa.trialNumber) == 0
                                PsychPortAudio('Start', pa.handleMiss);
                            end
                        end
                    end
                    pa.feedbackGiven = 1;
                    pa.waitTime = ds.vbl;
                    
                    
                elseif kb.responseGiven && ds.vbl >  pa.trialOnset + pa.targetMotionDuration + pa.itiLength && ds.vbl <  pa.trialOnset + pa.targetMotionDuration + pa.itiLength + pa.recenterScreenLength
                    
                    glPushMatrix;
                    recenteringDotCoords = [0, 0, -20 1]'; % does this do the right thing??? i'm just updating according to the head rotation info from beginning of frame... important to test explicitly somehow if the dot actually appears in a fixed straight ahead in world coords
                    moglDrawDots3D(ds.w, recenteringDotCoords, 40, [1 0 0 1], [], 2);
                    glPopMatrix;
                    
                    glPushMatrix;
                    invModelView = inv(eye.modelView);
                    moglDrawDots3D(ds.w, inv(modelView)*[pa.fixationVertexPos, 1.0]', pa.fixationDiameter, [1 1 1 1], [], 2); %drawing horizontal fixation dot
                    glPopMatrix;
                    
                elseif kb.responseGiven && ds.vbl >= pa.trialOnset + pa.targetMotionDuration + pa.itiLength + pa.recenterScreenLength % when you''ve responded and it's the last frame of the trial, start new trial
                    
                    [ds, pa, kb, oc] = SetupNewTrial(ds, pa, kb, oc); % setup new trial
                    
                elseif ~kb.responseGiven && ds.vbl >  pa.trialOnset + pa.targetMotionDuration + pa.itiLength && ds.vbl <  pa.trialOnset + pa.targetMotionDuration + pa.itiLength + pa.recenterScreenLength
                    
                    glPushMatrix;
                    recenteringDotCoords = [0, 0, -20 1]'; % does this do the right thing??? i'm just updating according to the head rotation info from beginning of frame... important to test explicitly somehow if the dot actually appears in a fixed straight ahead in world coords
                    if ds.simulated
                        moglDrawDots3D(ds.w, inv(rotationMatrixYawHomo)*recenteringDotCoords, 40, [1 0 0 1], [], 2); % I think this is close to right but still not completely right bc the recentering dot is being updated according to dot world translation or something FIND OUT> ATTN
                    else
                        moglDrawDots3D(ds.w, recenteringDotCoords, 40, [1 0 0 1], [], 2);
                    end
                    glPopMatrix;
                    
                    glPushMatrix;
                    invModelView = inv(eye.modelView);
                    moglDrawDots3D(ds.w, inv(modelView)*[pa.fixationVertexPos, 1.0]', pa.fixationDiameter, [1 1 1 1], [], 2); %drawing horizontal fixation dot
                    glPopMatrix;
                    
                elseif ~kb.responseGiven && ds.vbl >= pa.trialOnset + pa.targetMotionDuration + pa.itiLength + pa.recenterScreenLength % when you've not responded and it's the last frame of the trial, start new trial
                    
                    [ds, pa, kb, oc] = SetupNewTrial(ds, pa, kb, oc); % setup new trial
                    
                end
                
                %% End of experiment logic
                
                % Update view with keyboard responses
                glPopMatrix; % reset back to the origin
                moglDrawDots3D(ds.w, [inf inf inf], 10, [1 1 1 1], [], 2); % 'hack' to fix transparency issue - same one we've used in all 3D pong experiments - definitely works!
                
                % for mask (OLD - 1/f wall)
                %               glBindTexture(GL.TEXTURE_2D,ds.wall_texid);
                %               glCallList(ds.surroundTexture); % 1/f noise texture surround -  comes from CreateTexturesforSDK2.m
                %               glBindTexture(GL_TEXTURE_2D, 0);
                
                
                % Manually disable 3D mode before switching to other eye or to flip:
                Screen('EndOpenGL', ds.w);
                
                % Compute simulation time for this draw cycle:
                ds.tElapsed = (ds.vbl - tStart) * 1;
                
                % Repeat for renderPass of other eye
            end
        end
        if ds.fCount >= 71
           ds.fCount = 1;
%           keyboard;
        end
        
%         psX(k) = posX{1};
%         psY(k) = posY{1}; For plotting optic flow 
%         psZ(k) = posZ{1};
%         k = k + 1;
        for i = 1:shape.plane.numPlanes
            posX{i} = posX{i} + squeeze(vXw(pa.trialNumber, ds.fCount, i, :));
            posY{i} = posY{i} + squeeze(vYw(pa.trialNumber, ds.fCount, i, :));
        end
       
        % Head position tracked in the HMD?
        if ~isempty(ds.hmd)
            if ~bitand(state.tracked, 2) && ds.trackingFlag==1
                % Nope, user out of cameras view frustum. Tell it like it is:
                DrawFormattedText(ds.w, 'Vision based tracking lost\nGet back into the cameras field of view!', 'center', 'center', pa.red);
            end
        end
        
        % Stimulus ready. Show it on the HMD. We don't clear the color buffer here,
        % as this is done in the next iteration via glClear() call anyway:
        Screen('DrawingFinished', ds.w);
        
        Screen('Flip', ds.w,[],[],1);%, [], [], 2);%, ds.vbl + (1-0.5) * ds.ifi);
        % imageArray=Screen(‘GetImage’, windowPtr [,rect] [,bufferName] [,floatprecision=0] [,nrchannels=3])
        ds.vbl = GetSecs();
        ds.fCount = ds.fCount + 1;
    end
    
    while pa.trialNumber == pa.nTrials && ~pa.readyToBeginRun % confirm everything's ready to go
        % when you reach end of last trial
        for renderPass = 0:1 % loop over eyes
            ds.renderPass = renderPass;
            Screen('SelectStereoDrawBuffer',ds.w,ds.renderPass);
            Screen('BeginOpenGL',ds.w);
            
            % Setup camera position and orientation for this eyes view:
            glMatrixMode(GL.PROJECTION)
            glLoadMatrixd(ds.projMatrix{renderPass + 1});
            
            modelView = oc.initialState.modelView{ds.renderPass + 1}; % Use per-eye modelView matrices
            
            % N/A, redundant compared to modelView above. pa.modelViewSaveOutForFixed{ds.renderPass + 1} =  oc.initialState.modelView{ds.renderPass + 1};
            
            glLoadMatrixd(modelView);
            
            Screen('EndOpenGL', ds.w);
            if renderPass == 0
                %Screen('DrawText',ds.w,'TRIAL INFO/ CONDITION INFO',instructionsCenterxL,ds.textCoords(1)-50,[1 1 1]);
                Screen('DrawText',ds.w,'Ready to continue the experiment?',instructionsCenterxL,ds.textCoords(1)-200,[1 1 1]);
                Screen('DrawText',ds.w,'Press SPACE to confirm.',instructionsCenterxL,ds.textCoords(1)+50-200,[1 1 1]);
                
            elseif renderPass == 1
                %Screen('DrawText',ds.w,'TRIAL INFO/ CONDITION INFO',instructionsCenterxR,ds.textCoords(1)-50,[1 1 1]);
                Screen('DrawText',ds.w,'Ready to continue the experiment?',instructionsCenterxR,ds.textCoords(1)-200,[1 1 1]);
                Screen('DrawText',ds.w,'Press SPACE to confirm.',instructionsCenterxR,ds.textCoords(1)+50-200,[1 1 1]);
            end
        end
        
        Screen('DrawingFinished', ds.w);
        ds.vbl = Screen('Flip', ds.w);
        [kb.keyIsDown, kb.secs, kb.keyCode] = KbCheck(-1); % query the keyboard
        if kb.keyIsDown && kb.keyCode(kb.spacebarKey)
            pa.readyToBeginRun=1;
            %keyboard
            
            % save out data before you setup new run
            pa.dataFile = fullfile(pa.Desktop, 'stereoOFdata', [pa.subjectName '-' ds.experimentType '-' pa.date '-' num2str(pa.runNumber) '.mat']);
            save(pa.dataFile, 'pa', 'ds', 'kb','oc');
            
            
            % rerun setup parameters to create new run's trial sequence of
            % rotation velocities
            [ds, pa, kb, oc] = SetupNewRun(ds, pa, kb, oc);
        end
        
        
        
    end
    
    
end

% Calculate average framerate:
fps = ds.fCount / (ds.vbl - tStart), % uncomment to print out at end of run
% Done (or quit out). Save data (in pa.response) and other relevant parameters/variables, close screen and exit
PsychVRHMD('SetAutoClose', ds.hmd, 2);
PsychPortAudio('Close', pa.handleHit);
PsychPortAudio('Close', pa.handleMiss);
Priority(0);
ShowCursor(ds.screenId);
ListenChar(1); % Let's start listening to the keyboard again
% sca;
clear Screen
end
