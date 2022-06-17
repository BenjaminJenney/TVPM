function RunTVPMExperiment(subjectInitials,condition,monocularFlag,trainingFlag,mask, streamToHMDFlag)
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

% Add paths
addpath(genpath(pwd));

%% Init input variables

pa.subjectName     = subjectInitials; % initials of participant
ds.experimentType  = condition; % 'real', 'simulated', 'tvpmsd', 'tvpmmask'
ds.monocularFlag   = monocularFlag; % 1: monocular viewing, 2: binocular viewing
ds.oculusConnected = streamToHMDFlag; %1: show image in hmd, 0: show image in monitor
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
    PsychDebugWindowConfiguration([],.7); % display on monitor at half opacity in addition to in the HMD
end

%% Setup Psychtoolbox for OpenGL 3D rendering support
% and initialize the mogl OpenGL for Matlab/Octave wrapper:
global GL; % GL data structure needed for all OpenGL programs

InitializeMatlabOpenGL(1);
PsychDefaultSetup(2); % the input of 2 means: execute the AssertOpenGL command, execute KbName('UnifyKeyNames') routine, AND unifies the color mode and switches from 0-255 to 0.0-1.0 - color part only impacts the current function or script, not ones that are called

%% Initialize screen, experimental parameters, and keyboard
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
type = GL.TEXTURE_2D;
Screen('EndOpenGL', ds.w);
while ~readyToBegin % confirm everything's ready to go
    
    % Camera position when using head tracking + HMD: (according to SuperShapeDemo.m)
    globalPos = [0, 0, 0]; % x,y,z  % in meters - this is the starting position of head... should be at origin
    heading = 0; % yaw.
    ds.globalHeadPose = PsychGetPositionYawMatrix(globalPos, heading); % just an identity matrix to start with % initialize observer's start position to the default camera position specified above
    
    if isempty(ds.hmd) % Oculus not connected
        load DefaultHMDParameters.mat;
        % oc.defaultState = defaultState;
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
            
            modelView = oc.defaultState.modelView{ds.renderPass + 1}; % SET BACK to initialState.modelView for normal camera Use per-eye modelView matrices
            
            % N/A, redundant compared to modelView above. pa.modelViewSaveOutForFixed{ds.renderPass + 1} =  oc.initialState.modelView{ds.renderPass + 1};
            glMatrixMode(GL.MODELVIEW)
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
    Screen('DrawingFinished', ds.w);
    ds.vbl = Screen('Flip', ds.w);
    
    [kb.keyIsDown, kb.secs, kb.keyCode] = KbCheck(-1); % query the keyboard
    if kb.keyIsDown && kb.keyCode(kb.spacebarKey)
        readyToBegin=1;
    end
end


%% Participant is ready, so let's start first run of trials

ds.tElapsed = 0;
ds.fCount = 0;

[ds, pa, kb] = SetupNewTrial(ds, pa, kb);
ds.vbl = pa.trialOnset;
tStart = ds.vbl;
pa.experimentOnset = ds.vbl;

pa.block = 0;
kb.nextTrialKey = 0;


while (pa.runNumber <= pa.numRuns) && ~kb.keyCode(kb.escapeKey)
    if ds.tvpmsd
       [vXw, vYw, shape] = Preprocess(ds, pa, shape, GL);% m/f in world coords
    end
    
    while (pa.trialNumber < pa.nTrials) && ~kb.keyCode(kb.escapeKey) % wait until all of the trials have been completed or the escape key is pressed to quit out

        if (ds.tvpmsd) && ds.vbl <  pa.trialOnset + pa.targetMotionDuration  % if still presenting stim
            
%             startTime = GetSecs();
%             endTime = ds.vbl;
%             
%             if (currentTime - lastTime >= 1.0)
%                 disp('ms/frame', 1000.0/ds.fCount);
%                 ds.fCount = 0;
%                 lastTime = 
            ds.fCount = ds.fCount + 1;
            if ds.fCount == 1 % initialize the plaid positions with the initial positions generated by DiskPos.m which is called in SetupShapes.m
                clear posX; clear posY;
                for i = 1:shape.plane.numPlanes
                    posX{i} = shape.disk.X_m{i};
                    posY{i} = shape.disk.Y_m{i};
                    posZ{i} = ones(size(posX{i})) .* shape.plane.depths_m(i);
                end
            end
            if ds.tvpmsd
%                 % update the current plaid positions and modulate the
%                 % plaids in the virtual world
                [posX, posY] = updatePositionsWithOpticFlowAndModulate(shape.plane.numPlanes,... % flow is different for each plane since they are at different depths
                                                                       pa.trialNumber,... % current trial
                                                                       ds.fCount,...   % frame count
                                                                       shape.cyc_m,... % for modulation: cycle of sinwave -- plaid is two sinewaves interleaved
                                                                       posX, posY,...  % current plaid position.
                                                                       vXw, vYw,...    % how much the plaid is moved in VR world for current frame: displacement generated from optic flow (Preprocess.m).
                                                                       shape.disk.X_m, shape.disk.Y_m); % initial plaid poisitions.
            end
        end
        
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
                    
                    
                end
                
                Screen('SelectStereoDrawbuffer', ds.w, eye.eyeIndex); % Select 'eyeIndex' to render (left- or right-eye):
                
                modelView = oc.defaultState.modelView{renderPass+1};%eye.modelView;
                
                Screen('BeginOpenGL', ds.w); % Manually reenable 3D mode in preparation of eye draw cycle
                
                % Setup camera position and orientation for this eyes view:
                glMatrixMode(GL.PROJECTION)
                glLoadMatrixd(ds.projMatrix{renderPass + 1});
                
                glMatrixMode(GL.MODELVIEW);
                
                glLoadMatrixd(modelView); % updates openGL camera
                
                glClearColor(0.5, 0.5, 0.5, 0.5); % white background
                glClear(); % clear the buffers - must be done for every frame
                glColor3f(1,1,1);
                
                
                glPushMatrix;
                %% Experiment Logic
                if ds.vbl <  pa.trialOnset + pa.targetMotionDuration % if current time < present until time: draw target, 1 s target motion
                    
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
                    
                    
                    xDotDisplacement = pa.transSpeed.*sind(-1.*pa.heading(pa.trialNumber)).*elapsedTime; % CSB: should be in meters. check. % ATTN: take negative of heading angle, so it displays correctly. otherwise displayed heading is inverse of ground truth
                    % pa.rotVelocityRadPursuit(pa.trialNumber) = deg2radtrialNumber)).*elapsedTime; % CSB: should be in meters. check. % ATTN: take negative of heading angle, so it displays correctly. otherwise displayed heading is inverse of ground truth
                    zDotDisplacement = pa.transSpeed.*cosd(-1.*pa.heading(pa.trialNumber)).*elapsedTime;  % JT: make sure this trig is right.
                    
                    xDotDisplacementTraining = pa.transSpeedTraining.*sind(-1.*pa.heading(pa.trialNumber)).*elapsedTime;
                    zDotDisplacementTraining = pa.transSpeedTraining.*cosd(-1.*pa.heading(pa.trialNumber)).*elapsedTime;
                    
                    
                    
                    if ds.simulated
                        
                        %gluLookAt(-xDotDisplacement,0,-zDotDisplacement,0,0,shape.plane.depths_m(2),0,1,0);
                        for i = 1:shape.plane.numPlanes
                            
                            glPushMatrix;
                            glTranslatef(shape.plane.offsets_m(i), 0.0, 0.0);
                            
                            moglDrawDots3D(ds.w, [shape.disk.X_m{i} shape.disk.Y_m{i} shape.disk.Z_m{i}]', 10, [1 1 1 1], [], 2);
                            glPopMatrix;
                        end
                     
                        moglDrawDots3D(ds.w, [0 0 shape.plane.depths_m(2) 1]', 3*2, pa.red, [], 2); %drawing fixation dot in environ % CSB june 21 2021 uncomment
                    elseif ds.real
                        glPushMatrix;
                        glTranslatef(xDotDisplacementTraining,0,zDotDisplacementTraining); % shift the dot world to its position along its trajectory for this frame
                        moglDrawDots3D(ds.w, pa.sphere, 3, [1 1 1 1], [], 2); %drawing sphere
                        glPopMatrix;
                        
                        glPushMatrix;
                        glTranslatef(xDotDisplacementTraining,0,zDotDisplacementTraining);
                        moglDrawDots3D(ds.w, [0 0 -pa.cubeWidth/4 1]', 3*2, pa.red, [], 2); %drawing fixation dot in environ % CSB june 21 2021 uncomment
                        glPopMatrix;
                        if mask == 1, drawTheHoleyBag(ds, 1), end % draws the ptb mask. Params: drawTheHoleyBag(window,[openglEnbled = 1]);
                    elseif ds.tvpmsd
                        
                        drawFixationDot(GL, ds.w, pa.fixationDiameter, shape.plane.depths_m(2));
                        
                        drawPlaidsForTVPMSD(shape.disk.listIds,...
                                            shape.plane.numPlanes,...
                                            shape.disk.numDisksPerPlane,... 
                                            posX, posY); % draws the plaids bounded by the 3 planes (numPlanes = 3) these are the same 3 planes drawn in TVPMCD
                                          
                        drawMasksForTVPMSD(shape.mask.widths_m, shape.mask.listIds); % Draw the three masks for TVPMSD. The masks are projected .001 meters in front of the plaids (and as such are marginally different in size from the three planes drawn in TVPMCD). So there depths are shape.plane.depths_m - .001. See SetupShapes for reference. 
                        
                    elseif ds.tvpmcd
                        
                        % Draw fixation dot
                        drawFixationDot(GL, ds.w, pa.fixationDiameter, shape.plane.depths_m(2));
                        
                        % Draw the 3 plaid planes
                        drawPlaidPlanesForTVPMCD(shape.plane.widths_m, shape.plane.listIds, xDotDisplacement, zDotDisplacement, shape.plane.depths_m(2));
                        
                        % Draw the full screen mask
                        if renderPass == 0
                            %drawFullScreenMaskForTVPMCD(ds.w, shape.mask.fullWindowMaskLeftEye)
                        end
                        if renderPass == 1
                            %drawFullScreenMaskForTVPMCD(ds.w, shape.mask.fullWindowMaskRightEye)
                        end
                        
                    end
                    
                    
                elseif ~kb.responseGiven && ds.vbl >  pa.trialOnset + pa.targetMotionDuration && ds.vbl <=  pa.trialOnset + pa.targetMotionDuration + pa.itiLength  % show paddle and allow observers to adjust its position - no time constraints - they press the space bar to lock in their response and start a new trial
                    
                    if ds.real
                        glPushMatrix;
                        glTranslatef(xDotDisplacementTraining,0,zDotDisplacementTraining); % shift the target to its position along its trajectory for this frame
                        moglDrawDots3D(ds.w, pa.sphere, 3, [1 1 1 1], [], 2); %drawing sphere
                        moglDrawDots3D(ds.w, [0 0 -pa.cubeWidth/4], 3*2, pa.red, [], 2); %drawing fixation dot in environ % CSB june 21 2021 uncomment
                        glPopMatrix;
                    elseif ds.simulated
                        glPushMatrix;
                        gluLookAt(-xDotDisplacementTraining,0,-zDotDisplacementTraining,0,0,-15,0,1,0)
                        moglDrawDots3D(ds.w, pa.sphere, 3, [1 1 1 1], [], 2); %drawing sphere
                        glPopMatrix;
                        moglDrawDots3D(ds.w, [0 0 -pa.cubeWidth/4 1]', 3*2, pa.red, [], 2); %drawing fixation dot in environ % CSB june 21 2021 uncomment
                        
                        
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
                    
                    
                    if ds.stimulus || ds.real
                        glPushMatrix;
                        recenteringDotCoords = [0, 0, -20 1]'; % does this do the right thing??? i'm just updating according to the head rotation info from beginning of frame... important to test explicitly somehow if the dot actually appears in a fixed straight ahead in world coords
                        moglDrawDots3D(ds.w, recenteringDotCoords, 40, [1 0 0 1], [], 2);
                        glPopMatrix;
                        
                        glPushMatrix;
                        invModelView = inv(eye.modelView);
                        moglDrawDots3D(ds.w, inv(modelView)*[pa.fixationVertexPos, 1.0]', pa.fixationDiameter, [1 1 1 1], [], 2); %drawing horizontal fixation dot
                        glPopMatrix;
                    else
                        %Drawing fixation dot
                        glPushMatrix;
                        curModelViewNoTranslation1 = glGetFloatv(GL.MODELVIEW_MATRIX);
                        curModelViewNoTranslation2 = [curModelViewNoTranslation1(1:4),curModelViewNoTranslation1(5:8),curModelViewNoTranslation1(9:12),[0 0 0 1]']; % CSB, 4/12/2022- had to use the GL modelview without translation to make the fixation pt at right depth. it was broken b4 using eye.modelView. figure out why
                        moglDrawDots3D(ds.w, inv(curModelViewNoTranslation2)*[[0 0 shape.plane.depths_m(2) + .001], 1.0]', pa.fixationDiameter, [1 1 1 1], [], 2); %drawing horizontal fixation dot
                        glPopMatrix;
                        
                        glPushMatrix;
                        recenteringDotCoords = [0, 0, shape.plane.depths_m(2), 1]'; % does this do the right thing??? i'm just updating according to the head rotation info from beginning of frame... important to test explicitly somehow if the dot actually appears in a fixed straight ahead in world coords
                        moglDrawDots3D(ds.w, recenteringDotCoords, 40, [1 0 0 1], [], 2);
                        glPopMatrix;
                        
                        
                    end
                elseif kb.responseGiven && ds.vbl >= pa.trialOnset + pa.targetMotionDuration + pa.itiLength + pa.recenterScreenLength % when you''ve responded and it's the last frame of the trial, start new trial
                    
                    [ds, pa, kb, oc] = SetupNewTrial(ds, pa, kb, oc); % setup new trial
                    
                elseif ~kb.responseGiven && ds.vbl >  pa.trialOnset + pa.targetMotionDuration + pa.itiLength && ds.vbl <  pa.trialOnset + pa.targetMotionDuration + pa.itiLength + pa.recenterScreenLength
                    
                    if ds.simulated || ds.real
                        
                        recenteringDotCoords = [0, 0, -20, 1]'; % does this do the right thing??? i'm just updating according to the head rotation info from beginning of frame... important to test explicitly somehow if the dot actually appears in a fixed straight ahead in world coords
                        glPushMatrix;
                        if ds.simulated
                            moglDrawDots3D(ds.w, recenteringDotCoords, 40, [1 0 0 1], [], 2);
                        end
                        glPopMatrix;
                        
                        glPushMatrix;
                        invModelView = inv(eye.modelView);
                        moglDrawDots3D(ds.w, inv(modelView)*[pa.fixationVertexPos, 1.0]', pa.fixationDiameter, [1 1 1 1], [], 2); %drawing horizontal fixation dot
                        glPopMatrix;
                    else
                        %Drawing fixation dot
                        glPushMatrix;
                        curModelViewNoTranslation1 = glGetFloatv(GL.MODELVIEW_MATRIX);
                        curModelViewNoTranslation2 = [curModelViewNoTranslation1(1:4),curModelViewNoTranslation1(5:8),curModelViewNoTranslation1(9:12),[0 0 0 1]']; % CSB, 4/12/2022- had to use the GL modelview without translation to make the fixation pt at right depth. it was broken b4 using eye.modelView. figure out why
                        moglDrawDots3D(ds.w, inv(curModelViewNoTranslation2)*[[0 0 shape.plane.depths_m(2)+.001], 1.0]', pa.fixationDiameter, [1 1 1 1], [], 2); %drawing horizontal fixation dot
                        glPopMatrix;
                        
                        glPushMatrix;
                        recenteringDotCoords = [0, 0, shape.plane.depths_m(2), 1]';
                        moglDrawDots3D(ds.w, recenteringDotCoords, 40, [1 0 0 1], [], 2);
                        glPopMatrix;
                        
                        %                         glPushMatrix;
                        %                         invModelView = inv(eye.modelView);
                        %                         moglDrawDots3D(ds.w, inv(modelView)*[[0 0 shape.plane.depths_m(2)], 1.0]', pa.fixationDiameter, [1 1 1 1], [], 2); %drawing horizontal fixation dot
                        %                         glPopMatrix;
                    end
                    
                elseif ~kb.responseGiven && ds.vbl >= pa.trialOnset + pa.targetMotionDuration + pa.itiLength + pa.recenterScreenLength % when you've not responded and it's the last frame of the trial, start new trial
                    
                    [ds, pa, kb, oc] = SetupNewTrial(ds, pa, kb, oc); % setup new trial
                    
                end
                
                %% End of experiment logic
                
                % Update view with keyboard responses
                glPopMatrix; % reset back to the origin
                moglDrawDots3D(ds.w, [inf inf inf], 10, [1 1 1 1], [], 2); % 'hack' to fix transparency issue - same one we've used in all 3D pong experiments - definitely works!
                
                % Manually disable 3D mode before switching to other eye or to flip:
                Screen('EndOpenGL', ds.w);
                
                % Compute simulation time for this draw cycle:
                ds.tElapsed = (ds.vbl - tStart) * 1;
                
                % Repeat for renderPass of other eye
            end
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
        
        Screen('Flip', ds.w,[],[],1);
        ds.vbl = GetSecs();
        
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
            glMatrixMode(GL.MODELVIEW) % BJ: 3/17 added glMatMode call, if you load modelView it makes sense that we should be in that mode.
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
            
            % save out data before you setup new run
            pa.dataFile = fullfile(pa.Desktop, 'stereoOFdata', [pa.subjectName '-' ds.experimentType '-' pa.date '-' num2str(pa.runNumber) '.mat']);
            save(pa.dataFile, 'pa', 'ds', 'kb','oc');
            
            
            % rerun setup parameters to create new run's trial sequence of
            % rotation velocities
            [ds, pa, kb, oc] = SetupNewRun(ds, pa, kb, oc);
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
