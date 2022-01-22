function tvpmExperiment(subjectInitials, condition, monocularFlag, trainingFlag)
    %%Clean up workspace
    close all;
    sca;
    
    %% Init input variables
    pa.subjectName    = subjectInitials; % initials of participant
    ds.experimentType = condition;       % 'real', 'simulated', 'stabalized' %BJ to CS: Do we care about stabalized?
    ds.monocularFlag  = monocularFlag;   % 1: monocular viewing, 2: binocular viewing
    pa.TRAINING       = trainingFlag;    % 0: full expirement (all rotation velocities and headings) 1: training runs (will run two runs without rotation)
    
    global DEBUG_FLAG MONOCULAR AUDIO
    DEBUG_FLAG = 1; % 1: if you're debugging and want to print info to command line, and to show display on screen at half opacity ,  0: turn off those things
    MONOCULAR = ds.monocularFlag ;

    if pa.TRAINING == 1 AUDIO = 1; else AUDIO = 0; end
    AUDIO = 1; % to override and have feedback even when you're not training

    if DEBUG_FLAG == 1
        PsychDebugWindowConfiguration([],.4); % display on monitor at half opacity in addition to in the HMD
    end
    
    %% Setup Psychtoolbox for OpenGL 3D rendering support
    % and initialize the mogl OpenGL for Matlab/Octave wrapper:
    global GL; % GL struct for accessing OpenGL state options in multiple scripts (GL_TEXTURE2D --> GL.TEXTURE2D)
    InitializeMatlabOpenGL(1);
    PsychDefaultSetup(2);
    
    %% Initialize screen, experimental parameters, opengl shapes, and keyboard
    [ds,oc] = SetupDisplay(ds);       % set up the display, based on the DK2
    [ds,pa] = SetupParameters(ds,pa); % set up the experimental parameters for this session
    shape   = SetupShapes(ds, pa);    % set up plane and disk properties.
    kb      = SetupKeyboard(ds,DEBUG_FLAG); % get the keyboard info for the participant's responses
    
    %% Start the experiment - opening screen, getting the participant set

    pa.runNumber = 1; % initialize the run counter
    pa.readyToBeginRun = 0; % initialize the ready-to-begin-a-new-run variable

    readyToBegin = 0;
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

        Screen('DrawingFinished', ds.w);
        ds.vbl = Screen('Flip', ds.w);

        [kb.keyIsDown, kb.secs, kb.keyCode] = KbCheck(-1); % query the keyboard
        if kb.keyIsDown && kb.keyCode(kb.spacebarKey)
            readyToBegin=1;
        end
    end
    
end