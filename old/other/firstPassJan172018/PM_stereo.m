%% Phase Motion Stereo
% Charlie Burlingham
% January 17, 2018
% Modification of IK's original PTB script.


%% This program tests heading discrimination for flow fields made of
% plaids. The envelopes are stationary, the patterns move. We're using
% up/down staircases to find observers' su bjective forwards. Start
% values are randomly selected from the range of legal values.   
%Screen('Preference', 'SkipSyncTests', 1);
clear all;close all;clc;                                                    % Start with a blank slate
% HideCursor;                                         
% ListenChar(2); 

%try
    %-----Settings to Modify------------------
    initials                        = '';                                   % Subject initials used to save file
    subjectNumber                   = '';                                   % Number used to save eyelink file
    case_description                = 'translation';                        % Condition type (translation, rotation, tracking) used to save file
    n_trials                        = 25;                                   % Number of trials per staircase
    rotation                        = 0;                                    % deg/s
    fixateMovement                  = 0;                                    % 1 = red fixation dot moves, 0 = stationary red fixation dot
    eyeTracking                     = 0;                                    % 1 = eye tracking on, 0 = eye tracking off
    trainingMode                    = 0;                                    % 1 = gives feedback after each trial (training), 0 = no feedback (actual experiment)
    linearize                       = 0;                                    % Use calibrated LUT (do this when available)
    %data_path                       = '/e/2.3/p2/cburling/Desktop';         % Folder for saving data files
    data_path                       = '/Users/Charlie/Desktop/Optic-Flow-Stereo/firstPassJan172018';         % Folder for saving data files
    
   
    stimulus_duration               = 2;                                    % secs
    translation_speed               = 2;                                    % m/s, apparently 1.9 m/s is a brisk walking speed
    n_staircases                    = 2;                                    % Want to run multiple staircases so that the observer doesn't pick up on string of left/right heading values
    
    % Decides if experiment is simulated rotation or eye tracking and sets variables accordingly
    linear_velocity                 = (rotation * .0174533) * 30;           % m/s, calculated from rotation (angular velocity)
    if fixateMovement == 1
        rotation = 0;
    elseif fixateMovement == 0 
        linear_velocity = 0;
    end
     
    %-----Experiment Settings, Don't Change----
    %-----Movement Directions
    trans_axes                      = [3 1];                                % which way are we driving? 1=X,2=Y,3=Z first number is main direction, second is judged direction
    rot_axis                        = [0 1 0];                              % X (pitch),Y (yaw),Z (roll)
    
    %-----Array Settings
    plane_dist                      = [12.5 25];                            % m
    dot_density                     = [.16 .16];                            % dots/deg^2 on screen, can play around with this. Warren & Hannon was .22
    
    %-----Staircase Settings
    step                            = [4 2 1];
    limits                          = [-20 20];
    
    %-----Element Settings
    gabor_diam                      = 30;                                   % Arcmin
    gabor_sf                        = 2;                                    % c/deg
    gabor_contrast                  = 100;                                  % percent, for the pattern; each element will be half this
    spatial_envelope                = 2;                                    % 0 = disk, 1 = Gabor, 2 = raised cosine
    
    view_window                     = [55 32];                              % X,Y centered around fixation, in dva 
                                                                            % [36 27] for laptop, [55 32] for monitor
    exclude                         = [-20 -3 20 3];                        % Don't place elements where the FOE will be
    
    
    H_ecc_fix                       = 0;                                    % Horizontal fixation ecc (degs, neg is left)
    V_ecc_fix                       = 0;                                    % Vertical fixation ecc (degs, neg is up)
    
    experiment_id                   = 'pattern_heading';                    % Used in group filename
    fast                            = 1;                                    % Automatically trigger trials
    ITI                             = 1;                                    % Intertrial Inverval, Seconds
    background                      = 127;                                  % Grayscale Units
    fixate                          = 1;                                    % Present fixation spot during motion
    
    %-----Rig Settings----------------------
    focalLength                     = .57;                                  % m % CSB: need to change this for the oculus!
    IPD                             = 0.06285% /20;                              % meters. This is the US average btw men and woman
    scale_factor                    = 1.78;                                 % Arcmin/pixel
    frame_rate                      = 120;                                  % Screen frame rate (hz)
    linearize                       = 0;                                    % Use calibrated LUT (do this when available)
    
    %-----Housekeeping----------------------
    % Scale things based on viewing distance, and convert other stuff to
    % the units Psychtoolbox wants...
    nPlanes                         = length(plane_dist);
    tme                             = clock;
    logname = strcat(data_path,initials,'_',experiment_id,'_',case_description,'_trials');
    
    gabor_size                      = gabor_diam/scale_factor;
    stimulus_radius                 = round(gabor_size/2);
    H_ecc_fix                       = H_ecc_fix*60/scale_factor;
    V_ecc_fix                       = V_ecc_fix*60/scale_factor;
    mv_length                       = ceil(stimulus_duration*frame_rate);
    f                               = (gabor_sf*scale_factor/60)*2*pi;
    angle                           = 0;
    a                               = cos(angle)*f;
    b                               = sin(angle)*f;
    amplitude                       = background;
    
    %-----Spatial Envelope------------------
    [x,y]=meshgrid(-stimulus_radius:stimulus_radius,-stimulus_radius:stimulus_radius);
    bps = (stimulus_radius)*2+1;
    circle=((stimulus_radius)^2-(x.^2+y.^2));
    for i=1:bps; for j =1:bps; if circle(i,j) < 0; circle(i,j) = 0; else circle(i,j) = 1; end; end;
    end;
    if spatial_envelope == 1
        circle = (exp(-(((x)/(sqrt(2)*Gaussian_stdev/6)).^2)-((y/(sqrt(2)*Gaussian_stdev/2)).^2)).*circle);
    elseif spatial_envelope == 2
        R = (sqrt(x.^2 + y.^2) + eps).*circle;
        R = R/max(max(R));  
        cos2D = (cos(R*pi)+1)/2;
        circle = (cos2D.*circle);
    end
    circle = circle*255/2;
     
    %-----Open Screens----------------------
    %skip sync tests wasn't here before and these were'nt commented out -CB
    Screen('Preference', 'SkipSyncTests', 1); 
    
    % PsychDebugWindowConfiguration(0,.7);  % just for debugging on a single
    %screen, CSB added 
 
    %oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    %oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    screenNumber=max(Screen('Screens'));

    PsychImaging('PrepareConfiguration'); % CSB added 
    PsychImaging('AddTask', 'General', 'SideBySideCompressedStereo')  % CSB added 
    
    [w, windowRect] = PsychImaging('OpenWindow', screenNumber, background, [], [], [], 102); % CSB added, 102 is side-by-side compressed
    
    %w=Screen('OpenWindow',screenNumber,background,[],[],2); % old 
    
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % Set up alpha-blending for smooth (anti-aliased) drawing of dots:
    screen_rect = Screen('Rect',w);
   
    Screen('FillRect',w, background);
    Screen('Flip', w);
    Screen('FillRect',w, background);
    Screen('TextSize',w,20);
    
    if linearize
        calibration = load('0003_alta_141007.mat'); % CSB: must change this to the new gamma table Hormet and I recorded!
        table = calibration.calib.table;
        psychtoolbox_calib = repmat(table, 1, 3);
        Screen('LoadNormalizedGammaTable',screenNumber,psychtoolbox_calib);
    end
    
    %-----Screen Landmarks------------
    sr_hor = round(screen_rect(3)/2); % Middle of the screen, horizontally, in pixels
    sr_ver = round(screen_rect(4)/2); % Middle of the screen, vertically, in pixels
    fix_hor = sr_hor+H_ecc_fix;     % Horizontal location of fixation cross, in pixels
    fix_ver = sr_ver+V_ecc_fix;     % Vertical location of fixation cross, in pixels
    movie_rect= [0,0,bps,bps];
    
    %----- Set up EyeLink -------------
    if eyeTracking == 1
        mX = view_window(1);
        mY = view_window(2);
        dummymode = 0;
        
        % Provide Eyelink with details about the graphics environment
        % and perform some initializations. The information is returned
        % in a structure that also contains useful defaults
        % and control codes (e.g. tracker state bit and Eyelink key values).
        el=EyelinkInitDefaults(w);
        
        % Initialization of the connection with the Eyelink Gazetracker.
        % exit program if this fails.
        if ~EyelinkInit(dummymode, 1)
            fprintf('Eyelink Init aborted.\n');
            % **here execute any cleanup and shutdown functions**
            Screen('closeall');
            return;
        end
        
        [v vs]=Eyelink('GetTrackerVersion');
        fprintf('Running experiment on a ''%s'' tracker.\n', vs );
        
        fileName=sprintf('%s_RT.edf', subjectNumber); %NB: must be number, not string!!!
        Eyelink('Openfile',fileName);
        
        %configuration settings
        Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, mX-1, mY-1); %mX and mY are max x and y screen coordinates
        Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, mX-1, mY-1);
        
        % make sure that we get gaze data from the Eyelink
        Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
        
        % Calibrate the eye tracker
        EyelinkDoTrackerSetup(el);
        
        % do a final check of calibration using driftcorrection
        EyelinkDoDriftCorrection(el);
        
        % start recording eye position
        Eyelink('StartRecording');
        
        % record a few samples before we actually start displaying
        WaitSecs(0.1);
        
        % mark zero-plot time in data file
        Eyelink('Message', 'SYNCTIME');
        stopkey=KbName('space');
        eye_used = -1;
    end
    
    %-----Make the grating in all phases-----
    gabor = zeros(1,360);
    for i=1:360
        grating = round(((sin(a*x+b*y+i*pi/180)*amplitude)+background));
        gabor(i) = Screen('MakeTexture',w,cat(3,grating,circle));
    end
    
    %-----Set Up Conditions-----------
    n_conditions = 0;
    for j=1:length(rotation)
        for k=1:n_staircases
            n_conditions = n_conditions+1;
            
            cond(n_conditions).rot_index = j;
            
            cond(n_conditions).rotate = [0 0 0];
            cond(n_conditions).rotate = rot_axis*rotation(j)*pi/180;
            
            cond(n_conditions).contrast = gabor_contrast/100;
            cond(n_conditions).count = 0;
            
            cond(n_conditions).last_resp = 0;
            cond(n_conditions).n_flips = 0;
            
            cond(n_conditions).heading = limits(1)+ceil(rand*(limits(2)-limits(1)));
            cond(n_conditions).step = step(1);
            
            cond(n_conditions).resp_history = zeros(1,n_trials);
            cond(n_conditions).heading_history = zeros(1,n_trials);
        end
    end
    results = zeros(n_conditions,n_trials);
    
    %-----Randomize Trials------------
    total_trials = n_trials*n_conditions;
    perm = randperm(total_trials);
    perm = mod(perm,n_conditions)+1;
    duration_check = zeros(size(perm));
    
    Screen('SelectStereoDrawBuffer', w, 0); % left eye
    Screen('DrawText',w,'Heading discrimination: Pattern',500,240,250);
    Screen('DrawText',w,'Use Left/Right Arrows to respond',500,270,250);
    Screen('DrawText',w,'press SPACE BAR to start',500,300,250);
    
    Screen('SelectStereoDrawBuffer', w, 1); % right eye
    Screen('DrawText',w,'Heading discrimination: Pattern',500,240,250);
    Screen('DrawText',w,'Use Left/Right Arrows to respond',500,270,250);
    Screen('DrawText',w,'press SPACE BAR to start',500,300,250);
    
    Screen('Flip',w);
    
    FlushEvents('keyDown');
    validKey = 0;
    while ~validKey
        [secs, keyCode, deltaSecs] = KbWait(-1);
        if keyCode(KbName('space'))
            validKey = 1;
        end
    end
    
    for eye = 0:1 % eyes
        Screen('SelectStereoDrawBuffer', w, eye);
        Screen('FillRect', w, background);
    end
    
    Screen('Flip', w);
    tic;
    
    %-----Main experimental loop-----------------
    for trial=1:total_trials
        clear angle; clear pos; % CSB added. do this b/c these might have been used earlier. 
        
        if eyeTracking == 1
            % Divides eye tracking data by trial
            eyemsg = sprintf('trial %i', trials);
            Eyelink('Message', eyemsg);
        end
        
        aa = GetSecs;
        % Draw 
        % Set up the world
        translate = [0 0 0];
        translate(trans_axes(1)) = translation_speed*cosd(cond(perm(trial)).heading);
        translate(trans_axes(2)) = translation_speed*sind(cond(perm(trial)).heading);
        

        [dots2D_L{1}, dots2D_R{1}, dots3D_L{1}, dots3D_R{1}] = make_dot_plane_stereo2(dot_density(1), plane_dist(1), view_window, exclude, focalLength, IPD); % plane 1, for left and right eye
        [dots2D_L{2}, dots2D_R{2}, dots3D_L{2}, dots3D_R{2}] = make_dot_plane_stereo2(dot_density(2), plane_dist(2), view_window, exclude, focalLength, IPD); % plane 2, for left and right eye
        
        dots3D = {dots3D_L{1} dots3D_R{1}; dots3D_L{2} dots3D_R{2}};
        dots2D = {dots2D_L{1} dots2D_R{1}; dots2D_L{2} dots2D_R{2}};
        
        for e = 1:2 % left and right eyes
            clear tempo; % clear everything that will fuck up computation of the other eye's movie
            tempo.angle_all = [];
            tempo.speed_all = [];
            tempo.pos_all = [];
            
            for i=1:nPlanes
                tempo.velocity_field = calculate_plane_flow_stereo(dots3D{i,e}, dots2D{i,e}, translate, cond(perm(trial)).rotate, focalLength);
                
                tempo.nGabors(i) = size(dots2D{i,e},2);
                % now convert dot centers from m back to deg JUST in order to center rect on point:
                tempo.dots2D_deg{i,e} = atand(dots2D{i,e}/focalLength);
                tempo.posA = CenterRectOnPoint(movie_rect, tempo.dots2D_deg{i,e}(1,:)'*60/scale_factor+sr_hor, tempo.dots2D_deg{i,e}(2,:)'*60/scale_factor+sr_ver)'; %big Y numbers are lower on the screen in pixels so we have to flip
                tempo.pos = [tempo.posA tempo.posA];
                
                tempo.speed_ver = tempo.velocity_field(2,:);
                tempo.speed_hor = tempo.velocity_field(1,:);
                
                tempo.angle_ver = 270*(tempo.speed_ver>=0)+90*(tempo.speed_ver<0); % in deg
                tempo.angle_hor = 180*(tempo.speed_hor>0); % in deg
                tempo.speed_ver = abs(tempo.velocity_field(2,:));
                tempo.speed_hor = abs(tempo.velocity_field(1,:));
                
                tempo.speed = zeros(1,length(tempo.speed_ver)*2);
                tempo.angle = tempo.speed;
                tempo.speed = [tempo.speed_hor tempo.speed_ver];
                tempo.angle = [tempo.angle_hor tempo.angle_ver];
                tempo.angle_all = [tempo.angle_all tempo.angle];
                tempo.speed_all = [tempo.speed_all tempo.speed];
                tempo.pos_all = [tempo.pos_all tempo.pos];
            end
            
            tempo.angle = tempo.angle_all;
            angle{e} = tempo.angle; % store for later for displaying dots, CSB added
            tempo.speed = tempo.speed_all;
            tempo.pos = tempo.pos_all;
            pos{e} = tempo.pos; % store for later for displaying dots, CSB added
            tempo.phase_step = round(tempo.speed*gabor_sf*360/frame_rate)';
            
            cond(perm(trial)).count = cond(perm(trial)).count+1;
            tempo.phases = repmat(ceil(rand(size(tempo.speed,2),1)*360),1,mv_length)+repmat(tempo.phase_step,1,mv_length).*repmat(0:(mv_length-1),size(tempo.speed,2),1);
            tempo.phases = rem(tempo.phases,360);
            phasesCell{e} = tempo.phases + (tempo.phases==0)*360;  %CSB added
        end

        %Finish the ITI
        WaitSecs(ITI-(GetSecs-aa));
        
        % Draw the red fixation cross
        for eye = 0:1 % eyes
            Screen('SelectStereoDrawBuffer', w, eye);
            Screen('DrawDots', w, [fix_hor; fix_ver], 4, [250 0 0], [], 2);
        end 
        Screen('Flip',w); 
         
        FlushEvents('keyDown');
        priorityLevel=MaxPriority(w);
        Priority(priorityLevel);
        
        % Play the movie
        StimulusOnsetTime = zeros(1,mv_length);
        aa = GetSecs;
        tag = clock;
        for frame = 1:mv_length
            for eye = 1:2 % eyes
                Screen('SelectStereoDrawBuffer', w, eye-1);  % CSB added, needs to take 0 or 1 for left and right eyes
                if fixateMovement
                    posFixate = frame * linear_velocity;
                    Screen('DrawTextures', w, gabor(phasesCell{eye}(:,frame)),movie_rect,pos{eye},angle{eye});
                    if fixate
                        Screen('DrawDots', w, [fix_hor + posFixate; fix_ver], 4, [250 0 0], [], 2);
                    end
                    
                else
                    Screen('DrawTextures', w, gabor(phasesCell{eye}(:,frame)),movie_rect,pos{eye},angle{eye});
                    if fixate
                        Screen('DrawDots', w, [fix_hor; fix_ver], 4, [250 0 0], [], 2);
                    end
                end
            end
            [VBLTimestamp, StimulusOnsetTime(frame), FlipTimestamp, Missed, Beampos] = Screen('Flip',w);
        end
        
        duration_check(trial) = GetSecs-aa;
        
        for eye = 0:1 % eyes
            Screen('SelectStereoDrawBuffer', w, eye);
            Screen('FillRect',w, background);
        end
        
        if fixate
            for eye = 0:1 % eyes
                Screen('SelectStereoDrawBuffer', w, eye);
                Screen('DrawDots', w, [fix_hor; fix_ver], 4, [250 0 0], [], 2); % CSB: do I need to offset fixation cross to do stereo???
            end
        end
        Screen('Flip',w);
        
        % Get the response
        validKey = 0;
        while ~validKey
            [secs, keyCode, deltaSecs] = KbWait(-1);
            if keyCode(KbName('ESCAPE'))
                Screen('CloseAll');
                Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
                Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
                ListenChar(1);
                ShowCursor;
                break
            elseif keyCode(KbName('LeftArrow'))
                validKey = 1;
                resp = 0;
            elseif keyCode(KbName('RightArrow'))
                validKey = 1;
                resp = 1;
                results(perm(trial),cond(perm(trial)).count) = 1;
            else
                Beeper('low');
            end
            
            % Provides feedback during training phase
            if trainingMode
                if cond(perm(trial)).heading < 0 && resp == 0 
                    Screen('DrawText',w,'Correct',625,375,250);
                    Screen('Flip',w)
                    WaitSecs(ITI)
                elseif cond(perm(trial)).heading < 0 && resp == 1 
                    Screen('DrawText',w,'Incorrect',625,375,250);
                    Screen('Flip',w)
                    WaitSecs(ITI)
                elseif cond(perm(trial)).heading > 0 && resp == 1
                    Screen('DrawText',w,'Correct',625,375,250);
                    Screen('Flip',w)
                    WaitSecs(ITI)
                elseif cond(perm(trial)).heading > 0 && resp == 0
                    Screen('DrawText',w,'Incorrect',625,375,250);   
                    Screen('Flip',w)
                    WaitSecs(ITI)
                elseif cond(perm(trial)).heading == 0
                    Screen('DrawText',w,'Heading was 0',625,375 ,250);
                    Screen('Flip',w)
                    WaitSecs(ITI) 
                end
            end
        end
        Priority(0); 
        
        % Log the trial
        fid = fopen(logname,'a');
        fprintf(fid,'%s, %s, %d, %d, %d, %d, %d, %f, %d, %f, %f, %d, %f, %f, %d, %f, %f\n',initials,experiment_id, tag(1), tag(2), tag(3), tag(4), tag(5), tag(6), trial, translation_speed, rotation, fixateMovement, linear_velocity, cond(perm(trial)).heading, resp, stimulus_duration, duration_check(trial));        
        fclose(fid);
        
        % Tell the staircase what happened
        cond(perm(trial)).resp_history(cond(perm(trial)).count) = resp;
        cond(perm(trial)).heading_history(cond(perm(trial)).count) = cond(perm(trial)).heading;
        % Was this a flip? If so, decrease the step size
        if resp ~= cond(perm(trial)).last_resp && cond(perm(trial)).count > 1
            cond(perm(trial)).n_flips = cond(perm(trial)).n_flips + 1;
            if cond(perm(trial)).n_flips < length(step)
                cond(perm(trial)).step = step(cond(perm(trial)).n_flips+1);
            else
                cond(perm(trial)).step = step(length(step));
            end
        end
        
        % Set the new heading
        if resp
            cond(perm(trial)).heading = cond(perm(trial)).heading-cond(perm(trial)).step;
            if cond(perm(trial)).heading < limits(1)
                cond(perm(trial)).heading = limits(1);
            end
        else
            cond(perm(trial)).heading = cond(perm(trial)).heading+cond(perm(trial)).step;
            if cond(perm(trial)).heading > limits(2)
                cond(perm(trial)).heading = limits(2);
            end
        end
        
        % Set the new last value
        cond(perm(trial)).last_resp = resp;
        
    end
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    ListenChar(1);
    clc;
    
    % Print out results here.
    for i=1:n_conditions
       if rotation(cond(i).rot_index) > 0
           plot(cond(i).heading_history,'bo-');
           hold on;
       elseif rotation(cond(i).rot_index) < 0
           plot(cond(i).heading_history,'ro-');
           hold on;
       else
           plot(cond(i).heading_history,'ko-');
           hold on;
       end
    end
    filename = strcat(data_path,initials,'_',experiment_id,'_',int2str(tme(1)),'_',int2str(tme(2)),'_',int2str(tme(3)),'_',int2str(tme(4)),'_',int2str(tme(5)));
    fprintf('\n');
    time = toc/60;
    fprintf('Elapsed time  (minutes)  =  %4.1f\n',time);
    clear x y R grating cos2D circle
    save(filename);
    
    %----- Close EyeLink ---------
    if eyeTracking == 1
        
        % wait a while to record a few more samples
        WaitSecs(0.1);
        
        % finish up: stop recording eye-movements,
        % close graphics window, close data file and shut down tracker
        Eyelink('StopRecording');
        Eyelink('closefile');
        
        % download data file
        try
            fprintf('Receiving data file ''%s''\n', fileName);
            status=Eyelink('ReceiveFile',fileName);
            if status > 0
                fprintf('ReceiveFile status %d\n', status);
            end
            if 2==exist(fileName, 'file')
                fprintf('Data file ''%s'' can be found in ''%s''\n', fileName, pwd );
            end
        catch rdf
            fprintf('Problem receiving data file ''%s''\n', fileName );
            rdf;
        end
        
        cleanup;
        % Shutdown Eyelink:
        Eyelink('Shutdown');
        % Close window:
        sca;
    end
    %{
catch ME
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    ListenChar(1);
    ShowCursor;
    Screen('CloseAll');
    %these weren;t commented out - CB
    %Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    %Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    Priority(0);
end %try..catch
    %} 