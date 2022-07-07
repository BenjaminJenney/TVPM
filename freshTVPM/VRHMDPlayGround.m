function VRHMDPlayGround(streamToHMDFlag, projMatrixIdentityFlag, modelView)
 %{ 
                        BJ: Observances:
                        . 
                          stim: drew one dot at [x=0;y=0;z=-1]
                          params: projMatrix = identity (glLoadIdentity);
                                  modelView  = identity;
                                  modelView  = defaultState.modelViewDataLeft,...Right Same as old/JT_exp
                          perception: see two dots presented in each eye view. Closing one eye you can see the dot appears at something close to center of eye view. When using identity for Modelview headtracking does not take effect (as is expected) as such dot is glued to the headset. 
                          On monitor: When modelView is identity dot appears to be dead center of stereo display.
                                      When modelView is defaultState.modelView... dot does not appear at all.
 
                        . 
                            stim: drew one dot at [x=0;y=0;z=-1]
                            params: projMatrix = PsychVRHMD('GetStaticRenderParameters', hmd) same as JT_exp
                                  modelView  = identity;
                                  modelView  = defaultState.modelViewDataLeft,...Right Same as old/JT_exp
                            perception: when MV = identity, dot is fused
                                        at center of fov (no headtracking)
                                        When MV = defaultState... dot is
                                        seemingly fused at center of fov (hard to tell exactly with head tracking on)
                            monitor: When MV = identity dot appears vertically shifted up *in seemingly equal magnitude to our vertical shift* of each stereo display.
                                     When MV = defaultState... dot
                                               appears nearer center of
                                               screen (no more vertical
                                               shift?), still based on my
                                               measurements does not appear
                                               at vertical center of stereo
                                               display, shifted up a little
                                               this is also the case in
                                               JT_exp
                        .
                             stim: drew plane with width = Vangle(hfov=80, vfov=90)[x=0;y=0;z=-1]
                            params: projMatrix = PsychVRHMD('GetStaticRenderParameters', hmd) same as JT_exp
                                  modelView  = identity;
                                  modelView  = defaultState.modelViewDataLeft,...Right Same as old/JT_exp
                            perception: When MV = identity appears to take all of the vertical fov and nearly all of the
                                        horizontal fov, glued to headset.
                                        appears at dead center.
                                        When MV = defaultState... seams to
                                        be dead center as well and takes up all of vertival direction and most of horizontal though
                                        headtracking is enabled so hard to
                                        tell for sure.
                            monitor: When MV = identity, same vertical shift,
                            does not take up horizontal of stereo display,
                            but seems more exagerated then when the headset
                            is on.
                            When MV = defaultState... plane is not square,
                            likely issue with default camera position
                    %}

% GL data structure for cross script access to OpenGL state variables
    global GL;

    % Default PTB setup, 2 = full features. The input of 2 means: 
    % execute the AssertOpenGL command, execute KbName('UnifyKeyNames') routine, AND unifies the color mode and switches from 0-255 to 0.0-1.0 - color part only impacts the current function or script, 
    % not ones that are called -- not the color clamping only applies to ptb
    % functions.
    PsychDefaultSetup(2);

    screenid = max(Screen('Screens')); % Find the screen 
    movieframe_n = 1;

    try
        % Setup Psychtoolbox for OpenGL 3D rendering support and initialize the
        % mogl OpenGL for Matlab/Octave wrapper:
        InitializeMatlabOpenGL;
        PsychDebugWindowConfiguration([],.8);
        % Not sure why this needs to be called but program complains if
        % prepareConfiguration is not called before psychvrhmd
        PsychImaging('PrepareConfiguration'); %first to prepare the configuration phase!
        
        if streamToHMDFlag
        % LowPersistence TimeWarp FastResponse are all redundant for CV1. CV1
        % does these things automatically according to PTB Documentation
            hmd = PsychVRHMD('AutoSetupHMD', 'Tracked3DVR', 'LowPersistence TimeWarp FastResponse DebugDisplay', 0);

            if isempty(hmd)
                fprintf('No VR-HMD available, giving up!\n');
                return;
            end
        end
        %{
            The rect returned by PsychImaging(‘Openwindow’) and Screen(‘Rect’), as
            well as sizes returned by Screen(‘WindowSize’) define the net useable
            size of the window in display pixels. ***It is affected by all kind of
            PsychImaging operations, e.g., selection of stereo modes, high bit depth
            modes etc.***, but also by scaling on Retina displays in high res mode.}%
        %}

        % [windowPtr,rect]=Screen(‘OpenWindow‘,windowPtrOrScreenNumber [,color] [,rect] [,pixelSize] [,numberOfBuffers] [,stereomode] [,multisample][,imagingmode][,specialFlags][,clientRect]
        multisample = 8;
        stereomode = 4; % From PsychImaging Doc: If possible on your setup and OS,
                        % rather use a single window, spanning both stereo display outputs, and use
                        % stereomode 4 or 5 to display dual-display stereo.
        [win, winRect] = PsychImaging('OpenWindow', screenid, 0, [], [], [], stereomode, multisample); %Do we need multisample?
        if streamToHMDFlag
            hmdinfo = PsychVRHMD('GetInfo', hmd);
            doSeparateEyeRender = hmdinfo.separateEyePosesSupported;        % Yes: Ask the driver if separate passes would be beneficial, and
            % use them if the driver claims it is good for us: turns out
            % CV1 does not use seperate eye render, I think this is just an
            % optimization thing and isn't a big deal.
            

            if doSeparateEyeRender
                fprintf('Will use separate eye render passes for enhanced quality on this HMD.\n');
            else
                fprintf('Will not use separate eye render passes, because on this HMD they would not be beneficial for quality.\n');
            end
        end


        % Setup the OpenGL rendering context of the onscreen window for use by
        % OpenGL wrapper. After this command, all following OpenGL commands will
        % draw into the onscreen window 'win':
        Screen('BeginOpenGL', win);

        % Set viewport properly(according to a ptb VR demo):
        glViewport(0, 0, RectWidth(winRect), RectHeight(winRect));

        % Retrieve and set camera projection matrix for optimal rendering on the HMD:
        %{
            ‘projL’ is the 4x4 OpenGL projection matrix for the left eye rendering.
            ‘projR’ is the 4x4 OpenGL projection matrix for the right eye rendering.
            Please note that projL and projR are usually identical for typical rendering
            scenarios.
        %}
        
        if streamToHMDFlag
            % Doc for GetStaticRenderParameters is short: http://psychtoolbox.org/docs/PsychOculusVRCore-GetStaticRenderParameters
            %projL          %projR (I guess)
            [projMatrix{1}, projMatrix{2}] = PsychVRHMD('GetStaticRenderParameters', hmd);
        else
            glMatrixMode(GL.PROJECTION);
            %glLoadIdentity;
            load('dsHmdInfoCV1.mat');
            load('DefaultHMDParameters.mat')
            projMatrix = ds.projMatrix; % ds is a struct that comes from dsHmdInfoCV1.mat
            defaultState_ = defaultState; % defaultState is a struct that comes from DefaultHMDParameters.mat
        end
        % Setup modelview matrix: This defines the position, orientation and
        % looking direction of the virtual camera:
        glMatrixMode(GL.MODELVIEW);
        glLoadIdentity;

        % Set background clear color to 'black' (R,G,B,A)=(0,0,0,0):
        glClearColor(.5,.5,.5,0);

        % Clear out the backbuffer: This also cleans the depth-buffer for
        % proper occlusion handling: You need to glClear the depth buffer whenever
        % you redraw your scene, e.g., in an animation loop. Otherwise occlusion
        % handling will screw up in funny ways...
        glClear;

        centerDotxyz   = [0;0;-1];
        centerDotSize  = 10;
        centerDotColor = [1, 1, 1]; 

        drawCenterDot = glGenLists(1);
        glNewList(drawCenterDot, GL.COMPILE);
        moglDrawDots3D(win, centerDotxyz, centerDotSize, centerDotColor, [], 1);
        glEndList;
        
        drawPlaneList = createPlane(GL);
        
        Screen('EndOpenGL', win);
        while ~KbCheck
            % Camera position when using head tracking + HMD: (according to SuperShapeDemo.m)
            globalPos = [0, 0, 0]; % x,y,z  % in meters - this is the starting position of head... should be at origin
            heading = 0; % yaw.
            globalHeadPose = PsychGetPositionYawMatrix(globalPos, heading); % just an identity matrix to start with % initialize observer's start position to the default camera position specified above
            if streamToHMDFlag
                state = PsychVRHMD('PrepareRender', hmd);
                 defaultIPD = 0;%0.064; % in m
                viewingDistance = .80;%0.45;
                modelViewDataLeft = [1 0 0 -defaultIPD/2; 0 1 0 0; 0 0 1 viewingDistance; 0 0 0 1]; %+/- IPD/2
                modelViewDataRight = [1 0 0 defaultIPD/2; 0 1 0 0; 0 0 1 viewingDistance; 0 0 0 1];
                state.modelView = {modelViewDataLeft, modelViewDataRight};
                state.modelView{1}(13) = defaultIPD/2;
                state.modelView{1}(15) = -viewingDistance;
                state.modelView{2}(13) = -defaultIPD/2;
                state.modelView{2}(15) = -viewingDistance;
            else
                % When both defaultIPD and Viewing Distance are set to 0
                % both left and right eye 
               
                  
            end
            
            for renderPass = 0:1
                % Selected 'view' to render (left- or right-eye) equals the renderPass,
                % as order of rendering does not matter in this mode:
                Screen('SelectStereoDrawbuffer', win, renderPass); % You need to call this command after each
                                                                   % Screen(‘Flip’) command or after drawing to an offscreen window again in order to
                                                                   % reestablish your selection of draw buffer, otherwise the results of drawing
                                                                   % operations will be undefined and most probably not what you want.
                Screen('BeginOpenGL', win);
                
                if streamToHMDFlag
                    modelView = state.modelView{renderPass + 1};
                else
                    if renderPass == 0
                        modelView = defaultState_.modelViewDataLeft;
                    elseif renderPass == 1
                        modelView = defaultState_.modelViewDataRight;
                    end
                    
                end
                % Set per-eye projection matrix: This defines a perspective projection,
                % corresponding to the model of a pin-hole camera - which is a good
                % approximation of the human eye and of standard real world cameras --
                % well, the best aproximation one can do with 2 lines of code ;-)
                glMatrixMode(GL.PROJECTION);
                if streamToHMDFlag
                    glLoadMatrixd(projMatrix{renderPass+1});
                    %glLoadIdentity;
                else
                    %glLoadIdentity;
                    glLoadMatrixd(projMatrix{renderPass+1});
                end
                % Setup camera position and orientation for this eyes view:
                glMatrixMode(GL.MODELVIEW);
                if streamToHMDFlag
                    glLoadMatrixd(modelView);
                    %glLoadIdentity;
                else
                    glLoadMatrixd(modelView);
                    %glLoadIdentity;
                end
                
                glClear;

                %glCallList(drawCenterDot);
                glCallList(drawPlaneList);
                % Finish OpenGL rendering into PTB window. This will switch back to the
                % standard 2D drawing functions of Screen and will check for OpenGL errors.
                % [Basically we end th opengl context to use other Screen
                % functions like 'Flip'
                Screen('EndOpenGL', win);
                Screen('Flip', win, [], 1);
            end
        
         if movieframe_n < 10  % save first frame of each movie (not the rest for speed)
                 rect = [];
                 %added to save movie clip
                M = Screen('GetImage', win,rect,[],0,1);
                imwrite(M, fullfile(pwd, 'Movie', strcat('hmdFrame', num2str(movieframe_n),'.png')));
                movieframe_n = movieframe_n + 1;
         end

        % Next frame ...
        end

    catch
        sca;
        psychrethrow(psychlasterror);
        %When modelView set to identity dot appears vertically shifted up *in seemingly equal magnitude to our vertical shift* of each stereo display.
    end
