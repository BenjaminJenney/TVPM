function drawFixationDot(GL, window, fixationDotDiameter, depthOfMiddlePlane)
    % drawFixationDot draws the fixation dot.
    % 
    % INPUTS:
    %   * GL, a struct that allows us to access OpenGL Variables
    %     across scripts. We need this to access the ModelView matrix values explicitly.
    %   * window is an integer that references the display. This integer id is given by
    %     psychtoolbox (for either the 2d monitor, or headset) and defined in SetupDisplay in the ds struct (ds.w).
    %     We need a reference to the display since we are using psychtoolbox's moglDrawDots3D for drawing the fixation dot.
    %   * fixationDotDiameter is a variable that controls the size of the fixation dot
    %   * depthOfMiddlePlane, the fixation dot is projected at the same
    %     depth as the middle plane
    % 
    % ATTN: CHECK WITH CHARLIE TO MAKE SURE THIS IS CORRECT: We multiply the the fixation dot position with the inverse of the
    % ModelView matrix as to prevent the dot from tracking with your head.
    %
    % About curModelViewNoTranslation: Charlie had to remove the translation column from the ModelView matrix; 
    % something about the translation column was screwing up the dots projection.
    % the fix is hacky but it works and we haven't yet narrowed down the
    % culprit for this strange behavior. It doesn't seem to effect any
    % other aspects of the code. 
    
    
    glPushMatrix;
    curModelViewNoTranslation1 = glGetFloatv(GL.MODELVIEW_MATRIX);
    curModelViewNoTranslation2 = [curModelViewNoTranslation1(1:4),curModelViewNoTranslation1(5:8),curModelViewNoTranslation1(9:12),[0 0 0 1]']; % CSB, 4/12/2022- had to use the GL modelview without translation to make the fixation pt at right depth. it was broken b4 using eye.modelView. figure out why
    moglDrawDots3D(window, inv(curModelViewNoTranslation2)*[[0 0 depthOfMiddlePlane], 1.0]', fixationDotDiameter, [1 0 0], [], 2); % drawing horizontal fixation dot
    % Check here for input description and functionality explanation of
    % moglDrawDots3D: https://raw.githubusercontent.com/Psychtoolbox-3/Psychtoolbox-3/beta/Psychtoolbox/PsychOpenGL/moglDrawDots3D.m
    
    glColor3dv([1,1,1]); % Mogl draw dots calls glColor3dv to set the color of the fixation dot, this to change the opengl color state, have to set it back to white to not turn all textures red. Possibly a bug with MoglDrawDots
    
    glPopMatrix;
    
    
end