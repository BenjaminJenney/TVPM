function drawFullScreenMaskForTVPMCD(window, mask)
    % Draws the swiss cheese mask. This simulates having a blindfold on
    % with holes in it.
    %
    % INPUTS:
    %   * window is an integer that ids the display, neccessary for
    %     psychtoolbox calls to Screen.
    %   * mask is the full window mask, left eye and right eye.
    %
    % Screen is a psychtoolbox function, the mask is drawn through
    % psychtoolbox and is not actually projected into the virtual world,
    % but onto the physical screens inside the headset. This is essential
    % for the mask over your eyes effect we want to achieve but has some
    % unfortunate side effects: we can't actually record or print bitmaps
    % of the stimulus with the mask on.
    
    glPushMatrix;
    glLoadIdentity;
    Screen('EndOpenGL', window);
    Screen('DrawTexture', window, mask, [], [], [], [], []);
    Screen('BeginOpenGL', window);
    glPopMatrix;
end