function kb = SetupKeyboard()

% Setup universal Mac/PC keyboard and keynames - no need for hard-coding!

%ListenChar(2); % Stop making keypresses show up in the matlab scripts and
%command window - it's really annoying and causes all sorts of problems! %CSB debug
%HideCursor; %CSB debug

KbName('UnifyKeyNames');

% CSB: I will attempt to use aswz keys to control the camera matrices when
% there is no hmd attached (using monitor). June 5, 2018
kb.translateViewLeftKey = KbName('a'); % translates the cameras left
kb.translateViewRightKey = KbName('s'); % translates the cameras right
kb.translateViewForwardKey = KbName('w'); % translates the cameras toward the stimulus
kb.translateViewBackKey = KbName('z'); % translates the cameras back away from the stimulus
kb.translateViewUpKey = KbName('e'); % translates the cameras up 
kb.translateViewDownKey = KbName('x'); % translates the cameras down

kb.rotateViewRightKey = KbName('r'); % rotates the cameras up 
kb.rotateViewLeftKey = KbName('q'); % rotates the cameras down



kb.paddleClockWiseKey = KbName('LeftArrow'); % adjusts the paddle position clockwise along the invisible circular orbit
kb.paddleCounterClockWiseKey = KbName('RightArrow'); % adjusts the paddle position counter clockwise along the invisible circular orbit
kb.spacebarKey = KbName('space'); % intiates a new trial/locks in a paddle adjustment response

kb.reorientKey = KbName('UpArrow'); % resets the observer's view, in case of drift or other issues that cause the display to shift

kb.escapeKey = KbName('ESCAPE'); % quits out of the experiment before its completion

% Initialize KbCheck
[kb.keyIsDown, kb.secs, kb.keyCode] = KbCheck(-1);
KbReleaseWait; % Make sure all keys are released
kb.keyWasDown= 0;

end