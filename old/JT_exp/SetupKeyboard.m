function kb = SetupKeyboard(ds,DEBUG_FLAG)
% Setup universal Mac/PC keyboard and keynames - no need for hard-coding!

if ~DEBUG_FLAG
    HideCursor(ds.screenId);
    ListenChar(2); % Stop making keypresses show up in the matlab scripts and command window - it's really annoying and can cause all sorts of problems!
end

KbName('UnifyKeyNames');

kb.headingLeftKey = KbName('LeftArrow'); % 'left' means heading left of target
kb.headingRightKey = KbName('RightArrow'); % 'right' means heading left of target
kb.spacebarKey = KbName('space'); % 'space' starts experiment or end break

kb.escapeKey = KbName('j'); % quits out of the experiment before its completion

% Initialize KbCheck
[kb.keyIsDown, kb.secs, kb.keyCode] = KbCheck(-1);
KbReleaseWait; % Make sure all keys are released
kb.keyWasDown= 0;

end