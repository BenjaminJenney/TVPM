function [pa, kb] = GetResponse(pa,ds,kb)
% Saves out the last key pressed in the ITI.
 
[kb.keyIsDown, kb.secs, kb.keyCode] = KbCheck(-1);
kb.answer = -1; %default value for kb.answer signaling there was no response

if kb.keyIsDown
   % pa.rotationSpeed = pa.acceleratePaddle + pa.rotationSpeed;
    if kb.keyCode(kb.headingLeftKey) % left key
        kb.responseGiven = 1;
        kb.answer = 0; %% answer is 0 for left
        % record response and current parameters
            pa.responseTime = kb.secs - pa.trialOnset;
            pa.currentTime = ds.vbl - pa.experimentOnset;
            pa.response(pa.trialNumber,:) = [pa.trialNumber, pa.responseTime, pa.currentTime];
            %%% trial number, response time, current time.
            %%% (to link with head motion data later) 
    elseif kb.keyCode(kb.headingRightKey) %% right key
        kb.responseGiven = 1;
        kb.answer = 1; % answer is 1 for right
        % record response and current parameters
            pa.responseTime = kb.secs - pa.trialOnset;
            pa.currentTime = ds.vbl - pa.experimentOnset;
            pa.response(pa.trialNumber,:) = [pa.trialNumber, pa.responseTime, pa.currentTime];
            %%% trial number, response time, current time.
            %%% (to link with head motion data later) 
    elseif kb.keyCode(kb.headingLeftKey) == 0 && kb.keyCode(kb.headingRightKey) == 0 && sum(kb.keyCode) > 0 % any key that is not right or left 
        kb.responseGiven = 1;
        kb.answer = -1;
        % record response and current parameters
        pa.responseTime = kb.secs - pa.trialOnset;
        pa.currentTime = ds.vbl - pa.experimentOnset;
        pa.response(pa.trialNumber,:) = [pa.trialNumber, pa.responseTime, pa.currentTime];
    end
end

