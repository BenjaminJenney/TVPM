function [pa, kb] = GetResponse(pa,ds,kb)

% Adjusts the position of a "paddle" to the expected interception location
% with the approaching dot

% CSB: and also sets the camera matrix for each trial based on awsz keys...
% June 7th 2018. KeyCodes:
%{
kb.translateViewLeftKey = KbName('a'); % translates the cameras left
kb.translateViewRightKey = KbName('s'); % translates the cameras right
kb.translateViewForwardKey = KbName('w'); % translates the cameras toward the stimulus
kb.translateViewBackKey = KbName('z'); % translates the cameras back away from the stimulus
kb.translateViewUpKey = KbName('e'); % translates the cameras up 
kb.translateViewDownKey = KbName('x'); % translates the cameras down

kb.rotateViewRightKey = KbName('r'); % rotates the cameras up 
kb.rotateViewLeftKey = KbName('q'); % rotates the cameras down
%}

shiftAmt = .01; % CSB: for translating camera ridigly. Increase for more "sensitivity" to key press
theta = .01*pi; % CSB: for rotating camera .  Increase for more "sensitivity" to key press

rotationMatrixRollLeft = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
rotationMatrixRollRight = [cos(-theta) -sin(-theta) 0; sin(-theta) cos(-theta) 0; 0 0 1];

rotationMatrixYawRight = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta);];
rotationMatrixYawLeft = [cos(-theta) 0 sin(-theta); 0 1 0; -sin(-theta) 0 cos(-theta);];
 
[kb.keyIsDown, kb.secs, kb.keyCode] = KbCheck(-1);



if kb.keyIsDown
    pa.rotationSpeed = pa.acceleratePaddle + pa.rotationSpeed;
    
    if kb.keyCode(kb.spacebarKey) % response has been entered
        kb.responseGiven = 1;
        % record response and current parameters
        if pa.trialNumber > 0 && pa.trialNumber <= pa.nTrials % Do not write out data for the first and last dummy, trials
            pa.responseTime = kb.secs - pa.trialOnset;
            pa.currentTime = ds.vbl - pa.experimentOnset;
            pa.response(pa.trialNumber,:) = [pa.trialNumber, pa.xSpeed, pa.zSpeed, pa.targetContrast, pa.paddleAngleInitial, pa.paddleAngle(pa.trialNumber), ds.trackingLag, pa.timeToPaddle, pa.speedUpFlag, pa.responseTime, pa.currentTime];
            %%% trial number, speed of trajectory in x (includes
            %%% angle); speed of trajectory in z; target contrast value;
            %%% starting and final paddle angles (observer setting); lag in screen updating with head motion; time it would take the target to reach the paddle; amount the visual feedback was sped up if needed; time it took to adjust and lock
            %%% in response; current time since the onset of the experiment
            %%% (to link with head motion data later)
        end
        
    elseif kb.keyCode(kb.paddleClockWiseKey)  %% left arrow
        pa.paddleAngle(pa.trialNumber) = mod(pa.paddleAngle(pa.trialNumber) - pa.rotationSpeed,360);
    elseif kb.keyCode(kb.paddleCounterClockWiseKey)  %% right arrow
        pa.paddleAngle(pa.trialNumber) = mod(pa.paddleAngle(pa.trialNumber) + pa.rotationSpeed,360);
    elseif kb.keyCode(kb.translateViewLeftKey) % CSB. these control modelView, camera position, need to add orientation too... June 7th 2018
        pa.modelView(1,4) = pa.modelView(1,4) + shiftAmt; % translation x
    elseif kb.keyCode(kb.translateViewRightKey)
        pa.modelView(1,4) =  pa.modelView(1,4) - shiftAmt;
    elseif kb.keyCode(kb.translateViewForwardKey)
        pa.modelView(3,4) = pa.modelView(3,4) + shiftAmt;
    elseif kb.keyCode(kb.translateViewBackKey)
        pa.modelView(3,4) = pa.modelView(3,4) - shiftAmt;  % translation z
    elseif kb.keyCode(kb.translateViewUpKey)
        pa.modelView(2,4) = pa.modelView(2,4) - shiftAmt;
    elseif kb.keyCode(kb.translateViewDownKey)
        pa.modelView(2,4) = pa.modelView(2,4) + shiftAmt;  % translation y
    elseif kb.keyCode(kb.rotateViewRightKey)
        pa.modelView(1:3,1:3) = pa.modelView(1:3,1:3)*rotationMatrixYawRight; % rotation about y axis
    elseif kb.keyCode(kb.rotateViewLeftKey)
        pa.modelView(1:3,1:3) = pa.modelView(1:3,1:3)*rotationMatrixYawLeft;  
    end
else
    pa.rotationSpeed = pa.shiftPaddle;
end

