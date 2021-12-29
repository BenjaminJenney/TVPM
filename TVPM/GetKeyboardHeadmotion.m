function [pa, kb, eye, ds, state] = GetKeyboardHeadmotion(pa,ds,kb,eye,state)

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
   
    if kb.keyCode(kb.translateViewLeftKey) % CSB. these control modelView, camera position, need to add orientation too... June 7th 2018
        eye.modelView(1,4) = eye.modelView(1,4) + shiftAmt; % translation x
    elseif kb.keyCode(kb.translateViewRightKey)
        eye.modelView(1,4) =  eye.modelView(1,4) - shiftAmt;
    elseif kb.keyCode(kb.translateViewForwardKey)
        eye.modelView(3,4) = eye.modelView(3,4) + shiftAmt;
    elseif kb.keyCode(kb.translateViewBackKey)
        eye.modelView(3,4) = eye.modelView(3,4) - shiftAmt;  % translation z
    elseif kb.keyCode(kb.translateViewUpKey)
        eye.modelView(2,4) = eye.modelView(2,4) - shiftAmt;
    elseif kb.keyCode(kb.translateViewDownKey)
        eye.modelView(2,4) = eye.modelView(2,4) + shiftAmt;  % translation y
    elseif kb.keyCode(kb.rotateViewRightKey)
       
        eye.modelView(1:3,1:3) = eye.modelView(1:3,1:3)*rotationMatrixYawRight; % rotation about y axis
        %ds.projMatrix{1}(1:3,1:3) = ds.projMatrix{1}(1:3,1:3)*rotationMatrixYawRight; % rotation about y axis
        %ds.projMatrix{2}(1:3,1:3) = ds.projMatrix{2}(1:3,1:3)*rotationMatrixYawRight; % rotation about y axis
   
        %state.globalHeadPose(1:3,1:3) = state.globalHeadPose(1:3,1:3)*rotationMatrixYawRight;
        
    elseif kb.keyCode(kb.rotateViewLeftKey)
        
        eye.modelView(1:3,1:3) = eye.modelView(1:3,1:3)*rotationMatrixYawLeft;  
        %ds.projMatrix{1}(1:3,1:3) = ds.projMatrix{1}(1:3,1:3)*rotationMatrixYawLeft; % rotation about y axis
   %ds.projMatrix{2}(1:3,1:3) = ds.projMatrix{2}(1:3,1:3)*rotationMatrixYawLeft; % rotation about y axis
   
   %state.globalHeadPose(1:3,1:3) = state.globalHeadPose(1:3,1:3)*rotationMatrixYawLeft;
   
    end
else
    pa.rotationSpeed = pa.shiftPaddle;
end

