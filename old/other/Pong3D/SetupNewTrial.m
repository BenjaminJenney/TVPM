function [ds, pa, kb, oc] = SetupNewTrial(ds, pa, kb, oc)

% do this before every trial to get the next trial's starting info (e.g.,
% target motion parameters)


if pa.trialNumber>0 % if it's past the first trial, wait for the Up Arrow key to be pressed after the feedback to initiate a new trial
    [kb.keyIsDown, kb.secs, kb.keyCode] = KbCheck(-1);
    
    if kb.keyIsDown
        
        if kb.keyCode(kb.reorientKey)
  
            pa.trialNumber = pa.trialNumber + 1;
            
            if pa.trialNumber==pa.nTrials+1 % we want to repeat the first trial, in case the participant didn't see and just put a random setting - we'll throw out the first trial
                pa.thisTrial = 1;
                if strcmp(ds.condition,'mono')
                    pa.randeye(pa.trialNumber) = pa.randeye(1);
                end
            else
                pa.thisTrial = pa.trialNumber;
            end
            
            if pa.trialNumber <= pa.nTrials+1
                
                pa.xSpeed = pa.fullFactorial(pa.thisTrial,2);
                pa.zSpeed = pa.fullFactorial(pa.thisTrial,3);
                
                pa.timeToPaddle = (pa.paddleOrbitShift-pa.paddleHalfWidth*3.8) ./ sqrt(pa.xSpeed.^2 + pa.zSpeed.^2);
                
                % speed up the really really slow feedback - this happens
                % on maybe 10% of trials?
                if pa.timeToPaddle > 2
                    pa.speedUpFlag = pa.timeToPaddle / 2; % just make the feedback no longer than 2 s long
                else
                    pa.speedUpFlag = 1;
                end
                pa.timeToPaddle = (pa.paddleOrbitShift-pa.paddleHalfWidth*3.8) ./ sqrt((pa.speedUpFlag*pa.xSpeed).^2 + (pa.speedUpFlag*pa.zSpeed).^2);
                
                if pa.trialNumber==1
                    pa.paddleAngle(pa.thisTrial) = pa.paddleAngle; % deg
                    pa.paddleAngleInitial = pa.paddleAngle(pa.thisTrial);
                else
                    pa.paddleAngle(pa.thisTrial) = randi(361)-1; % deg  - assigns a random integer angle between 0 and 360
                    pa.paddleAngleInitial = pa.paddleAngle(pa.thisTrial);
                end
                
                pa.targetContrast = pa.fullFactorial(pa.thisTrial,1);
                
                kb.responseGiven = 0;
                pa.feedbackGiven = 0;
    
                pa.trialOnset = ds.vbl; % the trial starts....NOW
                
            end
        end
    end
    
else % if it is just the first trial, start right up after the subject presses 'space' to confirm that he/she is ready to go
    
    pa.trialNumber = pa.trialNumber + 1;
    
    if pa.trialNumber==pa.nTrials+1 % we want to repeat the first trial, in case the participant didn't see and just put a random setting - we'll throw out the first trial
        pa.thisTrial = 1;
        if strcmp(ds.condition,'mono')
            pa.randeye(pa.trialNumber) = pa.randeye(1);
        end
    else
        pa.thisTrial = pa.trialNumber;
    end
    
    if pa.trialNumber <= pa.nTrials+1
        
        pa.xSpeed = pa.fullFactorial(pa.thisTrial,2);
        pa.zSpeed = pa.fullFactorial(pa.thisTrial,3);
        
        pa.timeToPaddle = (pa.paddleOrbitShift-pa.paddleHalfWidth*3.8) ./ sqrt(pa.xSpeed.^2 + pa.zSpeed.^2);

                if pa.timeToPaddle > 2
                    pa.speedUpFlag = pa.timeToPaddle / 2; % just make the feedback no longer than 2 s long
                else
                    pa.speedUpFlag = 1;
                end
        
        pa.timeToPaddle = (pa.paddleOrbitShift-pa.paddleHalfWidth*3.8) ./ sqrt((pa.speedUpFlag*pa.xSpeed).^2 + (pa.speedUpFlag*pa.zSpeed).^2);
        
        
        pa.paddleAngle(pa.thisTrial) = pa.paddleAngle; % deg
        pa.paddleAngleInitial = pa.paddleAngle(pa.thisTrial);
        
        pa.targetContrast = pa.fullFactorial(pa.thisTrial,1);
        
        kb.responseGiven = 0;
        pa.feedbackGiven = 0;

        pa.trialOnset = ds.vbl; % the trial starts....NOW
    end
end

