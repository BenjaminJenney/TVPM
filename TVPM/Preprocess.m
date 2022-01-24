function [vXw, vYw] = Preprocess(ds, pa,shape, GL)

%hGratings = cell(360,1);
numPlanes = shape.plane.numPlanes;
sf =.1;% cycles per deg
sf = (ds.deg_per_px).*sf;%1/px_per_deg %convert from deg to
%af = 2*pi*sf; % angular freq.
nFrames = pa.nFrames; %
ndp = shape.disk.numDisksPerPlane;
numRenderPass = 2;
%hPhase = single(zeros(pa.nTrials,nFrames,2,numPlanes,180));%trialNum,t,renderPass+1,i,j
%vPhase = single(hPhase);
%ifi = ds.ifi;
tBase = linspace(0, pa.stimulusDuration_sec, nFrames);
pa.rotVelocityRad = deg2rad(0);
ly = pa.pursuitTargetDist; % intial distance (depth in Z) to fixation pt at time 0
s = pa.transSpeed; %m/s
vXw = single(zeros(pa.nTrials, pa.nFrames, shape.plane.numPlanes, shape.disk.numDisksPerPlane));
vYw = single(zeros(size(vXw)));
planeDist_m = shape.plane.depths_m;

for ii = 1:pa.nTrials
    for dd = 1:pa.numHeadingAngles
        if pa.headingAngleVecRandThisRun(ii) == dd
            pa.heading(ii) = pa.headingAngles(dd);
        end
    end
end

tstart = tic;
for ii = 1:pa.nTrials %pa.nTrials = 60 ; BJ: Number of frames per trial is dependent on pa.stimulusDuration as well as other constants see conditions on ln no. 723 of tvpmexpirement
    %time-varying heading based on fixating observer
    theta0 = pa.heading(ii); % initial heading angle in degrees
    
    rho = -1 * atand( ( (ly.*tand(theta0)) - ( ( (ly./cosd(theta0)) - s.*tBase).*sind(theta0) ) ) ...
        ./ ( ( (ly./cosd(theta0)) - s.*tBase) .*cosd(theta0) )  );  % deg
    
    dRho   = gradient(rho);
    omegaY = deg2rad(dRho)./ds.ifi; % rad/s
    h = rho+theta0; % deg, current heading angle
    
    T = [s.*sind(h); zeros(1,length(h)); s.*cosd(h)]; % heading angle = theta0+rho, m/s
    for t = 1:nFrames %90*.8
        %for renderPass = 0:1   
            elapsedTime = t*(1/90); % elapsed time in seconds since trial start (pressing space bar)
            xDotDisplacement = pa.transSpeed.*sind(-1.*pa.heading(pa.trialNumber)).*elapsedTime; % CSB: should be in meters. check. % ATTN: take negative of heading angle, so it displays correctly. otherwise displayed heading is inverse of ground truth
            % pa.rotVelocityRadPursuit(pa.trialNumber) = deg2radtrialNumber)).*elapsedTime; % CSB: should be in meters. check. % ATTN: take negative of heading angle, so it displays correctly. otherwise displayed heading is inverse of ground truth
            zDotDisplacement = pa.transSpeed.*cosd(-1.*pa.heading(pa.trialNumber)).*elapsedTime;  % JT: make sure this trig is right.
            
            glLoadIdentity;
            %                preMV = glGetFloatv(GL.MODELVIEW_MATRIX);
            %                preMV = reshape(preMV, [4,4]);
            
            gluLookAt(-xDotDisplacement,0,-zDotDisplacement,0,0,-pa.cubeWidth/4,0,1,0);
            %glTranslatef(-ds.planeWidths_m(2) + xDotDisplacement, 0.0, 0.0 + zDotDisplacement);
            curMV = glGetFloatv(GL.MODELVIEW_MATRIX);
            
            curMV = reshape(curMV, [4,4]);
            
            % for each render pass create the current texture frame for
            % each disk
            
            %                1) in preprocess.m, make optic flow for every frame, trial, eye.
            % 2) in preprocess.m, make vertical and horizontal phases for every frame, trial, eye  (these will be indices to index into grating array in runexpt..)
            % 3) in preprocess.m, make a 360-long struct containing one grating shifted by every possible phase.
            % 4) in RunTVPMExperiment.m, for every frame, make the image for every plaid by indexing into the right phase for each component grating, then transpose one of them that's vertical, and sum them together, and normalize. then repmat to get three channels and permute [if this is slowing everything down, then make the pre-generated grating 3 dimensional first]
            % 5) bind textures and render.
            
%             hPhaseTmp = 0;
%             vPhaseTmp = 0;
            for i = 1:numPlanes
                [vX, vY] = opticFlow(ds, pa, shape.disk.xpos_deg{i}, shape.disk.ypos_deg{i}, shape.plane.planes(i), 1050, curMV, t, omegaY, T); % deg/s

                vXw_s = 2.*-planeDist_m(i).*tand(vX./2); % m/s in world coords
                vXw(ii, t, i, :) = vXw_s.*(pa.stimulusDuration_sec./pa.nFrames); % m/f in world coords
                
                vYw_s = 2.*-planeDist_m(i).*tand(vY./2); % m/s in world coords
                vYw(ii, t, i, :) = vYw_s.*(pa.stimulusDuration_sec./pa.nFrames); % m/f in world coords
%                 v = (1/ds.deg_per_px).*vY; % vX and vY are vectors of optic flow for all 180 dots per plane
%                 u = (1/ds.deg_per_px).*vX;
%                 
%                 hPhaseTmp = hPhaseTmp + ((2 .* pi .* (v .* (1./1000) .* sf)) ./ ifi);
%                 vPhaseTmp = vPhaseTmp + ((2 .* pi .* (u .* (1./1000) .* sf)) ./ ifi);
%                 
%                 hPhase(ii,t,renderPass+1,i,:) = single(round(wrapTo360(rad2deg(hPhaseTmp)))); %index for gratings
%                 vPhase(ii,t,renderPass+1,i,:) = single(round(wrapTo360(rad2deg(vPhaseTmp)))); %index for gratings  
            end           
        %end      
    end    
end
tElapsed = toc(tstart)

%for i = 1:360; imagesc(hGratings{i}); pause(ds.ifi); colormap grey; shg; end
% for i = 1:360
%     hGratings{i} = (.5 + .5 * cos(af * shape.disk.texture.y - deg2rad(i)));
% end

 
   %keyboard


%plot example grating over time to make sure tvpm looks correct
% figure;
% for ii = 1:60
%     vCurPhase = vPhase(1,ii,1,1,50); if vCurPhase == 0; vCurPhase = 1; end
%     hCurPhase = hPhase(1,ii,1,1,50); if hCurPhase == 0; hCurPhase = 1; end
%     image = .5.*hGratings{vCurPhase} + .5.*hGratings{hCurPhase}';
%     imagesc(image);
%     colormap gray;
%     pause(1/90);
% end

end