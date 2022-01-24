function [vX, vY] = opticFlow(ds, pa, xpos_deg, ypos_deg, plane, numPoints, modelView, t, omegaY, T) % add arugment for plane centers and input below
    kernel = .2 * plane.width_m; % this is to make up for the edge nan problem when interpolating
    Y = rand(numPoints,1)* plane.height_m - (plane.height_m / 2); % ACTION ITEM: add points along the edges, currently not gaurenteed to always work.
    X = rand(numPoints,1)* (plane.width_m + kernel) - ((plane.width_m + kernel) / 2) + plane.offset; % change width/2 to something that centers planes at their 3 respective centers
    Z = ones(numPoints,1)* -plane.depth_m;
    W = ones(numPoints,1);
    %keyboard
    initPos = [X,Y,Z,W]'; % initial dot positions

    TCur = T(:,t) .* ds.ifi; % current translation velocity, m/s -> m/f
    omega = [0; omegaY(t) .* ds.ifi; 0]; % current rotation velocity, rad/s -> rad/f

    %rotMatY = [cosd(h) 0 sind(h); 0 1 0; -sind(h) 0 cosd(h)]; % Don't remove
    %curPositions  = rotMatY*initPos; % don't remove

    curPositions = modelView'*initPos; % update dot positions based on current observer position and viewing angle based on gluLookAt 

    x = ds.focalLength * curPositions(1,:) ./ curPositions(3,:);
    y = ds.focalLength * curPositions(2,:) ./ curPositions(3,:);
    p = 1 ./ curPositions(3,:);

    %keyboard
    for i = 1:size(x,2)
        A = [-ds.focalLength, 0, x(i); 0, -ds.focalLength, y(i)];
        B =  [(x(i)*y(i))/ds.focalLength, -(ds.focalLength+x(i)^2/ds.focalLength), y(i);...
            ds.focalLength+y(i)^2/ds.focalLength, -(x(i)*y(i))/ds.focalLength, -x(i)];
        velocityField1(:,i) = p(i).*A*TCur+B*omega; % m/f on imaging plane
    end

    velocityField2 = (velocityField1.*pa.nFrames)./pa.stimulusDuration_sec; % m/s

    velocityField = 2.*atand(velocityField2./(2.*ds.focalLength)); % m/s -> deg/s.  % figure out issue with this units - ask david and bas

    xdeg = 2.*atand(x./(2.*ds.focalLength)); % m -> deg
    ydeg = 2.*atand(y./(2.*ds.focalLength)); % m -> deg

    %% interpolate flow field at desired retinal locations

    % specify the lattice of image poisitons we want to sample the flow
    % field at (in degrees):
    %numPlanes = 3

    %these are the 2d retinal positions, they are passed to optic flow, and 
    %xPos = linspace(-ds.hFOV_perPersonAvg/2,ds.hFOV_perPersonAvg/2,numPoints/6); % something like this but replace with correct bounds for each plane
    %yPos = linspace(-ds.vFOV_perPersonAvg/2,ds.vFOV_perPersonAvg/2,numPoints/6);

    % xPos = linspace(-4,1,numPoints/6); % temporary hardcoded
    % yPos = linspace(-2,12,numPoints/6);

    %coords = allComb(xPos,yPos);

    coords = [xpos_deg,ypos_deg];
    % interpolte horizontal and vertical components of optic flow
    % separately (two scalar fields)
    %keyboard
    vX = griddata(xdeg,ydeg,velocityField(1,:),coords(:,1),coords(:,2));
    vY = griddata(xdeg,ydeg,velocityField(2,:),coords(:,1),coords(:,2));
    
    
%         figure; subplot(2,1,1); quiver(xdeg,ydeg,velocityField(1,:),velocityField(2,:),0);
%         subplot(2,1,2); quiver(coords(:,1),coords(:,2),vX,vY,0) % plot to check
%         hold on; scatter(coords(:,1),coords(:,2))
    
  end