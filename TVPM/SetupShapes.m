function shape = SetupShapes(ds, pa)
    InitializeMatlabOpenGL(1);
    global GL;
    Screen('BeginOpenGL', ds.w);
    type = GL.TEXTURE_2D;
    %% Set up texture and dimensions for the three planes
    
    planeTextureWidth_px   = (ds.hFOV_perPersonAvg/3) * ds.px_per_deg;
    planeTextureHeight_px  = (ds.vFOV_perPersonAvg) * ds.px_per_deg;
      shape.plane.TextureWidth_px  = planeTextureWidth_px;
    shape.plane.TextureHeight_px = planeTextureHeight_px;
    xHalfTextureSize = planeTextureWidth_px/2;
    yHalfTextureSize = floor(planeTextureHeight_px/2);

    [x,y] = meshgrid(-xHalfTextureSize+1:xHalfTextureSize,-yHalfTextureSize+1:yHalfTextureSize);

    midPlaneAlpha = ones(size(x));
    
    % All three plaid are painted with the same texture data
    % (staticPliadData),
    % except midPlane needs a hole in its alpha chanel (alpha chanel = 4th index)
    fixationDotSize_px = 15;
    shape.plane.textureData  = staticPlaidData(ds, pa.sf, x, y);
    midPlaneAlpha = min(midPlaneAlpha, sqrt(x.^2 + y.^2) > fixationDotSize_px);%Cut out hole for fixation dot
    midPlaneTexData = shape.plane.textureData;
    midPlaneTexData(4,:,:) = shiftdim(255 .* midPlaneAlpha', -1);
    midPlaneTexture = Texture(GL, GL_TEXTURE_2D, xHalfTextureSize*2, yHalfTextureSize*2, uint8(midPlaneTexData));

    nearPlaneTexture = Texture(GL, GL_TEXTURE_2D, xHalfTextureSize*2, yHalfTextureSize*2, uint8(shape.plane.textureData));
    farPlaneTexture  = Texture(GL, GL_TEXTURE_2D, xHalfTextureSize*2, yHalfTextureSize*2, uint8(shape.plane.textureData));
    shape.plane.textures = [nearPlaneTexture, midPlaneTexture, farPlaneTexture];
    % End texture for the three planes

    % Setup vertices for the three planes
    
    depths_m  = [-.5, -1, -2]; % near, mid, far.
    numPlanes = length(depths_m);
    widths_m  = zeros(1, numPlanes);
    heights_m = zeros(1, numPlanes);
     for i = 1:numPlanes
        widths_m(i)  = 2 * -depths_m(i) * tand((ds.hFOV_perPersonAvg/numPlanes)/2);
        heights_m(i) = 2 * -depths_m(i) * tand(ds.vFOV_perPersonAvg/2);
     end
     %keyboard;
    numVertices = 4;
    corners = [0 0;
        1 0;
        1 1;
        0 1];

    shape.plane.listIds = zeros(1, numPlanes);
    shape.plane.vertices = {};
    for i = 1:numPlanes

        ithPlaneVertices =[-widths_m(i)/2 -heights_m(i)/2 depths_m(i) ;... 
                            widths_m(i)/2 -heights_m(i)/2 depths_m(i) ;...
                            widths_m(i)/2  heights_m(i)/2 depths_m(i) ;...
                           -widths_m(i)/2  heights_m(i)/2 depths_m(i) ]';
        shape.plane.vertices{i} = ithPlaneVertices;
        %ds.initPositions{i} = opticFlow(shape.plane.widths_m(i), shape.plane.heights_m(i), shape.plane.depths_m(i), 50);
        shape.plane.listIds(i) = glGenLists(1);
        glNewList(shape.plane.listIds(i), GL.COMPILE);
            shape.plane.textures(i).bind;
            glBegin(GL.POLYGON);
            for j = 1:numVertices
                glTexCoord2dv(corners(j,:));
                glVertex3dv(ithPlaneVertices(:,j));
            end
            glEnd;
            shape.plane.textures(i).unbind;
        glEndList(); 
    end
    
    shape.plane.depths_m  = depths_m;
    shape.plane.widths_m  = widths_m;
    shape.plane.heights_m = heights_m; 
    shape.plane.near.width_m  = widths_m(1);
    shape.plane.near.height_m = heights_m(1);
    shape.plane.near.depth_m  = depths_m(1);
    shape.plane.mid.width_m   = widths_m(2);
    shape.plane.mid.height_m  = heights_m(2);
    shape.plane.mid.depth_m   = depths_m(2);
    shape.plane.far.width_m   = widths_m(3);
    shape.plane.far.height_m  = heights_m(3);
    shape.plane.far.depth_m   = depths_m(3);
    
    shape.plane.near.offset = -shape.plane.near.width_m; 
    shape.plane.mid.offset  = 0; 
    shape.plane.far.offset  = shape.plane.far.width_m;
    shape.plane.offsets_m   = [-shape.plane.near.width_m,0,shape.plane.far.width_m];
    shape.plane.offsets_deg = [-ds.hFOV_perPersonAvg/3, 0, ds.hFOV_perPersonAvg/3]; 
    shape.plane.numPlanes   = numPlanes;
    shape.plane.planes      = [shape.plane.near, shape.plane.mid, shape.plane.far];
    shape.plane.numVertices = numVertices;
    %% Set up texture and dimensions for the Disks and apertures
    
    diskTextureWidth   = 32;
    diskTextureHeight  = 32;
    halfDiskTexWidth   = diskTextureWidth/2;
    halfDiskTexHeight  = diskTextureHeight/2;
    
    [shape.disk.texture.x,shape.disk.texture.y] = meshgrid(-halfDiskTexWidth+1:halfDiskTexWidth,-halfDiskTexHeight+1:halfDiskTexHeight);
    
    diskSize_px               = 32;
    apertPlaneSize_px         = diskSize_px*2;
    diskSize_deg              = diskSize_px*ds.deg_per_px; %2.0714
    apertPlaneSize_deg        = apertPlaneSize_px*ds.deg_per_px; 
    shape.disk.size_px        = diskSize_px; %not totally sure what units this is. glPointSize draws a square with equal sides in pixels, supposedly.
    shape.disk.size_m         = (1/ds.pixelsPerM)*sqrt(diskSize_px^2 + diskSize_px^2);
    shape.disk.size_deg       = diskSize_deg;
    shape.disk.aptPlane_px    = apertPlaneSize_px;
    shape.disk.aptPlane_m     = (1/ds.pixelsPerM)*sqrt(apertPlaneSize_px^2 + apertPlaneSize_px^2);
    shape.disk.aptPlane_deg   = apertPlaneSize_deg;
    shape.disk.texture.width  = diskTextureWidth;
    shape.disk.texture.height = diskTextureHeight;
    
    shape.disk.numDisksPerPlane = 10;
    shape.disk.numDisks = shape.disk.numDisksPerPlane * shape.plane.numPlanes;
    shape.disk = diskPos(ds, shape.disk, shape.plane);
    
    sf        = pa.sf; %.1; % cycles per deg
    sf        = (ds.deg_per_px).*sf; %1/px_per_deg %convert from deg to 
    af        = 2*pi*sf;
    xGrating  = (.5 + .5 * sin(af * shape.disk.texture.x - 0)); % vertically oriented grating
    yGrating  = (.5 + .5 * cos(af * shape.disk.texture.y - 0)); % horiz grating
    plaidData = 255.*((xGrating + yGrating)/2);
    plaidData = repmat(plaidData,[ 1 1 3 ]);
    plaidData = permute(plaidData,[ 3 2 1 ]);
    
    textureid = glGenTextures(1);
    glBindTexture(type, textureid(1));
        glTexParameterfv(type, GL.TEXTURE_WRAP_S, GL.REPEAT);
        glTexParameterfv(type, GL.TEXTURE_WRAP_T, GL.REPEAT);
        glTexParameterfv(type, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
        glTexParameterfv(type, GL.TEXTURE_MIN_FILTER, GL.LINEAR);
        glTexImage2D(type, 0, GL.RGB, shape.disk.texture.width, shape.disk.texture.width, 0, GL.RGB, GL.UNSIGNED_BYTE, uint8(plaidData));
        glTexEnvi(GL.POINT_SPRITE, GL.COORD_REPLACE, GL.TRUE);
    glBindTexture(type,0);
   Screen('EndOpenGL', ds.w);
   
    %% Masks
    b = .067; % inter-pupillary distance in meters TODO: Figure out if there is a way to get IP from the headset instead of hardcoding.
    Z = shape.plane.depths_m(2); % depth of middle plane.
    d_deg = 2.*( atand((b./2)./(2.*Z)) - atand((-b./2)./(2.*Z)) ); % Should fix disparity of the mask.
    d_px  = abs(d_deg*ds.px_per_deg);

    % Get the size of the on screen window
    [screenXpixels, screenYpixels] = Screen('WindowSize', ds.w);
    white = WhiteIndex(0);
    black = BlackIndex(0);

    % Make a gaussian aperture with the "alpha" channel
    gaussDim = pa.apertureDia_px; %100 px? TODO: DOES NOT SEEM TO CHANGE SIZE. 
    gaussSigma = pa.apertureDia_px; %gaussDim /4; % csb: controls how big the apertures are
    s1 = screenXpixels;
    s2 = screenYpixels;
    [xm, ym] = meshgrid(-(s2/2)+1:s2/2, -(s1/2)+1:s1/2);
    %ym   = randi([-s2/2, s2/2], s2/2, s2/2);
    %gauss = exp(-(((xm .^2) + (ym .^2)) ./ (2 * gaussSigma^2)));

    f = abs(cos(1/15.*[-s2/2:s2/2])).^8; f2 = f'*f;
     gauss = f2(1:s1,1:s2);
    fixationHole_px = 15 ;
    gauss = 1-gauss;
    shape.mask.gauss = min((gauss), sqrt((xm).^2 + (ym).^2)>fixationHole_px);
    
    bag = ones(screenYpixels, screenXpixels) .* black;
    bag = repmat(bag,[ 1 1 3 ]);
    
    % LEFT EYE
    mask = ones(s1, s2, 1);
    mask(:, :, 2) = gauss; % Original Raised cosine calculation
    mask(:,:,2) = round(mask(:,:,2),3,'significant');
    endCols = mask(:,end-round(d_px/2)+1:end,2);
    mask(:,:,2) = 1;
    mask(:,1:round(d_px/2),2) = endCols;
    bag(:,:,4) = circshift(mask(:,:,2)', 50);
    shape.mask.fullWindowMaskLeftEye = Screen('MakeTexture', ds.w, bag);
    
    %RIGHT EYE
    mask(:, :, 2) = gauss; % restart Alpha with original raised cosine calc.
    frontCols = mask(:,1:end-round(d_px/2)+1,2); % Take the columns up to end-disparity shift)
    mask(:,:,2) = 1; % reset Alpha to completely opaque
    mask(:,round(d_px/2):end,2) = frontCols; % paste to the end of the Alpha chanel
    bag(:,:,4) = circshift(mask(:,:,2)', 50);
    shape.mask.fullWindowMaskRightEye = Screen('MakeTexture', ds.w, bag);

    %ds.masktex = Screen('MakeTexture', ds.w, mask);

    %black = ones(screenYpixels, screenXpixels) .* black;
    %black(:, : 4) = mask;
     
    %bag = permute(bag,[ 3 2 1 ]);
    
    %bag(628:668,:,4) = -1;
    % Make a grey texture to cover the full window
    %shape.mask.fullWindowMask = Screen('MakeTexture', ds.w, bag);

    % Make coordinates in which to draw the apertures into our full screen mask
    [xg, yg] = meshgrid(-gaussDim:gaussDim, -gaussDim:gaussDim); % csb: and this controls how many apertures there are
    spacing = gaussDim * 2; % csb: i think this controls how far apart the apertures ares
    xg = xg .* spacing + screenXpixels / 2;
    yg = yg .* spacing + screenYpixels / 2; % hole y position%
    xg = reshape(xg, 1, numel(xg));
    yg = reshape(yg, 1, numel(yg));


    % Make the destination rectangles for the gaussian apertures
    shape.dstRects = nan(4, numel(xg));
    for i = 1:numel(xg)
    shape.dstRects(:, i) = CenterRectOnPointd([0 0 s1, s2], xg(i), yg(i));
    end
    
    % APERTURES:
    %{
    [x,y]         = meshgrid(-1*(apertPlaneSize_px/2):(apertPlaneSize_px/2)-1, -1*(apertPlaneSize_px/2): (apertPlaneSize_px/2)-1);
    %rmin_bg    = 45.6874;% pixels 
    apertureSize_px = 20; %137.7631;% pixels  Not sure what dimension 
    %pa.rstrip     = 11.6268;% this cuts a strip into the fixation disk that has a height the size of the paddle height
    aperturePlane = zeros(apertPlaneSize_px,apertPlaneSize_px).*255;
    aperture      = ones(size(x'));
    aperture      = min(aperture, ((sqrt((x').^2+(y').^2) > apertureSize_px))); % | ((abs(y') > pa.rstrip) & sqrt((x').^2+(y').^2) < pa.rmin_bg)));
    aperturePlane = repmat(aperturePlane,[ 1 1 4 ]);
    aperturePlane = permute(aperturePlane,[ 3 2 1 ]);
    aperturePlane(4,:,:) = shiftdim(255 .* aperture, -1);
    
    glBindTexture(type, textureid(2));
        glTexParameterfv(type, GL.TEXTURE_WRAP_S, GL.REPEAT);
        glTexParameterfv(type, GL.TEXTURE_WRAP_T, GL.REPEAT);
        glTexParameterfv(type, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
        glTexParameterfv(type, GL.TEXTURE_MIN_FILTER, GL.LINEAR); 
        glTexImage2D(type, 0, GL.RGBA, diskTextureWidth*2, diskTextureHeight*2, 0, GL.RGBA, GL.UNSIGNED_BYTE, uint8(aperturePlane));
        glTexEnvi(GL.POINT_SPRITE, GL.COORD_REPLACE, GL.TRUE);
    glBindTexture(type,0)
    
    %}
    shape.disk.texture.id     = textureid;
    %shape.disk.texture.aptrId = textureid(2);
    
end