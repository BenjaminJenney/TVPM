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

    
    shape.plane.textureData  = staticPlaidData(ds, x, y);

    shape.plane.texture = Texture(GL, GL_TEXTURE_2D, xHalfTextureSize*2, yHalfTextureSize*2, shape.plane.textureData);
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
%         shape.plane.listIds(i) = glGenLists(1);
% 
%         glNewList(shape.plane.listIds(i), GL.COMPILE);
%             shape.plane.texture.bind;
%             glBegin(GL.POLYGON);
%             for j = 1:numVertices
%                 glTexCoord2dv(corners(j,:));
%                 glVertex3dv(ithPlaneVertices(:,j));
%             end
%             glEnd;
%             shape.plane.texture.unbind;
%         glEndList(); 
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
    
    sf        = 2; %.1; % cycles per deg
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
    %% Masks
    
    
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
    Screen('EndOpenGL', ds.w);
end