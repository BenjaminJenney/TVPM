function shape = SetupShapes(ds, pa)
InitializeMatlabOpenGL(1);
global GL;
Screen('BeginOpenGL', ds.w);
type = GL.TEXTURE_2D;
%% Set up texture and dimensions for the three planes

planeTextureWidth_px   = (ds.hFOV_perPersonAvg/3) * ds.hor_px_per_deg;
planeTextureHeight_px  = (ds.vFOV_perPersonAvg) * ds.ver_px_per_deg;
shape.plane.TextureWidth_px = planeTextureWidth_px;
shape.plane.TextureHeight_px  = planeTextureHeight_px;
xHalfTextureSize = planeTextureWidth_px/2;
yHalfTextureSize = floor(planeTextureHeight_px/2);

[x,y] = meshgrid(-xHalfTextureSize+1:xHalfTextureSize,-yHalfTextureSize+1:yHalfTextureSize);

midPlaneAlpha = ones(size(x));

% All three plaid are painted with the same texture data
% (staticPliadData),
% except midPlane needs a hole in its alpha chanel (alpha chanel = 4th index)
fixationDotSize_px = 15; % ATTN: convert to deg!
shape.plane.textureData  = staticPlaidData(ds, pa.sf, x, y);
midPlaneAlpha   = min(midPlaneAlpha, sqrt(x.^2 + y.^2) > fixationDotSize_px);%Cut out hole for fixation dot
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
    heights_m(i) = (2 * -depths_m(i) * tand(ds.vFOV_perPersonAvg/2));
end
%keyboard;
numVertices = 4;
corners = [0 0;
    1 0;
    1 1;
    0 1];

shape.plane.listIds  = zeros(1, numPlanes);
shape.plane.vertices = {};
for i = 1:numPlanes
    
    ithPlaneVertices =[-widths_m(i)/2 -heights_m(i)/2 depths_m(i) ;...
        widths_m(i)/2 -heights_m(i)/2 depths_m(i) ;...
        widths_m(i)/2  heights_m(i)/2 depths_m(i) ;...
        -widths_m(i)/2  heights_m(i)/2 depths_m(i) ]';
    shape.plane.vertices{i} = ithPlaneVertices;
    shape.plane.listIds(i)  = glGenLists(1);
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

shape.plane.depths_m      = depths_m;
shape.plane.widths_m      = widths_m;
shape.plane.heights_m     = heights_m;
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
diskSize_deg       = 2.5;
diskTextureWidth   = diskSize_deg * ds.ver_px_per_deg; 
diskTextureHeight  = diskTextureWidth;

halfDiskTexWidth   = diskTextureWidth/2;
halfDiskTexHeight  = diskTextureHeight/2;
[shape.disk.texture.x,shape.disk.texture.y] = meshgrid(-halfDiskTexWidth+1:halfDiskTexWidth,-halfDiskTexHeight+1:halfDiskTexHeight);
shape.disk.texture.width  = diskTextureWidth;
shape.disk.texture.height = diskTextureHeight;
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

shape.disk.texture.id = textureid(1);

numDisksPerPlane = 10;

for i = 1:numPlanes
    plaidWidths_m(i)  = 2 * -depths_m(i) * tand(diskSize_deg/2);
    plaidHeights_m(i) = (2 * -depths_m(i) * tand(diskSize_deg/2));
end

shape.disk.listIds  = zeros(1, numPlanes);
shape.disk.vertices = {};

for i = 1:numPlanes
    ithDiskVertices =[-plaidWidths_m(i)/2 -plaidHeights_m(i)/2 depths_m(i) ;...
                        plaidWidths_m(i)/2 -plaidHeights_m(i)/2 depths_m(i) ;...
                        plaidWidths_m(i)/2  plaidHeights_m(i)/2 depths_m(i) ;...
                        -plaidWidths_m(i)/2  plaidHeights_m(i)/2 depths_m(i) ]';
        shape.disk.vertices{i} = ithDiskVertices;
        shape.disk.listIds(i)  = glGenLists(1);
        glNewList(shape.disk.listIds(i), GL.COMPILE);
        glBindTexture(type, textureid(1));
        glBegin(GL.POLYGON);
        for j = 1:numVertices
            glTexCoord2dv(corners(j,:));
            glVertex3dv(ithDiskVertices(:,j));
        end
        glEnd;
        glBindTexture(type,0);
        glEndList();
end


shape.disk.numDisksPerPlane = numDisksPerPlane;
shape.disk.numDisks = shape.disk.numDisksPerPlane * shape.plane.numPlanes;
shape.disk.diskSize_deg = diskSize_deg;
shape.disk = diskPos(ds, shape.disk, shape.plane);


Screen('EndOpenGL', ds.w);

%% Full Window Mask
% Unlike the other shapes the Full Window Mask is drawn with high level
% PsychToolBox functions since the Full Window Mask essentially has no
% depth and can be expressed with only 2 dimensions.
apertureDia_deg = diskSize_deg/1.5;
apertureRadius_px  = (apertureDia_deg/2)*ds.hor_px_per_deg;%sqrt((apertureDia_deg^2 + apertureDia_deg^2)*ds.hor_px_per_deg);

%{ Force disparity of the Full Window Mask to match the disparity of the middle plane %}
b = .067; % inter-pupillary distance in meters TODO: Figure out if there is a way to get IP from the headset instead of hardcoding.
Z = abs(shape.plane.depths_m(2)); % depth of middle plane.
d_deg = -2.*( atand((b./2)./(2.*Z)) - atand((-b./2)./(2.*Z)) ); % Should fix disparity of the mask.
d_px  = abs(d_deg*ds.hor_px_per_deg);
d_px = 4*d_px;  %ATTN: THIS IS WRONG AND NEEDS TO BE FIXED
screenXpixels = ds.screenRenderWidthMonocular_px;
screenYpixels = ds.screenRenderHeightMonocular_px;
black = BlackIndex(0);

%{ Make a raised Cosine aperture for the "alpha" channel of the mask %}
%     gaussDim = pa.apertureDia_px; %100 px? TODO: DOES NOT SEEM TO CHANGE SIZE.
%     gaussSigma = pa.apertureDia_px; %gaussDim /4; % csb: controls how big the apertures are
s1 = screenXpixels;
s2 = screenYpixels;
[xm, ym] = meshgrid(-(s2/2)+1:s2/2, -(s1/2)+1:s1/2); %%ben flipped this, originally it was meshgrid(-(s2/2)+1:s2/2, -(s1/2)+1:s1/2);
raisedCos = ones(size(xm));

      for i = 1:shape.plane.numPlanes
        x_px = shape.disk.xpos_deg{i} .* ds.hor_px_per_deg;
        y_px = shape.disk.ypos_deg{i} .* ds.hor_px_per_deg;
        for j = 1:shape.disk.numDisksPerPlane
            raisedCos = min(raisedCos, sqrt((xm + y_px(j)).^2 + (ym + x_px(j)).^2) > apertureRadius_px);
        end
      end
% Create hole for fixation dot
%raisedCos = min(raisedCos, sqrt((xm + 0).^2 + (ym + 0).^2) > pa.apertureDia_px);
%       figure; imagesc(raisedCos)
%x_px(j)
% to make a grid,uncomment:
%     f = abs(cos(1/15.*[-s2/2:s2/2])).^8; f2 = f'*f
%      raisedCos = f2(1:s1,1:s2);
%     fixationHole_px = 15 ;
%     raisedCos = 1-raisedCos;



%{ Create texture image data for Full Window Mask %}
bag = ones(s2, s1) .* black; % Ben flipped this originally it was:  ones(screenYpixels, screenXpixels) .* black;
bag = repmat(bag,[ 1 1 3 ]);

% REMEMBER: d_px = d_px * 4; THIS IS A HACK
%Force Disparity Right Eye
alphaMask        = ones(s1, s2, 1); % init alpha channel  TPB Standard 0:completely transluscent -- 255:completely opaque.
alphaMask(:,:,1) = raisedCos;%raisedCos; % Original Raised cosine calculation

alphaMask(:,:,1) = round(alphaMask(:,:,1),3,'significant');
%bagBeforeDisparity(:,:,4) = flipud(rot90(alphaMask(:,:,1),1));
endCols          = alphaMask(round(abs(d_px)/2)+1:end,:,1);
alphaMask(:,:,1) = 1;
alphaMask(1:end - round(abs(d_px)/2),:,1) = endCols;
bag(:,:,4)       = fliplr(flipud(alphaMask(:,:,1)'));
shape.mask.fullWindowMaskLeftEye = Screen('MakeTexture', ds.w, bag);

%Left Eye
alphaMask(:,:,1) = raisedCos;%raisedCos; % restart Alpha with original raised cosine calc.
frontCols = alphaMask(1:end-round(abs(d_px)/2)+1,:,1); % Take the columns up to end (disparity shift)
alphaMask(:,:,1) = 1; % reset Alpha to completely opaque
alphaMask(round(abs(d_px)/2):end,:,1) = frontCols; % paste to the end of the Alpha chanel
bag(:,:,4) = fliplr(flipud(alphaMask(:,:,1)'));% circshift(mask(:,:,2)', 50,2);
shape.mask.fullWindowMaskRightEye = Screen('MakeTexture', ds.w, bag);

%% Masks for TVPM Full
Screen('BeginOpenGL', ds.w);

shape.disk = diskPos(ds, shape.disk, shape.plane); % Get plaid positions for TVPM-Full, bounded by shape.plane

corners = [0 0;
    1 0;
    1 1;
    0 1];

maskTexHalfWidth  = floor(shape.plane.TextureWidth_px/2);
maskTexHalfHeight = floor(shape.plane.TextureHeight_px/2);


[x,y]  = meshgrid(-maskTexHalfWidth+1:maskTexHalfWidth, -maskTexHalfHeight+1:maskTexHalfHeight);
opaque = ones(size(x));
maskWidths_m  = zeros(1, shape.plane.numPlanes);
maskHeights_m = zeros(1, shape.plane.numPlanes);
maskDepths_m  = shape.plane.depths_m + .001;
for i = 1:shape.plane.numPlanes
    maskWidths_m(i)  = 2 * -maskDepths_m(i) * tand((ds.hFOV_perPersonAvg/shape.plane.numPlanes)/2);
    maskHeights_m(i) = 2 * -maskDepths_m(i) * tand(ds.vFOV_perPersonAvg/2);
end

shape.mask.widths_m = maskWidths_m;
shape.mask.heights_m = maskHeights_m;

for i = 1:shape.plane.numPlanes
    
    blackTexData = zeros(maskTexHalfWidth*2, maskTexHalfHeight*2);
    blackTexData = repmat(blackTexData',[ 1 1 3 ]);
    blackTexData = permute(uint8(blackTexData),[ 3 2 1 ]);
    
    x_px = shape.disk.X_px{i};
    y_px = shape.disk.Y_px{i};
    for j = 1:shape.disk.numDisksPerPlane
        opaque = min(opaque, sqrt((x + x_px(j)).^2 + (y + y_px(j)).^2)>apertureRadius_px);
    end
    
    %planeAlphas{i} = opaque;
    ithMaskVertices = [-maskWidths_m(i)/2 -maskHeights_m(i)/2 maskDepths_m(i);...
        maskWidths_m(i)/2 -maskHeights_m(i)/2 maskDepths_m(i);...
        maskWidths_m(i)/2  maskHeights_m(i)/2 maskDepths_m(i);...
        -maskWidths_m(i)/2  maskHeights_m(i)/2 maskDepths_m(i)]';
    
    blackTexData(4,:,:) = shiftdim(255 .* opaque', -1);%
    
    shape.mask.texture(i) = Texture(GL, type, 2*maskTexHalfWidth, 2*maskTexHalfHeight, uint8(blackTexData));
    shape.mask.listIds(i) = glGenLists(1);
    glNewList(shape.mask.listIds(i), GL.COMPILE);
    shape.mask.texture(i).bind
    glBegin(GL.POLYGON);
    for j = 1:shape.plane.numVertices
        glTexCoord2dv(corners(j,:));
        glVertex3dv(ithMaskVertices(:,j));
    end
    glEnd;
    shape.mask.texture(i).unbind;
    glEndList();
    opaque = ones(size(x));
end
Screen('EndOpenGL', ds.w);

%ds.masktex = Screen('MakeTexture', ds.w, mask);

%black = ones(screenYpixels, screenXpixels) .* black;
%black(:, : 4) = mask;

%bag = permute(bag,[ 3 2 1 ]);

%bag(628:668,:,4) = -1;
% Make a grey texture to cover the full window
%shape.mask.fullWindowMask = Screen('MakeTexture', ds.w, bag);

%     % Make coordinates in which to draw the apertures into our full screen mask
%     [xg, yg] = meshgrid(-gaussDim:gaussDim, -gaussDim:gaussDim); % csb: and this controls how many apertures there are
%     spacing = gaussDim * 2; % csb: i think this controls how far apart the apertures ares
%     xg = xg .* spacing + screenXpixels / 2;
%     yg = yg .* spacing + screenYpixels / 2; % hole y position%
%     xg = reshape(xg, 1, numel(xg));
%     yg = reshape(yg, 1, numel(yg));
%
%     shape.mask.gauss = gauss;
%
%     % Make the destination rectangles for the gaussian apertures
%     shape.dstRects = nan(4, numel(xg));
%     for i = 1:numel(xg)
%     shape.dstRects(:, i) = CenterRectOnPointd([0 0 s1, s2], xg(i), yg(i));
%     end

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
        glTex Image2D(type, 0, GL.RGBA, diskTextureWidth*2, diskTextureHeight*2, 0, GL.RGBA, GL.UNSIGNED_BYTE, uint8(aperturePlane));
        glTexEnvi(GL.POINT_SPRITE, GL.COORD_REPLACE, GL.TRUE);
    glBindTexture(type,0)
    
%}
%shape.disk.texture.id     = textureid;
%shape.disk.texture.aptrId = textureid(2);

%     figure; imagesc(bagBeforeDisparity(:,:,4)); axis equal;
%      planeAlpha = [planeAlphas{1} planeAlphas{2} planeAlphas{3}]';
%      figure; imagesc(planeAlpha); axis equal;


% make phase modulation limits in m for each plane
for ii =1:shape.plane.numPlanes
    shape.cyc_m(ii) =  2*-depths_m(ii)*tand((1/pa.sf)/2);
end

end