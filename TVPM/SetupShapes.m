function shape = SetupShapes(ds, pa)

% Some explanation of OpenGL variables:
%
% . listIds is an array of integers (Call List Ids) that correspond to a block of code in SetupShapes responsible for setting up
%   the vertices and shape type to be rendered to the screen. As such the listIds are like identifiers for macros.
%
% . textureId is an integer that corresponds to the plaid texture to be
%   drawn on the plaids. Each plaid has the same texture, so there is only
%   one texture id.
%
% . openglTextureType is a variable global to SetupShapes and should always
%   be set to GL.Texture_2D since we are drawing 2d textures (that just happen to be projected into virtual 3d space).

InitializeMatlabOpenGL(1);
global GL;
Screen('BeginOpenGL', ds.w);
type = GL.TEXTURE_2D;

%% Set up texture and dimensions for the three planes
planeWidth_deg  = ds.hFOV_perPersonAvg/3; % 
planeHeight_deg = ds.vFOV_perPersonAvg;

planeTextureWidth_px  = planeWidth_deg * ds.hor_px_per_deg;
planeTextureHeight_px = planeHeight_deg * ds.ver_px_per_deg;

shape.plane.TextureWidth_px = planeTextureWidth_px;
shape.plane.TextureHeight_px  = planeTextureHeight_px;

xHalfTextureSize = planeTextureWidth_px/2;
yHalfTextureSize = floor(planeTextureHeight_px/2);

[x,y] = meshgrid(-xHalfTextureSize+1:xHalfTextureSize,-yHalfTextureSize+1:yHalfTextureSize);

midPlaneAlpha = ones(size(x));

% All three plaid are painted with the same texture data
% (staticPliadData),
% except midPlane needs a hole in its alpha chanel (alpha chanel = 4th index)
fixationDotSize_px = pa.apertureDia_px;%pa.fixationDiameter; % ATTN!!!!!: convert to deg! WHAT IS GOING ON WITH THIS!
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
    shape.plane.textures(i).bind; % applies the appropriate texture to the geometry set up by listIds.
    glBegin(GL.POLYGON);
    for j = 1:numVertices
        glTexCoord2dv(corners(j,:));
        glVertex3dv(ithPlaneVertices(:,j));
    end
    glEnd;
    shape.plane.textures(i).unbind;
    glEndList();
end

shape.plane.offsets_deg = [-planeWidth_deg, 0, planeWidth_deg];

shape.plane.depths_m    = depths_m;
shape.plane.widths_m    = widths_m;
shape.plane.heights_m   = heights_m;
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
shape.plane.offsets_m   = [-widths_m(1),0,widths_m(3)];
shape.plane.offsets_deg = [-ds.hFOV_perPersonAvg/3, 0, ds.hFOV_perPersonAvg/3];
shape.plane.numPlanes   = numPlanes;
shape.plane.planes      = [shape.plane.near, shape.plane.mid, shape.plane.far];
shape.plane.numVertices = numVertices;
%% Set up texture and dimensions for the plaids and apertures for TVPMS

diskSize_deg       = 2.5;
diskTextureWidth   = ceil(diskSize_deg * ds.hor_px_per_deg); 
diskTextureHeight  = ceil(diskTextureWidth);

halfDiskTexWidth   = diskTextureWidth/2;
halfDiskTexHeight  = diskTextureHeight/2;
[shape.disk.texture.x,shape.disk.texture.y] = meshgrid(-halfDiskTexWidth+1:halfDiskTexWidth,-halfDiskTexHeight+1:halfDiskTexHeight);
shape.disk.texture.width  = diskTextureWidth;
shape.disk.texture.height = diskTextureHeight;

 
numDisksPerPlane = 4;
shape.disk.numDisksPerPlane = numDisksPerPlane;
shape.disk.numDisks = shape.disk.numDisksPerPlane * shape.plane.numPlanes;
shape.disk.diskSize_deg = diskSize_deg;

shape.disk = diskPos(ds, shape.disk, shape.plane);

sf        = pa.sf; % cycles per deg
sf        = (ds.deg_per_px).*sf; %cyc/deg * deg/px = cyc/px
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
glTexEnvfv(GL.TEXTURE_ENV, GL.TEXTURE_ENV_MODE, GL.MODULATE);
glBindTexture(type,0);

shape.disk.texture.id = textureid(1);



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
        glBindTexture(type, textureid(1)); % applies the appropriate texture to the geometry set up by listIds.
        glBegin(GL.POLYGON);
        for j = 1:numVertices
            glTexCoord2dv(corners(j,:));
            glVertex3dv(ithDiskVertices(:,j));
        end
        glEnd;
        glBindTexture(type,0);
        glEndList();
end




Screen('EndOpenGL', ds.w);

%% Full Window Mask for TVPMCD
% Unlike the other shapes the Full Window Mask is drawn with high level
% PsychToolBox functions since the Full Window Mask essentially has no
% depth and can be expressed with only 2 dimensions.
apertureDia_deg = 1;
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
      raisedCos = min(raisedCos, sqrt(xm.^2 + ym.^2) > apertureRadius_px);
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
bag(:,:,4)       = fliplr(alphaMask(:,:,1)');
shape.mask.fullWindowMaskLeftEye = Screen('MakeTexture', ds.w, bag);

%Left Eye
alphaMask(:,:,1) = raisedCos;%raisedCos; % restart Alpha with original raised cosine calc.
frontCols = alphaMask(1:end-round(abs(d_px)/2)+1,:,1); % Take the columns up to end (disparity shift)
alphaMask(:,:,1) = 1; % reset Alpha to completely opaque
alphaMask(round(abs(d_px)/2):end,:,1) = frontCols; % paste to the end of the Alpha chanel
bag(:,:,4) = fliplr(alphaMask(:,:,1)');% circshift(mask(:,:,2)', 50,2);
shape.mask.fullWindowMaskRightEye = Screen('MakeTexture', ds.w, bag);

%% Masks for TVPMSD
Screen('BeginOpenGL', ds.w);
maskOpacity = 175; % for debugging purposes, opacity range = [0,255], where 255 is full opacity, 0 is completely translucent.

% Get plaid positions for TVPMSD, bounded by shape.plane

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
for i = 1:shape.plane.numPlanes % Use visual angle formula to convert from degrees to meters.
    maskWidths_m(i)  = 2 * -maskDepths_m(i) * tand((ds.hFOV_perPersonAvg/shape.plane.numPlanes)/2);
    maskHeights_m(i) = 2 * -maskDepths_m(i) * tand(ds.vFOV_perPersonAvg/2);
end

shape.mask.widths_m = maskWidths_m;
shape.mask.heights_m = maskHeights_m;

for i = 1:shape.plane.numPlanes
    % blackTexData is the image of the mask and initially has 3 RGB channels set to black
    % we then add an alpha channel opaque to blackTexData, opaque is alpha
    % values that are either 0 opacity (the holes) and 255 full opacity.
    % use imagesc(opaque) to see the where the 'holes' appear in the alpha
    % channel
    blackTexData = zeros(maskTexHalfWidth*2, maskTexHalfHeight*2);
    blackTexData = repmat(blackTexData',[ 1 1 3 ]);
    blackTexData = permute(uint8(blackTexData),[ 3 2 1 ]);
    
    x_px = shape.disk.X_px{i};
    y_px = shape.disk.Y_px{i};
    for j = 1:shape.disk.numDisksPerPlane
        opaque = min(opaque, sqrt((x + x_px(j)).^2 + (y - y_px(j)).^2)>apertureRadius_px);
    end
    
    if i == 2
        opaque = min(opaque, sqrt(x.^2 + y.^2) > apertureRadius_px) ; % cut out hole in middle mask for fixation dot
    end
    
    ithMaskVertices = [-maskWidths_m(i)/2 -maskHeights_m(i)/2 maskDepths_m(i);...
        maskWidths_m(i)/2 -maskHeights_m(i)/2 maskDepths_m(i);...
        maskWidths_m(i)/2  maskHeights_m(i)/2 maskDepths_m(i);...
        -maskWidths_m(i)/2  maskHeights_m(i)/2 maskDepths_m(i)]';
    
    blackTexData(4,:,:) = shiftdim(maskOpacity .* opaque', -1);%
    
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


% make phase modulation limits in m for each plane
for ii =1:shape.plane.numPlanes
    shape.cyc_m(ii) =  2*-depths_m(ii)*tand((1/pa.sf)/2);
end

end