%this script contains code for the 3D mask, which we currently are not
%using.

% %% The Holey Bag - swiss cheese plane i.e. a black surround texture with holes in it
% % m_per_px = 0.0002645833; %m/px
% % px_per_m = 1/m_per_px; %px/m
% % screenX_m   = m_per_px*ds.windowRect(3); %0.3556 meters
% % screenY_m   = m_per_px*ds.windowRect(4); %0.4233 meters
% 
% screenX_deg = ds.hFOV_perPersonAvg;
% screenY_deg = ds.vFOV_perPersonAvg;
% px_per_deg  = screenRenderWidthMonocular_px/screenX_deg; %~15.45 this seems pretty good check out: https://www.roadtovr.com/understanding-pixel-density-retinal-resolution-and-why-its-important-for-vr-and-ar-headsets/
% deg_per_px  = 1/px_per_deg;
% focalLength = 1.2; % meters % This is according to Bas
% 
% halfTextureSize = screenRenderWidthMonocular_px/2;
% [x,y] = meshgrid(-halfTextureSize+1:halfTextureSize,-halfTextureSize+1:halfTextureSize);
% 
% blackTexData = ones(screenRenderWidthMonocular_px, screenRenderWidthMonocular_px);
% blackTexData = repmat(blackTexData,[ 1 1 3 ]);
% blackTexData = permute(uint8(blackTexData),[ 3 2 1 ]);
% 
% xoffset = -400:1:400; % hole x Position
% yoffset = randi([-400, 400], 1, length(xoffset)); % hole y position
% 
% pa.rmax_bg = 2;  % px
% 
% % this code pokes out the transparent aperture
% opaque = ones(size(x'));
% 
% for i = 1:length(xoffset)
%     opaque = min(opaque, ((sqrt((x'+xoffset(i)).^2+(y'+yoffset(i)).^2) > pa.rmax_bg)));
% end
% blackTexData(4,:,:) = shiftdim(255 .* opaque, -1); 
% 
% ds.maskTexture = Texture(GL, GL_TEXTURE_2D, halfTextureSize*2, halfTextureSize*2, blackTexData);
% 
% corners = [0 0;
%     1 0;
%     1 1;
%     0 1];
% 
% z = -.5; % meters
% 
% xMask_m = 2*z*tand(ds.hFOV_perPersonAvg/2); 
% 
% v=[-xMask_m -xMask_m  z ;... 
%     xMask_m -xMask_m  z ;...
%     xMask_m  xMask_m  z ;...
%    -xMask_m  xMask_m z ]';
% 
% ds.mask = glGenLists(1);
% glNewList(ds.mask, GL.COMPILE);
% 
% glBegin(GL.POLYGON);
% for i = 1:4
%     glTexCoord2dv(corners(i,:));
%     glVertex3dv(v(:,i));
% end
% glEnd;
% glEndList();