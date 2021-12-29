function shape = SetupShapes(ds, pa)
    InitializeMatlabOpenGL(1);
    global GL;
    Screen('BeginOpenGL', ds.w);
    %Set up texture for the three planes
    planeTextureWidth_px   = (ds.hFOV_perPersonAvg/3) * ds.px_per_deg;
    planeTextureHeight_px  = (ds.vFOV_perPersonAvg) * ds.px_per_deg;

    xHalfTextureSize = planeTextureWidth_px/2;
    yHalfTextureSize = floor(planeTextureHeight_px/2);

    [x,y] = meshgrid(-xHalfTextureSize+1:xHalfTextureSize,-yHalfTextureSize+1:yHalfTextureSize);

    spatialFrequency = 1;% cycles/deg % sf of plaid
    shape.plane.textureData  = staticPlaidData(ds, x, y);

    shape.plane.texture = Texture(GL, GL_TEXTURE_2D, xHalfTextureSize*2, yHalfTextureSize*2, shape.plane.textureData);
    % End texture for the three planes

    % Setup vertices for the three planes
    
    depths_m     = [-.5, -1, -2]; % near, mid, far.
    numPlanes    = length(depths_m);
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
    
    for i = 1:numPlanes

        ithPlaneVertices =[-widths_m(i)/2 -heights_m(i)/2 depths_m(i) ;... 
                            widths_m(i)/2 -heights_m(i)/2 depths_m(i) ;...
                            widths_m(i)/2  heights_m(i)/2 depths_m(i) ;...
                           -widths_m(i)/2  heights_m(i)/2 depths_m(i) ]';

        %ds.initPositions{i} = opticFlow(shape.plane.widths_m(i), shape.plane.heights_m(i), shape.plane.depths_m(i), 50);
        shape.plane.listIds(i) = glGenLists(1);

        glNewList(shape.plane.listIds(i), GL.COMPILE);
            shape.plane.texture.bind;
            glBegin(GL.POLYGON);
            for j = 1:numVertices
                glTexCoord2dv(corners(j,:));
                glVertex3dv(ithPlaneVertices(:,j));
            end
            glEnd;
            shape.plane.texture.unbind;
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
    
    [shape.disk.texture.x,shape.disk.texture.y] = meshgrid(-128+1:128,-128+1:128);
    %u = 1; v = 1;
    %[shape.disk.movieData, shape.disk.numFrames] = moviePlaidData(ds,pa,x,y,u,v);

    shape.disk.texture.width  = 256;
    shape.disk.texture.height = 256;
    
    %shape.disk.coords   = randn(3,100);
    shape.disk.size_px  = 23; %not totally sure what units this is. glPointSize draws a square with equal sides in pixels, supposedly.
    shape.disk.size_deg = shape.disk.size_px * ds.deg_per_px; %1.4888
    
%     for i = 1:shape.disk.numFrames
%         shape.disk.textureArray(i) = PointTexture(GL, GL_TEXTURE_2D, pointTextureWidth, pointTextureHeight, shape.disk.movieData{i}); 
%     end
    shape.disk.numDisksPerPlane = 60;
    shape.disk.numDisks = shape.disk.numDisksPerPlane * shape.plane.numPlanes;
    
    shape.disk = diskPos(ds, shape.disk, shape.plane);
    
    
    Screen('EndOpenGL', ds.w);
end