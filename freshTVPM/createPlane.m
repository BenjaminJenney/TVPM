function drawPlaneList = createPlane(GL)
    %h_cm = 193.59/10; % height of my dell g5 5587 monitor
    %d_cm = 78; % distance from my monitor in centimeters

    %v_px = 1080; % vertical resolution of my dell g5 5587 monitor

    %s_px = 1920; % width of stimulus: the three planes


    %deg_per_px = rad2deg(atan2(.5*h_cm, d_cm))/(.5*v_px)

    hfov = 80;
    vfov = 90;

    type = GL.TEXTURE_2D;
    
    renderWidth_px = 1344*2;
    renderHeight_px = 1600;
    planeWidth_px  = renderWidth_px/3;
    planeHeight_px = renderHeight_px;
    halfPlaneWidth_px = planeWidth_px/2;
    halfPlaneHeight_px = planeHeight_px/2;

    [x,y] = meshgrid(-halfPlaneWidth_px+1:halfPlaneWidth_px,-halfPlaneHeight_px+1:halfPlaneHeight_px);
    
    textureData = zeros(size(x));
    textureData = repmat(textureData,[ 1 1 3 ]);
    textureData = permute(uint8(textureData),[ 3 2 1 ]);
    
    textureid = glGenTextures(1);
    glBindTexture(type, textureid(1));
    glTexParameterfv(type, GL.TEXTURE_WRAP_S, GL.REPEAT);
    glTexParameterfv(type, GL.TEXTURE_WRAP_T, GL.REPEAT);
    glTexParameterfv(type, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
    glTexParameterfv(type, GL.TEXTURE_MIN_FILTER, GL.LINEAR);
    glTexImage2D(type, 0, GL.RGB, planeWidth_px, planeHeight_px, 0, GL.RGB, GL.UNSIGNED_BYTE, textureData);
    glTexEnvfv(GL.TEXTURE_ENV, GL.TEXTURE_ENV_MODE, GL.MODULATE);
    glBindTexture(type,0);

    planeDepth_m = -1;
    planeWidth_m  = 2 * -planeDepth_m * tand(hfov/2);
    planeHeight_m = 2 * -planeDepth_m * tand(vfov/2);
    
    corners = [0 0;
    1 0;
    1 1;
    0 1];
    
    planeVertices =[-planeWidth_m/2 -planeHeight_m/2 planeDepth_m;...
                    planeWidth_m/2 -planeHeight_m/2 planeDepth_m;...
                    planeWidth_m/2  planeHeight_m/2 planeDepth_m;...
                    -planeWidth_m/2 planeHeight_m/2 planeDepth_m]';

    drawPlaneList = glGenLists(1);

    glNewList(drawPlaneList, GL.COMPILE);

    glBegin(GL.POLYGON);
    for i = 1:4
        glVertex3dv(planeVertices(:,i));
    end
    glEnd;
    glEndList();
    
end


