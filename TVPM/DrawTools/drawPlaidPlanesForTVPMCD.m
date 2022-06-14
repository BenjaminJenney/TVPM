function drawPlaidPlanesForTVPMCD(plaidPlaneWidths_m, listIds, xDisplacement, zDisplacement, middlePlaneDepth)
    % draws the three plaid planes set up in SetupShape.
    %
    % INPUT:
    %   * plaidPlaneWidths_m is the widths of the three planes in meters.
    %   * listIds references code that generates the geometry for the
    %     three planes and wraps (binds) the appropriate texture to them.
    %   * xDisplacement, zDisplacement: Moves the planes incrementally on
    %     each frame. Simulates Optic Flow with changing disparity using
    %     Opengl default.
    
    glPushMatrix;
    glTranslatef(-plaidPlaneWidths_m(1), 0.0, 0.0);
    gluLookAt(-xDisplacement,0,-zDisplacement,0,0,middlePlaneDepth,0,1,0)
    glCallList(listIds(1));
    glPopMatrix;

    glPushMatrix;
    glTranslatef(0.0, 0.0, 0.0);
    gluLookAt(-xDisplacement,0,-zDisplacement,0,0,middlePlaneDepth,0,1,0)
    glCallList(listIds(2));
    glPopMatrix;

    glPushMatrix;
    glTranslatef(plaidPlaneWidths_m(3), 0.0, 0.0);
    gluLookAt(-xDisplacement,0,-zDisplacement,0,0,middlePlaneDepth,0,1,0)
    glCallList(listIds(3));
    glPopMatrix;

    
    
    
        
end