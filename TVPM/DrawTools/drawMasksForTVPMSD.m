function drawMasksForTVPMSD(maskWidths_m, maskListIds)
    % Draw the three masks for TVPMSD. 
    % INPUTS: 
    %   * maskWidths_m, an array of the widths of the masks in meters.
    %   * maskListIds, an array of integer ids that reference code blocks
    %       that set up the geometry of the masks and apply the textures to
    %       that geometry. These ids are like the names of macros that are
    %       activated by glCallList.
    %
    % The masks are projected .001 meters in front of the plaids (and as such are marginally different in size from the three planes drawn in TVPMCD),
    % so there depths are shape.plane.depths_m - .001. See SetupShapes for reference.
    %
    % The three masks, like the planes, are projected in the virtual world
    % at different depths, yet each appears as the same size on the retina
    % as long as your head is positioned at the origin of the VR Sensors.
    % Remember there is an origin exogenous to the virtual world dictated
    % by the VR sensors in the physical room that they exist. One only has
    % to move their head back and forth left and right while the headset is on to see this
    % phenomena in action! I emphasis this because ones head has to be at a
    % precise position to see that each of the three masks are in fact the same size on the
    % retina.
    
    glPushMatrix;
    glTranslatef(-maskWidths_m(1), 0.0, 0.0);
    glCallList(maskListIds(1));
    glPopMatrix;

    glPushMatrix;
    glTranslatef(0.0, 0.0, 0.0);
    glCallList(maskListIds(2));
    glPopMatrix;

    glPushMatrix;
    glTranslatef(maskWidths_m(3), 0.0, 0.0);
    glCallList(maskListIds(3));
    glPopMatrix;
end