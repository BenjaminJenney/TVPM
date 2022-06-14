function drawPlaidsForTVPMSD(listIds, numPlanes, numPlaidsPerPlane, posX, posY)
% drawPlaidForTVPMSD draws the plaids at (x,y) coordinates dictated by posX, posY. 
% INPUTS: 
%   * listIds, listIds is an array of integers (Call List Ids) that correspond to a block of code in SetupShapes responsible for setting up
%     the vertices, shape type and texture to be rendered to the screen. As such the listIds are like identifiers for rendering macros.
%     The (x,y) positions are bounded by each plane and are projected at the respective depths of each plane (shape.plane.depths_m).
%     There are always 3 planes, and the number of plaids drawn is controlled by shape.disk.numDisksPerPlane defined in SetupShapes.m.
%   * numPlanes is the number of planes to be presented. In TVPMSD we are sampling plaid positions from the 3 planes which can be seen in
%     TVPMCD.
%   * numPlaidsPerPlane corresponds with shape.disk.numDisksPerPlane defined in SetupShapes. 
%     Controls for how many Plaids in total are drawn in the VR world. 
%   * posX, posY: positions of the plaids based on optic flow, bounded by the three planes.
    
    for i = 1:numPlanes
        for j = 1:numPlaidsPerPlane
            glPushMatrix; % Since we are translating the plaids in the virtual world, we are multiplying against the ModelView matrix. We want to keep the ModelView matrix isolated as to not be used in proceeding matrix multiplications.
            glTranslatef(posX{i}(j), posY{i}(j), 0.0) % The plaids are already projected at the appropriate depths (shape.plane.depths_m). We just have to move them from the origin (center of FOV) to the correct position in the virtual world. 
            glCallList(listIds(i)) % Calls a code block in SetupShapes that sets up the plaid geometry for each plaid and wraps the appropriate texture to that geometry.
            glPopMatrix; % We are done with the current translation and want to go back to the default ModelView loaded in Opengl.
        end
    end

end
                      