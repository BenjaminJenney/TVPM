function drawPlaidForTVPMSD(openglTextureType, textureId, listIds, numPlanes, numPlaidsPerPlane, posX, posY)
% drawPlaidForTVPMSD draws the plaids at (x,y) coordinates dictated by posX, posY. 
% The (x,y) positions are bounded by each plane and are projected at the respective depths of each plane (shape.plane.depths_m).
% There are always 3 planes, and the number of plaids drawn is controlled by shape.disk.numDisksPerPlane.
% 
% Some explanation of OpenGL variables: 
% . listIds is an array of integers (Call List Ids) that correspond to a block of code in SetupShapes responsible for setting up
% the vertices and shape type to be rendered to the screen. As such the listIds are like identifiers for macros.
%
% . textureId is an integer that corresponds to the plaid texture to be
% drawn on the plaids. Each plaid has the same texture, so there is only
% one texture id.
%
% . openglTextureType is a variable global to runTVPMExperiment and should always
% be set to GL_Texture_2D since we are drawing 2d textures (that just happen to be projected into virtual 3d space).

    
    %glBindTexture(openglTextureType, textureId); 
    
    for i = 1:numPlanes
        for j = 1:numPlaidsPerPlane
            glPushMatrix; % Since we are translating the plaids in the virtual world, we are multiplying against the ModelView matrix. We want to keep the ModelView matrix isolated as to not be used in proceeding matrix multiplications.
            glTranslatef(posX{i}(j), posY{i}(j), 0.0) % The plaids are already projected at the appropriate depths (shape.plane.depths_m). We just have to move them from the origin (center of FOV) to the correct position in the virtual world. 
            glCallList(listIds(i)) % Calls a code block in SetupShapes that sets up the plaid geometry for each plaid.
            glPopMatrix; % We are done with the current translation and want to go back to the last loaded ModelView.
        end
    end
    
    %glBindTexture(openglTextureType, 0); % Tells opengl we are finished applying the plaid textures for TVPMSD on this frame

end
                      