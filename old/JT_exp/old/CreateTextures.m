function [ds,pa] = CreateTexturesfor(ds,pa)

% OLD

%{
InitializeMatlabOpenGL(1);

Screen('BeginOpenGL', ds.w); % Setup the OpenGL rendering context
glEnable(GL.TEXTURE_2D); % Enable 2D texture mapping

mag_filter = GL_LINEAR; %GL_NEAREST; %_LINEAR;
min_filter = GL_LINEAR; %GL_NEAREST; %GL_LINEAR_MIPMAP_NEAREST; %_LINEAR

% Initialize textures
ds.wall_texid = glGenTextures(1); % this will be the surround texture with the fixation disk and fixation lines embedded
largepaddle_texid = glGenTextures(1); % this will be the texture for the large paddle faces
smallpaddle_texid = glGenTextures(1); % this will be the texture for the small paddle faces
%ds.floor_texid = glGenTextures(1); % this will be the floor texture
%ds.ceiling_texid = glGenTextures(1); % this will be the ceiling texture
%ds.roomwall_texid = glGenTextures(1); % this will be the wall texture

%% Suuround Texture - the fixation plane surround with fixation disk and fixation lines embedded


halfTextureSize = 1024;%ds.xc;  % needs to be in pixels for meshgrid

[x,y] = meshgrid(-halfTextureSize+1:halfTextureSize,-halfTextureSize+1:halfTextureSize);

noysSlope = 1.0; %1.5;
noys = 255.*oneoverf(noysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
noys=repmat(noys,[ 1 1 3 ]);
noys=permute(uint8(noys),[ 3 2 1 ]);

xoffset = 0;
yoffset = 0;
pa.rmin_bg = 45.6874;% pixels 
pa.rmax_bg = 350.7631;% pixels  
pa.rstrip = 11.6268;% this cuts a strip into the fixation disk that has a height the size of the paddle height

%{
% this code pokes out the transparent aperture
opaque = ones(size(x'));
 
for i = 1:length(xoffset)
    opaque = min(opaque, ((sqrt((x'+xoffset(i)).^2+(y'+yoffset(i)).^2) > pa.rmax_bg)  | ((abs(y'+yoffset(i)) > pa.rstrip) & sqrt((x'+xoffset(i)).^2+(y'+yoffset(i)).^2) < pa.rmin_bg)));
end
noys(4,:,:) = shiftdim(255 .* opaque, -1); 
%}

glBindTexture(GL_TEXTURE_2D, ds.wall_texid);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag_filter);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min_filter);
%glTexParameterfv(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
%glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, halfTextureSize*2, halfTextureSize*2, 0, GL_RGBA, GL_UNSIGNED_BYTE, noys);
glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

corners = [0 0;
    1 0;
    1 1;
    0 1];

canvas = max(ds.halfHeight, ds.halfWidth); % pick the larger of the two dimensions, since they are not equal

v=[-canvas -canvas 0 ;... 
    canvas -canvas 0 ;...
    canvas canvas 0 ;...
    -canvas canvas 0]';

ds.surroundTexture = glGenLists(1);
glNewList(ds.surroundTexture, GL.COMPILE);

%{
glBegin(GL.POLYGON);
for i = 1:4
    glTexCoord2dv(corners(i,:));
    glVertex3dv(v(:,i));
end
glEnd;
%}


% Add the fixation lines to the list as well using glLines instead of PTB's DrawLines - more efficient coding


glLineWidth(2);

glBegin(GL.LINES);
glColor3f(0, 0, 0); % black diagonal line left
glVertex3f(-pa.fixationHalfSquare-pa.fixationLineLength, -pa.fixationHalfSquare-pa.fixationLineLength, 0);
glVertex3f(-pa.fixationHalfSquare, -pa.fixationHalfSquare, 0);

glColor3f(1, 0, 0); % red diagonal line left
glVertex3f(-pa.fixationHalfSquare-pa.fixationLineLength, pa.fixationHalfSquare+pa.fixationLineLength, 0);
glVertex3f(-pa.fixationHalfSquare, pa.fixationHalfSquare, 0);

glColor3f(1, 0, 0); % red diagonal line right
glVertex3f(pa.fixationHalfSquare+pa.fixationLineLength, -pa.fixationHalfSquare-pa.fixationLineLength, 0);
glVertex3f(pa.fixationHalfSquare, -pa.fixationHalfSquare, 0);

glColor3f(0, 0, 0); % black diagonal line right
glVertex3f(pa.fixationHalfSquare+pa.fixationLineLength, pa.fixationHalfSquare+pa.fixationLineLength, 0);
glVertex3f(pa.fixationHalfSquare, pa.fixationHalfSquare, 0);
glEnd();

glEndList(); % 1/f noise texture is complete



% Close the OpenGL rendering context
Screen('EndOpenGL', ds.w);

%}

end