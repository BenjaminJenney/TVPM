function diskTexture = DiskTextureold(GL, type, width, height, data)

diskTexture.type   = type;
diskTexture.width  = width;
diskTexture.height = height;

diskTexture.id = glGenTextures(1);
%glEnable(GL.TEXTURE_2D); % Enable 2D texture mapping

glBindTexture(type, diskTexture.id);

glTexParameterfv(type, GL.TEXTURE_WRAP_S, GL.REPEAT);
glTexParameterfv(type, GL.TEXTURE_WRAP_T, GL.REPEAT);
glTexParameterfv(type, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
glTexParameterfv(type, GL.TEXTURE_MIN_FILTER, GL.LINEAR);
glTexImage2D(type, 0, GL.RGB, width, height, 0, GL.RGB, GL.UNSIGNED_BYTE, uint8(data));
glTexEnvi(GL.POINT_SPRITE, GL.COORD_REPLACE, GL.TRUE);
%glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

%glBindTexture(type, 0);

end