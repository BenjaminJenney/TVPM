function diskTexture = DiskTexture(GL, type, width, height, data)

diskTexture.type   = type;
diskTexture.width  = width;
diskTexture.height = height;


%glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
glBindTexture(type, 0)

%keyboard


end