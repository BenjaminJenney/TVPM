function diskTexture = DiskTexture(GL, type, width, height, data)

diskTexture.type   = type;
diskTexture.width  = width;
diskTexture.height = height;

diskTexture.id = glGenTextures(2);
%glEnable(GL.TEXTURE_2D); % Enable 2D texture mapping

glBindTexture(type, diskTexture.id(1));

glTexParameterfv(type, GL.TEXTURE_WRAP_S, GL.REPEAT);
glTexParameterfv(type, GL.TEXTURE_WRAP_T, GL.REPEAT);
glTexParameterfv(type, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
glTexParameterfv(type, GL.TEXTURE_MIN_FILTER, GL.LINEAR);
glTexImage2D(type, 0, GL.RGB, width, height, 0, GL.RGB, GL.UNSIGNED_BYTE, uint8(data));
glTexEnvi(GL.POINT_SPRITE, GL.COORD_REPLACE, GL.TRUE);
%glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
glBindTexture(type, 0)
s1 = 512; s2 = 512;
aperture = zeros(s1,s2);
f = abs(cos(1/15.*[-s2/2:s2/2])).^8; f2 = f'*f;
gauss = f2(1:s1,1:s2);
aperture = repmat(aperture,[ 1 1 4 ]);
aperture = permute(aperture,[ 3 2 1 ]);
%keyboard
aperture(4,:,:) = gauss;
glBindTexture(type, diskTexture.id(2));
glTexParameterfv(type, GL.TEXTURE_WRAP_S, GL.REPEAT);
glTexParameterfv(type, GL.TEXTURE_WRAP_T, GL.REPEAT);
glTexParameterfv(type, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
glTexParameterfv(type, GL.TEXTURE_MIN_FILTER, GL.LINEAR);
glTexImage2D(type, 0, GL.RGBA, width, height, 0, GL.RGBA, GL.UNSIGNED_BYTE, uint8(aperture));
glTexEnvi(GL.POINT_SPRITE, GL.COORD_REPLACE, GL.TRUE);
glBindTexture(type, 0)
end