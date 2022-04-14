classdef Texture
   properties
       id uint32;
       width int32;
       height int32;
       type;
       data;
   end
   methods
       function self = Texture(GL, type, width, height, data)
           if nargin < 4
               error("expected 4 arguments: texture_id, width, height, texture_data");
               return
           end
           
           self.type   = type;
           self.id     = glGenTextures(type);
           self.width  = width;
           self.height = height;
           self.data   = data;
           
           glBindTexture(type, self.id);
           
           glTexParameterfv(type, GL.TEXTURE_WRAP_S, GL.REPEAT);
           glTexParameterfv(type, GL.TEXTURE_WRAP_T, GL.REPEAT);
           glTexParameterfv(type, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
           glTexParameterfv(type, GL.TEXTURE_MIN_FILTER, GL.LINEAR);
           glTexImage2D(type, 0, GL.RGBA, width, height, 0, GL.RGBA, GL.UNSIGNED_BYTE, data);
           glTexEnvfv(GL.TEXTURE_ENV, GL.TEXTURE_ENV_MODE, GL.MODULATE);
           
           
           glBindTexture(type, 0);
       end
       
       function bind(self)
           glBindTexture(self.type, self.id);
       end
       
       function unbind(self)
           glBindTexture(self.type, 0);
       end
   end
end