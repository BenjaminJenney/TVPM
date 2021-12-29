classdef PointTexture
   properties
       id;
       width;
       height;
       type;
       data;
       GL;
   end
   methods
       function self = PointTexture(GL, type, width, height, data)
%            if nargin < 5
%                error("expected 5 arguments: texture_id, width, height, texture_data");
%                return
%            end
           
           self.type   = type;
           self.id     = glGenTextures(type);
           self.width  = width;
           self.height = height;
           self.data   = data;
           self.GL     = GL;
           
           glBindTexture(type, self.id);
           
           glTexParameterfv(type, GL.TEXTURE_WRAP_S, GL.REPEAT);
           glTexParameterfv(type, GL.TEXTURE_WRAP_T, GL.REPEAT);
           glTexParameterfv(type, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
           glTexParameterfv(type, GL.TEXTURE_MIN_FILTER, GL.LINEAR);
           glTexImage2D(type, 0, GL.RGB, width, height, 0, GL.RGB, GL.UNSIGNED_BYTE, data);
           glTexEnvi(GL.POINT_SPRITE, GL.COORD_REPLACE, GL.TRUE);
           glBindTexture(type, 0);
       end
       
       function bind(self)
           glBindTexture(self.type, self.id);
       end
       
       function unbind(self)
           glBindTexture(self.type, 0);
       end
       function draw(self, coords, size)
           [rows, cols] = size(coords);
           glEnable(self.GL.POINT_SPRITE);
           glEnable(self.GL.POINT_SMOOTH);
           glPointSize(20);
           self.bind;
           
           glBegin(self.GL.POINTS);
               for i = 1:cols 
                   glVertex3dv(coords(:,i));
               end
           glEnd;
           
           self.unbind;
           glDisable(self.GL.POINT_SPRITE);
       end
   end
end