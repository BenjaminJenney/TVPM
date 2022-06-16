function disk = diskPos(ds, disk, plane)

    hFOV   = ds.hFOV_perPersonAvg;
    vFOV   = ds.vFOV_perPersonAvg;
    T      = plane.offsets_deg;
    n      = plane.numPlanes;
    ndp    = disk.numDisksPerPlane;
    Z      = plane.depths_m;
    % Width of plane = 29^o, if we focus at the center of the screen (0,0),
    % the line x=0 intersects our gaze forming a line of symmetry to the
    % image plane, thus the x-dimension of the middle plane stimulus is [-14.5^o,14.5^o]
    a  = -hFOV/(2*n); 
    b  = hFOV/(2*n); % x range = [a,b]
    a_ = -vFOV/2;
    b_ = vFOV/2; % y range = [a_, b_]
    m_per_px = 0.0002645833; %m/px
    disk.xpos_deg = {};
    disk.ypos_deg = {};
    disk.X_m = {};
    disk.Y_m = {};
    disk.Z_m = {};
    
    for i = 1:n
        radius_deg = sqrt(disk.diskSize_deg^2 + disk.diskSize_deg^2);
        [xpos_deg, ypos_deg] = GetPointsRandom(ndp, b*2-3, b_*2-3, radius_deg);
        xpos_deg = xpos_deg + T(i);
        
        
        %xpos_deg = (b - a).* rand(ndp,1) + a + T(i); % rand range formula + offset of the plane (so the planes don't overlap)
        %ypos_deg = (b_ - a_).* rand(ndp,1) + a_;
        
        %%positions of plaids in degrees, these get passed to optic flow
        disk.xpos_deg{i} =  xpos_deg;
        disk.ypos_deg{i} =  ypos_deg;
        
        %%positions for holes in the masks for TVPMSD
        disk.X_px{i}  = (-1).*((disk.xpos_deg{i} + (2-i)*hFOV/3) .* ds.hor_px_per_deg); % convert to coordinate system of each individual plane in pixels, by going to center of each plane
        disk.Y_px{i}  = (disk.ypos_deg{i} .* ds.hor_px_per_deg);
        disk.Xpx_m = (-1).*(2 .* Z(i).* tand((disk.xpos_deg{i} + (2-i)*hFOV/3)./2));
        
        %%positions for holes in the full window mask
        disk.Xfwm_px{i} = disk.xpos_deg{i} .* ds.hor_px_per_deg;
        disk.Yfwm_px{i} = disk.ypos_deg{i} .* ds.ver_px_per_deg;
        
        %%initial positions in the world
        disk.X_m{i}  =  2 .* -Z(i).* tand(xpos_deg./2); % must flip back bc brain flips retinal image
        disk.Y_m{i}  = 2 .* -Z(i).* tand(ypos_deg./2);
        disk.Z_m{i}  = ones(ndp,1).*Z(i);%disk.Z_m{i} = ones(180,1).*Z(i);
        
% uncomment next 4 lines are for testing.
%          axis equal
%          plot3(plane.vertices{i}(1,:),plane.vertices{i}(2,:),plane.vertices{i}(3,:),'*')
%          hold on; plot3(disk.X_m{i}, disk.Y_m{i}, disk.Z_m{i},'o')
%          shg
    end
      
%     disk = load('disk.mat');
%     disk = disk.disk;
% %     
 end