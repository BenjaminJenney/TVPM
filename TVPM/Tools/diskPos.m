function disk = diskPos(ds, disk, plane)

    hFOV   = ds.hFOV_perPersonAvg;
    vFOV   = ds.vFOV_perPersonAvg;
    T      = plane.offsets_deg;
    n      = plane.numPlanes;
    ndp    = disk.numDisksPerPlane;
    Z      = plane.depths_m;
    disk.xpos_deg = {};
    disk.ypos_deg = {};
    disk.X_m = {};
    disk.Y_m = {};
    disk.Z_m = {};
    
    for i = 1:n             
        %[xpos_deg, ypos_deg] = GetPointsRandom(ndp,  hFOV/(n) - hFOV/(2*n) + T(i),  vFOV - vFOV/2, disk.size_deg + .1);% .* hFOV/(n) - hFOV/(2*n) + T(i); 
        % the gap is the plane width divided by two in degrees. hFOV/(2*n) is the plane width divided two: centers the coordinate system at screen zero before translating
        xpos_deg = rand(ndp,1) .* hFOV/(n) - hFOV/(2*n) + T(i);
        ypos_deg = rand(ndp,1) .* vFOV - vFOV/2;
        disk.xpos_deg{i} = xpos_deg;
        disk.ypos_deg{i} = ypos_deg;
        disk.X_m{i} = 2 .* Z(i).* tand(xpos_deg./2);
        disk.Y_m{i} = 2 .* Z(i).* tand(ypos_deg./2);
        disk.Z_m{i} = ones(ndp,1).*Z(i);%disk.Z_m{i} = ones(180,1).*Z(i);
        %disk.apertureZ_m{i} = disk.Z_m{i} + .001; %place aperture .001 meters in front of each disk.
        %keyboard
    end
    disk.diskPos_m = [cat(1, disk.X_m{:}), cat(1,disk.Y_m{:}), cat(1, disk.Z_m{:})];
    
 end