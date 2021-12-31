function disk = diskPos(ds, disk, plane)

    hFOV   = ds.hFOV_perPersonAvg;
    vFOV   = ds.vFOV_perPersonAvg;
    T      = plane.offsets_deg;
    n      = plane.numPlanes;
    Z      = plane.depths_m;
    disk.xpos_deg = {};
    disk.ypos_deg = {};
    disk.X_m = {};
    disk.Y_m = {};
    disk.Z_m = {};
    for i = 1:n             
        xpos_deg = rand(disk.numDisks,1) .* hFOV/(n) - hFOV/(2*n) + T(i); 
        % the gap is the plane width divided by two in degrees. hFOV/(2*n) is the plane width divided two: centers the coordinate system at screen zero before translating
        ypos_deg = rand(disk.numDisks,1) .* vFOV - vFOV/2;
        disk.xpos_deg{i} = xpos_deg;
        disk.ypos_deg{i} = ypos_deg;
        disk.X_m{i} = 2 .* Z(i).* tand(xpos_deg./2);
        disk.Y_m{i} = 2 .* Z(i).* tand(ypos_deg./2);
        disk.Z_m{i} = ones(disk.numDisks,1).*Z(i);
        disk.apertureX_m{i} = disk.X_m{i};
        disk.apertureY_m{i} = disk.Y_m{i};
        disk.apertureZ_m{i} = disk.Z_m{i} + .01; %place aperture .001 meters in front of each disk.
        keyboard
    end
end