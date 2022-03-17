function diskPosNew()
    %declare range
    b = 14.5; %b = hFOV/(2*numplanes), width of the plane = hFOV/3
    a = -14.5;
    b_ = 42;
    a_ = -42;
    depths_m = -.5;
    xpos_deg = (b - a).*rand(100, 1) + a - 29;
    ypos_deg = (b_ - a_) .*rand(100,1) + a_;
    xpos_m   = 2 .* depths_m .* tand(xpos_deg./2);
    ypos_m   = 2 .* depths_m .* tand(ypos_deg./2);
    
    widths_m = 2 * -depths_m * tand((87/3)/2);
    heights_m = 2 * -depths_m * tand(84/2);
    %ypos_m = 2 .* Z(i).* tand(xpos_deg./2);
    ithPlaneVertices =[-widths_m/2 -heights_m/2 depths_m ;... 
                            widths_m/2 -heights_m/2 depths_m ;...
                            widths_m/2  heights_m/2 depths_m ;...
                           -widths_m/2  heights_m/2 depths_m ]';
    ithPlaneVertices = ithPlaneVertices - [-widths_m; 0; 0;];
axis equal;
plot3(ithPlaneVertices(1,:), ithPlaneVertices(2,:), ithPlaneVertices(3,:), '*'); hold on; plot3(xpos_m, ypos_m, ones(100, 1).*(depths_m), '.');


end