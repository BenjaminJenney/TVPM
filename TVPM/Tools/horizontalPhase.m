function yGrating = horizontalGrating(v, dty, ifi, deg_per_px)
    sf =.1; % cycles per deg
    sf = (deg_per_px).*sf; %1/px_per_deg %convert from deg to 
    %Vx = (1/deg_per_px).*u; %(1/mmPerPixel)*(tand(u/2)*2*viewingDistance); % (dg/s) -> px/s
    Vy = (1/deg_per_px).*v;%1/mmPerPixel)*(tand(v/2)*2*viewingDistance); 

    %pa.dtx = pa.dtx + ((2 * pi * (Vx * (1/1000) * sf)) / ifi); % px/s * s/ms * cyc/px * ms/frame = cyc/frame
    dty =  dty + ((2 * pi * (Vy * (1/1000) * sf)) / ifi);
    yGrating = (.5 + .5 * cos(af * y - dty));
end