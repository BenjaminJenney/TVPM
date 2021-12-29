function plaidData = staticPlaidData(ds, x, y)
    white = 1.0;
    grey  = white/2;
    inc   = white - grey;
    
    sf =.1; % cycles per deg
    sf = (ds.deg_per_px).*sf; %1/px_per_deg %convert from deg to 
    af = 2*pi*sf;
    
    xGrating = (grey + inc * sin(af * x)); % vertically oriented grating
    yGrating = (grey + inc * cos(af * y)); % horiz grating
    plaidData = 255.*((xGrating + yGrating)/2);
    plaidData = repmat(plaidData,[ 1 1 3 ]);
    plaidData = permute(uint8(plaidData),[ 3 2 1 ]);
    plaidData(4,:,:) = 255;
        
    plaidData(4,:,628:668) = 0;
    
    %{
    figure; imagesc(255.*((xGrating + yGrating)/2))
    axis equal
    colormap gray
    shg
    keyboard
%}
end