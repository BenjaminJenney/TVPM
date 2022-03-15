function drawTheHoleyBag(fullWindowMask, ds, openglIsEnabled)
    if nargin == 1
        openglIsEnabled = 1; 
    end
    
    if openglIsEnabled == 1
        Screen('EndOpenGL', ds.w);
        %Screen('DrawTextures', ds.fullWindowMask, ds.masktex, [], ds.dstRects); % see setupDisplay for rest of code doing this
        Screen('DrawTexture', ds.w, fullWindowMask, [], [], [], [], []);
%         Screen('DrawTextures', ds.fullWindowMask, ds.masktex, [], ds.dstRects); % see setupDisplay for rest of code doing this
%         Screen('DrawTexture', ds.w, ds.fullWindowMask, [], [], [], [], 1);
        Screen('BeginOpenGL', ds.w);
        return;
    end
    
    Screen('DrawTextures', fullWindowMask, ds.masktex, [], ds.dstRects); % see setupDisplay for rest of code doing this
    Screen('DrawTexture', ds.w, fullWindowMask, [], [], [], [], 1);  
    Screen('DrawTextures', fullWindowMask, ds.masktex, [], ds.dstRects); % see setupDisplay for rest of code doing this
    Screen('DrawTexture', ds.w, fullWindowMask, [], [], [], [], 1);
end