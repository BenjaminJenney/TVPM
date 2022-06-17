if movieframe_n == 1  % save first frame of each movie (not the rest for speed)
     rect = [];
     % added to save movie clip
     M = Screen('GetImage', window, rect, [], 0, 1);
     imwrite(M,[moviepath, '/Image_',num2str(movieframe_n),'.png']);
     movieframe_n = movieframe_n + 1;
end
