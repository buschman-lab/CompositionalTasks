function MakeMovieFromFrames(mvFrame,fps,saveMovieStr)
% create move from the frames 
% inputs 
% MvFrame       Frames captured with     mvFrame(i) = getframe(gcf);
%fps            Frame per second
%saveMovieStr   Save Path 
% use     mvFrame(cur_win) = getframe(gcf);

saveMovieStr = [saveMovieStr, '.avi'];

writerObj = VideoWriter(saveMovieStr, 'Motion JPEG AVI');
writerObj.FrameRate = fps;
open(writerObj); 
fprintf('writing movie...\n')
mvFrame = mvFrame(~cellfun(@isempty,{mvFrame.cdata})); %remove dropped frames
try
    for i =1:length(mvFrame)
        if isempty(mvFrame(i).cdata)
            mvFrame(i) = [];
        end
    end
catch
end
writeVideo(writerObj, mvFrame);
close(writerObj);