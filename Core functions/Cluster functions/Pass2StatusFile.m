%%% this function creates a status file which shows where we are in our
%%% runs in the server
function Pass2StatusFile(StatusFileName,SNum,RecDate,Ch,IsFinished)
global StatusMatFile
if isempty(StatusMatFile)
    StatusMatFile=matfile([StatusFileName '.mat'],'Writable',true);
end
StatusMatFile.ChannelDone(1,SNum)=Ch;
StatusMatFile.RecDate(1,SNum)=str2num(RecDate);
StatusMatFile.IsFinished(1,SNum)=IsFinished;
fprintf('\nUpdated file status...')
