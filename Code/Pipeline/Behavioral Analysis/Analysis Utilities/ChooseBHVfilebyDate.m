function bhv=ChooseBHVfilebyDate(InitialDate,StopDate,opts)
cd(opts.PATH);
DIR=dir;
StartDate = datetime(InitialDate);
StopDate = datetime(StopDate);

N=1;
for HH=1:length(DIR)
    [~,FileName,Ext]=fileparts(DIR(HH).name);
    if strcmp(Ext,'.mat') & (strcmp(FileName(1:5),'Chico') || strcmp(FileName(1:5),'Silas'))
        ThisFileDate=datetime(2000+str2double(FileName(7:8)),str2double(FileName(9:10)),str2double(FileName(11:12)));
        if ThisFileDate>=StartDate && ThisFileDate<=StopDate
            %  if DIR(HH).date>=StartDate && DIR(HH).date<=StopDate && strcmp(Ext,'.mat') && ~isempty(findstr(DIR(HH).name,'Chico'))
            disp(['Loading ....' char(DIR(HH).name)])
            bhv{N}=load(DIR(HH).name,'bhv');
            bhv{N}.bhv.FigNum=N;
            bhv{N}.Date=DIR(HH).date;
            N=N+1;
        end
    end
    close all
end
cd(opts.CodePath);