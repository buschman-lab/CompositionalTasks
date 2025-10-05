function [AreaName,ChennelNumbs]=FindChannelArea(TableFile,Date,Channel)

Table=readtable(TableFile);
DateRaw=find(strcmp(Table.Date,Date));
TableCell=table2cell(Table);
ChannelGroups=2:10;

for i=1:length(ChannelGroups)
    CH=ChannelGroups(i);
  %  ChennelNumbs{i}=arrayfun(@(x) (TableCell(DateRaw,x)),CH,'UniformOutput',false);
    Chans=arrayfun(@(x) (cell2mat(x)),TableCell(DateRaw,CH),'UniformOutput',false);
    ChennelNumbs{i}=[];
    for j=1:length(Chans)
        if ischar(Chans{j})
            ChennelNumbs{i}=cat(2,ChennelNumbs{i},str2num(Chans{j}));
        else
            ChennelNumbs{i}=cat(2,ChennelNumbs{i},Chans{j});            
        end
    end
end
AreaInd=cellfun(@(x) (sum(x==Channel)),ChennelNumbs, 'UniformOutput', 1);
if sum(AreaInd)~=0
    AreaName=cell2mat(Table.Properties.VariableNames(find(AreaInd>0)+1));
else
    AreaName='None';
end


