

 for RecDat={'061321','062321','062521','070421','070921','071121','071321','072521','073021','080121','080421','080621','081021','081721','081921',...
        '060818','061118','061418','061518','061818','062018','062618','062718','063018','070418','070518','070818','071018','071118','071218','071718','071818','072018','072118'}; % Chico and Silas dataset
    
     RecDate=RecDat{1};
     fprintf('\n Copying Data from Rec %s',RecDate)
    if str2num(RecDate(end-1:end))==18 %if we are processing Silas data
        SummeryFile=['Rec' RecDate '_TimingData.mat'];
        SourcePath=['Z:\Projects\Rule_Representation\Data\Silas_Recording\Neural_Recording\' RecDate '\ANALYSIS\'];
    elseif str2num(RecDate(end-1:end))==21 % if we are processing Chico data
        SummeryFile=['Rec' RecDate '_TimingData.mat'];
        SourcePath=['Z:\Projects\Rule_Representation\Data\Chico_Recording\Neural_Recording\' RecDate '\ANALYSIS\'];
    end
    PathTarg='Z:\Projects\Rule_Representation\Data\Neural Summery Files\';
    copyfile([SourcePath SummeryFile],[PathTarg SummeryFile])
end