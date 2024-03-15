% sparses file name and find out it charctristics
function Specs=SparseFileName(Path)

[Specs.path,Specs.Name,Specs.Ext]=fileparts(Path);

RecInd=findstr(Specs.Name,'Rec');
ChInd=findstr(Specs.Name,'Channel');
ThInd=findstr(lower(Specs.Name),'th');
ManulInd=findstr(Specs.Name,'Manu');

Specs.RecDate=GrabString(Specs.Name,RecInd,3,6);
Specs.Ch=extractNumFromStr(GrabString(Specs.Name,ChInd,8,3));
Specs.Th=extractNumFromStr(GrabString(Specs.Name,ThInd,2,3));
Specs.ManulSpk=~isempty(ManulInd);


end
function str=GrabString(String,StrInd,Dent,Length)
if ~isempty(StrInd)
    str=String(StrInd+Dent:StrInd+Dent+Length-1);
else
    str='';
end
end