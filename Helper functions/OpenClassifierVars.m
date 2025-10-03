%opens important classifier variables  
Vars2Open={'TrainDataAllFactors','Train2DataAllFactors','Train3DataAllFactors',...
     'TestDataAllFactors','Test2DataAllFactors','Test3DataAllFactors'};
cellfun(@(x) openvar(x),Vars2Open);