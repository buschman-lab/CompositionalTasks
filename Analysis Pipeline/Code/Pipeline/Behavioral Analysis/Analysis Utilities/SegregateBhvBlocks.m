function [BlockNums,BlockIden]=SegregateBhvBlocks(BlockSeq)

Blk1or3Ind=find(BlockSeq~=2);
Blk1or3=BlockSeq(BlockSeq~=2);
BlkChng=find(diff(Blk1or3)~=0);


if isempty(BlkChng)
        
    BlockNums{1}=1:length(BlockSeq);
        
else

        for i=1:length(BlkChng)
         %   Rule=BlockSeq(Blk1or3Ind(BlkChng(i)));

            if i==1
                BlockNums{i}=[1:Blk1or3Ind(BlkChng(i))];
            else
                BlockNums{i}=[Blk1or3Ind(BlkChng(i-1))+1:Blk1or3Ind(BlkChng(i))];            
            end
            
            IdenTemp=BlockSeq(BlockNums{i});
            BlockIden(i)=unique(IdenTemp(IdenTemp~=2));

        end

        if Blk1or3Ind(BlkChng(i))<(length(BlockSeq)-1)

         %  Rule=BlockSeq(Blk1or3Ind(BlkChng(i)));   

                BlockNums{i+1}=[Blk1or3Ind(BlkChng(i))+1:length(BlockSeq)]; 
                
                IdenTemp=BlockSeq(BlockNums{i});
               BlockIden(i+1)=unique(IdenTemp(IdenTemp~=2));

        end
end