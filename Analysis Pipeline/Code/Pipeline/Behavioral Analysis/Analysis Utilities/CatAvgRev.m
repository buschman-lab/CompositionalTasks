

function Ifin=CatAvgRev(Mat)
    for i=1:size(Mat,1)
        for ii=1:size(Mat,2)
            S(i,ii)=length(Mat{i,ii});
        end
    end
 
    MaxL=max(S,[],2);Ifin{1}=[];Ifin{2}=[];Ifin{3}=[];
    
    for i=1:size(Mat,1)
        for ii=1:size(Mat,2)
            tempMat=padcat(Mat{i,ii},ones(1,MaxL(i)));
            A=tempMat(1,:);
            A=[A(isnan(A)==1) A(isnan(A)~=1)];
            Ifin{i}=[Ifin{i};A];
         end
    end
end