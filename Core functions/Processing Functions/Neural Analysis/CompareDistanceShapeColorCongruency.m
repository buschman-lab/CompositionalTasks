%Compares distance for shape and color in congruency
Var2Look={'ShapeDistCong','ShapeDistInCong','ColorDistCong','ColorDistInCong'};

% find these variables 
IndVars=cellfun(@(x) find(strcmp(EncodingDistVars,x)),Var2Look);

EncodingDistThisVars=EncodingDist(IndVars);

EncodingDistThisVarsRS=cellfun(@(x) obj.ManData.ReshapeCell2Mat(x,62),EncodingDistThisVars,'UniformOutput',0);

EncodingDistThisVarsRSMean=cellfun(@(x) mean(x(:,Time2LookDim{1}),2),EncodingDistThisVarsRS,'UniformOutput',0);
figure
hold on
Col=copper(length(EncodingDistThisVarsRSMean{1}));
arrayfun(@(x) plot(EncodingDistThisVarsRSMean{1}(x),EncodingDistThisVarsRSMean{3}(x),'*','MarkerSize',15,'Color',Col(x,:)),1:length(EncodingDistThisVarsRSMean{1}));
arrayfun(@(x) plot(EncodingDistThisVarsRSMean{2}(x),EncodingDistThisVarsRSMean{4}(x),'d','MarkerSize',15,'Color',Col(x,:)),1:length(EncodingDistThisVarsRSMean{1}));
axis square
xlabel('Distance Shape')
ylabel('Distance Color')
%legend('Congruent','Incongruent')


