function AngleAnalysis
% simulate classifier angles
HeldoutRatio=0.2;
% create 4 objects with firing rate coming from
% X axis is color and Y axis is shape
%R=1 G=-1 B=1 T=-1
RB=[1 1]';RT=[1 -1]';GB=[-1 1]';GT=[-1 -1]';
ObjTxt={'RB','RT','GB','GT'};
ObjSize=[10,5,10,5];
ObjCol={'r','r','g','g'};
ObjsAll=[RB RT GB GT];
% define some math functions
AxisRotMat=@(theta) [cosd(theta) -sind(theta);sind(theta) cosd(theta)];
AxisScaleMat=@(x,y) [x 0;0 y];
XYTransMat=@(theta,x,y) [AxisRotMat(theta)]'*AxisScaleMat(x,y)*AxisRotMat(theta);
%'CongAxis',@(x,y) [0 -x;-y,0],0:0.1:2,0:0.1:2
CompCond={'CompColor','CompShape','CompColorDecompShape','CogRotAxis45','CogRotAxis20','CogRotAxisFixedComp'}; % what is the condition of compression
CompMat={@(x,y) [x 0;0,1],@(x,y) [1 0;0 y],@(x,y) [x 0;0 y],@(x,y) XYTransMat(45,x,y),@(x,y) XYTransMat(40,x,y),@(Rot,x,y) XYTransMat(Rot,x,y)      };
RngX={0.1:-0.01:-0.1    ,ones(1,21)       ,0.2:-0.01:0     ,ones(1,21)               ,ones(1,21)               ,ones(1,19)};
RngY={[ones(1,21)]      ,0.1:-0.01:-0.1    ,0:0.01:0.2     ,0:0.01:0.2                  ,0:0.01:0.2                  ,0.25*ones(1,19)   };
RngRot={[],[],[],[],[],0:5:90};

ManData=ManipulateData;

% now create compression along Color and Shape axis independatly.
% Train a classifier to discriminate shape and color and
% measure the angle of the classifier as you chnage the scaling
NSamp=100; % 100 sample per object
NoiseLvl=0.1;
NSampTest=round(NSamp*HeldoutRatio);
TrainResponse=[repmat(RB,[1 NSamp]) repmat(RT,[1 NSamp]) repmat(GB,[1 NSamp]) repmat(GT,[1 NSamp])];
validationResponse =[repmat(RB,[1 NSampTest]) repmat(RT,[1 NSampTest]) repmat(GB,[1 NSampTest]) repmat(GT,[1 NSampTest])];
ObjNum=[ones(1,NSampTest) 2*ones(1,NSampTest) 3*ones(1,NSampTest) 4*ones(1,NSampTest)];

TrainPredictorBuff=TrainResponse+NoiseLvl*randn(2,4*NSamp);
validationPredictorBuff=validationResponse+NoiseLvl*randn(2,4*NSampTest);
NCol=6;
NRow=3;
%% loop on x and y for every condition
figure
sgtitle(sprintf('Noise STD:%0.2f',NoiseLvl))
k=5;kk=1;
XLblTxt=[];
for k=1:6;%[1 2 3 4 5]
    clear Accuracy_Color Accuracy_Shape classificationSVM_Color classificationSVM_Shape AngCol AngShp AngShpCol
    for s=1:length(RngX{k})

        if isempty(RngRot{k})
            X=RngX{k}(s);Y=RngY{k}(s);
            %% Apply compression
            TrainPredictor=CompMat{k}(X,Y)*TrainPredictorBuff+NoiseLvl*randn(2,4*NSamp);
            validationPredictor=CompMat{k}(X,Y)*validationPredictorBuff+NoiseLvl*randn(2,4*NSampTest);
        else
            X=RngX{k}(s);Y=RngY{k}(s);Rot=RngRot{k}(s);
            %% Apply compression
            TrainPredictor=CompMat{k}(Rot,X,Y)*TrainPredictorBuff+NoiseLvl*randn(2,4*NSamp);
            validationPredictor=CompMat{k}(Rot,X,Y)*validationPredictorBuff+NoiseLvl*randn(2,4*NSampTest);
        end

        XLblTxt=[XLblTxt {sprintf('X%0.2fY%0.2f',X,Y)}];

        TrainPredictorBuff=TrainResponse;
        validationPredictorBuff=validationResponse;

        %% Show train and test
        subplot(NRow,NCol,k);cla
        hold on
        xlim([-2 2]);ylim([-2 2])
        arrayfun(@(x) text(validationPredictor(1,ObjNum==x),validationPredictor(2,ObjNum==x),ObjTxt{x},'color',ObjCol{x},'FontSize',ObjSize(x)),1:4);
        arrayfun(@(x) text(ObjsAll(1,x),ObjsAll(2,x),ObjTxt{x},'color',ObjCol{x},'FontSize',20),1:4);

        axis square; xlabel('Color');ylabel('Shape');title(CompCond{k})

        % train abd test on color and shape
        [Accuracy_Color(s),classificationSVM_Color{s}]=TrainTestClassifier(TrainPredictor,TrainResponse(1,:),validationPredictor,validationResponse(1,:));
        [Accuracy_Shape(s),classificationSVM_Shape{s}]=TrainTestClassifier(TrainPredictor,TrainResponse(2,:),validationPredictor,validationResponse(2,:));

        subplot(NRow,NCol,NCol+k);cla;
        hold on
        plot(1:s,Accuracy_Color,'r');
        plot(1:s,Accuracy_Shape,'b');
        if k==1
            legend('Col','Shp','Location','best')
        end

        xticks(1:2:s);
        xticklabels(XLblTxt(1:2:s));
        xtickangle(90);
        title('Classifier Accuracy');
        xlabel('Compression Step')
        ylabel('Accuracy')
        ylim([0 1])

        AngCol(s)=ManData.GetAngleBetVectors(classificationSVM_Color{s}.Beta,classificationSVM_Color{1}.Beta);
        AngShp(s)=ManData.GetAngleBetVectors(classificationSVM_Shape{s}.Beta,classificationSVM_Shape{1}.Beta);
        AngShpCol(s)=ManData.GetAngleBetVectors(classificationSVM_Color{s}.Beta,classificationSVM_Shape{s}.Beta);

        subplot(NRow,NCol,2*NCol+k);cla
        hold on

        plot(1:s,AngCol,'r');
        plot(1:s,AngShp,'b');
        plot(1:s,AngShpCol,'g');
        if k==1
            legend('Col','Shp','ColShp','Location','best')
        end

        xticks(1:2:s);
        xticklabels(XLblTxt(1:2:s));
        xtickangle(90);
        title('Classifier angle');
        xlabel('Compression Step')
        ylabel('Angle')

        mvFrame(kk) = getframe(gcf);
      %  pause
        kk=kk+1;
    end
end

end

function [Accuracy,classificationSVM,validationScores]=TrainTestClassifier(TrainPredictor,TrainResponse,validationPredictor,validationResponse)

classificationSVM=fitclinear(TrainPredictor',TrainResponse,'Learner','logistic','Lambda',1/60,'Regularization','ridge','FitBias',true,'PostFitBias',true);
[validationPredictions, validationScores] = predict(classificationSVM,validationPredictor');

% Compute validation accuracy
correctPredictions = (validationPredictions' == validationResponse);
isMissing = isnan(validationResponse);
correctPredictions = correctPredictions(~isMissing);
Accuracy = sum(correctPredictions)/length(correctPredictions);

end
