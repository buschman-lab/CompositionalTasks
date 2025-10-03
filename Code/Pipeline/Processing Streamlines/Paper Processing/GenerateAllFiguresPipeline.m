function GenerateAllFiguresPipeline(Analysis,Animal,AreaNum,RunMode,PSTHBin,TestGroup,TestName,StpStrTime,varargin) %pipeline to generate all of the figures and data for these figures for analysis
%% for classifier 
% @RunMode:0 dont run shuffle, run main and plot; 1= run missing files main and shuffle and then plot
% 2=plot only 3=Run main, no shuffle no plot 4=run missing main only no shuffle or plot (in subspace it is run shuffle only)
% 5: plot comparision only 6: collect comparision results
% 7: GLM: plot summery for each area 8: concatinate remaining files 9: run shuffle only  
% 10: runs Targeted Dim reduction analysis
% all of these figure configurations go together with

%% for GLM 
% Runmode 0: fit each neuron without shuffle (PS11)
% Runmode 1: fit each neuron with it's shuffle (PS11), concatinate (PS15) it and create a summery(PS14) for each neuron(No Plots) 
% Runmode 2: Plot compact singlecell cahractristics with GLM(PS10)
% Runmode 5: run model comparision(PS12)
% Runmode 6 plot model comparision results(PS13)
% Runmode 7 plot all GLM result summery for each area(PS16)
% Runmode 8 plot concatinate GLM result summery for each area(PS14)
% Runmode 9  concatinate GLM result summery for each neuron(PS15)

% @TestGroup: a cell array indicating the test group
% 'Decoding','CrossDecoding','Learning','Comparision','FineTune'
% StpStrTime cell array with values SACCADE_START SAMPLE_ON
% PSTHBin can be 50 100 150
% Analysis: can be 'GLM', 'Classifier','Subspace'
global AnalysisOpts
CL=ClusterFuncsTemp;
FigParams=fig_params;
SetAnalysisOptions_RuleRepresentation([varargin,'RunMode',RunMode]) % set all of the analysis functions

if AnalysisOpts.ThisIsSinaPC
    input(['Are you sure you want to use ' AnalysisOpts.ClassifierFunctiononClust],'s');
end

if contains(AnalysisOpts.ClassifierFunctiononClust,'Temp')
    NeuAnaFunc=NeuralAnalysisFuncsTemp;    
else
    NeuAnaFunc=NeuralAnalysisFuncs;
end

if ispc
    addpath(' Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Input Output Data\Core Data');
else
    addpath('/Volumes/buschman/Projects/Rule_Representation/ElecPhys_Analysis/Rule Representation Project/Analysis Pipeline/Input Output Data/Core Data');
end
if isempty(AnalysisOpts.SpockPass)
    if ~AnalysisOpts.ThisIsSinaPC
        % get the credentials
        UsrNm=inputdlg({'EnterUserName(a:Adel,s:Sina,q:Qinpu,p:Polina'},'Spock Username');
        if contains(UsrNm,'a','IgnoreCase',1);AnalysisOpts.username='aa1749';
        elseif contains(UsrNm,'s','IgnoreCase',1);AnalysisOpts.username='tafazoli';
        elseif contains(UsrNm,'q','IgnoreCase',1);AnalysisOpts.username='qh8777';
        elseif contains(UsrNm,'p','IgnoreCase',1);AnalysisOpts.username='pi7127';
        end
    else
        AnalysisOpts.username='tafazoli';UsrNm='tafazoli';
    end
    if contains(UsrNm,'s','IgnoreCase',1) | AnalysisOpts.ThisIsSinaPC
        AnalysisOpts.SpockPass='PoqttirwpirPU2084!';
    else
        AnalysisOpts.SpockPass=passwordEntryDialog('PasswordLengthMax',200,'WindowName','Enter your spock password');
    end
end
%% prepare variables
if isempty(StpStrTime);StpStrTime={'SAMPLE_ON','SACCADE_START'};end
if ~isempty(TestName);varargin=[varargin {sprintf('PaperFigs.%s.TestName',TestGroup{1}),[TestName]}];end

switch Analysis
     %% %%%%%%%%%%%%%%%%%%%%%%generate GLM results%%%%%%%%%%%%%%%%%%%%%%%%
    case 'GLM'
        switch RunMode
            case {[0],[1],[2],[5],[6],[7],[8],[9]}  % fit GLM models and their shuffles
                TestGroup={'GLM'};
                arrayfun(@(spkstr) arrayfun(@(Area) GenerateAllPaperFigures_SingleCellAnalysis(Animal,Area,StpStrTime(spkstr),PSTHBin,TestGroup,RunMode,varargin{:}),AreaNum),1:length(StpStrTime));
        end
    case 'Classifier'  % or TDR analysis
        %% %%%%%%%%%%%%%%%%%%generate classifier results%%%%%%%%%%%%%%%%%%v
        fprintf('\nGenerating classifier data and figures...')
        switch RunMode
            case {[0],[1],[2],[3],[4],[8],[9],[10]}  % genrate data or plot data for groups of conditions
                if isempty(TestGroup);TestGroup={'Decoding','CrossDecoding','FineTune'};end
                arrayfun(@(spkstr) arrayfun(@(Area) GenerateAllPaperFigures_Classifier(Animal,Area,StpStrTime(spkstr),PSTHBin,TestGroup,RunMode,varargin{:}),AreaNum),1:length(StpStrTime));
            case 5 % plot comparision
                arrayfun(@(spkstr) arrayfun(@(Area) GenerateAllPaperFigures_Classifier(Animal,Area,StpStrTime(spkstr),PSTHBin,{'Comparision'},RunMode,varargin{:}),1),1:length(StpStrTime));
            case 6   %% collate all of the figures generated from comparision and put them in a big image file
                for spkstr=1:length(StpStrTime)
                    GenerateAllPaperFigures_Classifier(Animal,1,StpStrTime(spkstr),PSTHBin,{'Comparision'},RunMode,varargin{:});
                    AnalysisOpts.Animal=Animal;
                    [FigFileNames,OutFileName]=NeuAnaFunc.GetClassifierAnaTestOutputImgFileName('Classifier',...
                        [AnalysisOpts.PaperFigs.Comparision.TestName],...
                        1,[],'tif','Comparision',AnalysisOpts.PaperFigs.Comparision.SpkCntStartFieldName{1},AnalysisOpts.PaperFigs.Comparision.PSTHbin);
                    FigParams.CollectFiguresandSave(FigFileNames,OutFileName(1:end-4),'tif');
                end
            case 7 %collate all figures for a specific analysis and put the figures together 
                % define AnalysisOpts.Page2SaveNum as page numbers to be
                % taken from each analysis
                AnalysisOpts.Animal=Animal;
                for spkstr=1:length(StpStrTime)                  
                    for tst=1:length(TestName)
                        [FigFileNames,OutFileName]=NeuAnaFunc.GetClassifierAnaTestOutputImgFileName('Classifier',...
                            TestName(tst),AreaNum,AnalysisOpts.Page2SaveNum{tst},'tif',TestGroup{1},StpStrTime{spkstr},PSTHBin);
                        FigParams.CollectFiguresandSave(FigFileNames,OutFileName(1:end-4),'tif',AnalysisOpts.Page2SaveNum{tst});
                    end
                end
        end
         %% %%%%%%%%%%%%%%%%%%%%%%generate Subspace results%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Subspace'
        switch RunMode
            case {[0],[1],[2],[3],[4]}                
                arrayfun(@(spkstr) arrayfun(@(Area) GenerateAllPaperFigures_SubspaceGeometry(Animal,Area,StpStrTime(spkstr),PSTHBin,TestGroup,RunMode,varargin{:}),AreaNum),1:length(StpStrTime));
        end        
end
end

