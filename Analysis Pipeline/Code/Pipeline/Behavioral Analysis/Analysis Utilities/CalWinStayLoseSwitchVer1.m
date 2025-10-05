function I_WinStayLoseSwitch=CalWinStayLoseSwitchVer1(i,NSwitch,Trials,Nback,TrialInterval,fh)
% looks at all of the trials of a block and calculates the win stay lose
% switch index during a block%   corrects for next trial's identiry
if isempty(fh); fh = figure(10000); end
figure(fh)
set(0, 'currentfigure', fh);
hold on;
sm_kern = ones(1, Nback); sm_kern = sm_kern./sum(sm_kern);

INDOBJ_AllTrials=(([Trials(TrialInterval).StopCondition] == 1 | [ Trials(TrialInterval).StopCondition] == -1));
INDOBJ_CorrectTrials=([Trials(TrialInterval).StopCondition] == 1);
INDOBJ_InCorrectTrials=([Trials(TrialInterval).StopCondition] == -1);
Response=([Trials(TrialInterval).Response]);


CorrectMat=INDOBJ_CorrectTrials(INDOBJ_AllTrials); 
InCorrectMat=INDOBJ_InCorrectTrials(INDOBJ_AllTrials);
ResponseMat=Response(INDOBJ_AllTrials);

CorrectMatWin=CorrectMat(1:end-1);
CorrectMatWinP1=CorrectMat(2:end);
ConjucateMatWS=(CorrectMatWin==1 & CorrectMatWinP1==0);
RespWin=ResponseMat(ConjucateMatWS);
RespWinP1=ResponseMat(find(ConjucateMatWS==1)+1);

eqmatWinStay= RespWin ==RespWinP1;
I_WinStay=convn(eqmatWinStay,sm_kern,'valid');

%%

InCorrectMatLose=InCorrectMat(1:end-1);
InCorrectMatLoseP1=InCorrectMat(2:end);
ConjucateMatLS=(InCorrectMatLose==1 & InCorrectMatLoseP1==1);
RespLose=ResponseMat(ConjucateMatLS);
RespLoseP1=ResponseMat(find(ConjucateMatLS==1)+1);

eqmatLoseSwitch= RespLose ~=RespLoseP1;

I_LoseSwitch=convn(eqmatLoseSwitch,sm_kern,'valid');
%%

 

eqmat_WinStayLoseSwitch=zeros(1,length(InCorrectMat)-1);
eqmat_WinStayLoseSwitch(ConjucateMatWS)=eqmatWinStay;
eqmat_WinStayLoseSwitch(ConjucateMatLS)=eqmatLoseSwitch;

I_WinStayLoseSwitch=convn(eqmat_WinStayLoseSwitch,sm_kern,'valid');

  