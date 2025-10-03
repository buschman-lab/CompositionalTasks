function I_WinStayLoseSwitch=CalWinStayLoseSwitch(i,NSwitch,Trials,Nback,TrialInterval,fh)
% looks at all of the trials of a block and calculates the win stay lose
% switch index during a block
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
RespWin=ResponseMat(CorrectMatWin);
RespWinP1=ResponseMat(find(CorrectMatWin)+1);
eqmatWinStay= RespWin ==RespWinP1;
I_WinStay=convn(eqmatWinStay,sm_kern,'valid');

%%
InCorrectMatLose=InCorrectMat(1:end-1);
RespLose=ResponseMat(InCorrectMatLose);
RespLoseP1=ResponseMat(find(InCorrectMatLose)+1);
eqmatLoseSwitch= RespLose ~=RespLoseP1;
I_LoseSwitch=convn(eqmatLoseSwitch,sm_kern,'valid');
 %%

eqmat_WinStayLoseSwitch=zeros(1,length(InCorrectMat)-1);
eqmat_WinStayLoseSwitch(CorrectMatWin)=eqmatWinStay;
eqmat_WinStayLoseSwitch(InCorrectMatLose)=eqmatLoseSwitch;

I_WinStayLoseSwitch=convn(eqmat_WinStayLoseSwitch,sm_kern,'valid');

  