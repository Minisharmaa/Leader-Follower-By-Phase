function [CoherenceValues]=LeaderFollowerByPhase(signal1,signal2,varargin)

%Argument order:
%signal1: NX1 first signal
%signal2: NX1  second signal
%lowFreq: The default value is 0.01 Hz
%highFreq: The default value is 0.5 Hz
%phaseRange: The default value is 90 degrees
%Threshold: The default value is 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Argument Description:
%The two signals must be equal.
%lowFreq and highFreq: between 0 to 1 Hz.
%phaseRange: range of phases for all condition (synch, signal1 leads,
%signal2 leads, negative synch). The range of value is between 0 to 90. 
%For example: If phaseRange is set to 60, then the phase range for synchronize will be between -30 to 30; The phase range for
%signal1 leads -120 to -60; The phase range for signal2 leads 60 to 120; The phase range for
%negative synch 150 to 180 and  -180 to -150. 
%Threshold: The cut of color in wtc graph (minimum Rsq value).The range is
%between 0 to 1. It is equivalent to R^2. Hence, if you want calculate the R value you should perform sqrt.
%For example: If threshold is set to 0.49, that means the R value is 0.7. 
%
%Examples how to implement this function:
%
%In this examples values are set by the user:
%cohervalues=LeaderFollowerByPhase(signal1,signal2,lowFreq,highFreq,phaseRange,threshold)
%cohervalues=LeaderFollowerByPhase(signal1,signal2,0.002,0.1,60,0.49)
%
%In these examples the function uses default values:
%cohervalues=LeaderFollowerByPhase(signal1,signal2);value of all variables
%except for signal1 and signal2 will be determined according to their default values.
%cohervalues=LeaderFollowerByPhase(signal1,signal2,[],[],60,0.3); the values of lowFreq and highFreq will be determined according to their default values.
%
%The output of the function include:
%5 graphs :
% (1) Box plot (presenting median, min and max values of Rsq for all interacition [synch,signal1 leads,
%     signal2 leads,negative synch])
% (2) Bar plot (presenting mean, max and median values of Rsq for all interacition [synch,signal1 leads,
%     signal2 leads,negative synch])
% (3) Scatter plot (Rsq by time according to the different types of interacition [synch,signal1 leads,
%     signal2 leads,negative synch])  
% (4) Pie plot (percentages of samples by type of interaction [synch,signal1 leads,
%     signal2 leads,negative synch])

%The output also contains table with statistic values of Rsq according to type of interaction (synch,signal1 leads,
%signal2 leads,negative synch) (i.e., datatable). 


%------------------------------------------------------------------------------------------------------------------

% min and max variable in varargin option
narginchk(2,6)

% Default values in cases of empty values or missing values
if numel(varargin)==0
    varargin{1}=[];
    varargin{2}=[];
    varargin{3}=[];
    varargin{4}=[];
end
if numel(varargin)==1
    varargin{2}=[];
    varargin{3}=[];
    varargin{4}=[];
end
if numel(varargin)==2
    varargin{3}=[];
    varargin{4}=[];
end
if numel(varargin)==3
    varargin{4}=[];
end
    if isempty(varargin{1})
        varargin{1}=0.01;
    end
    if isempty(varargin{2})
        varargin{2}=0.5;
    end
    if isempty(varargin{3})
        varargin{3}=90;
    end
    if isempty(varargin{4})
        varargin{4}=0;
    end
lowFreq=varargin{1};
highFreq=varargin{2};
phaseRange=varargin{3};
threshold=varargin{4};

% check list for the signals variables
    if size(signal1,1)==1
        signal1=signal1';
    end
    if size(signal2,1)==1
        signal2=signal2';
    end
    if length(signal1)~=length(signal2)
        error("The length of the two signals must to be equal")
    end
%check list for variable at position 1 in varargin
    if numel(varargin{1})~=1
        error("The dimention of lowFreq variable must to be (1,1)")
    end
    if varargin{1}<0 || varargin{1}>1
        error("The value of lowFreq variable must to be between 0 to 1")
    end
%check list for variable at position 2 in varargin
    if numel(varargin{2})~=1
        error("The dimension of highFreq variable must to be (1,1)")
    end
    if varargin{2}<0 || varargin{2}>1
        error("The value of highFreq variable must to be between 0 to 1")
    end
    if varargin{2}<varargin{1}
        error("The value of highFreq variable must to bigger from the value of lowFreq variable")
    end
%check list for variable at position 3 in varargin
     if numel(varargin{3})~=1
        error("The dimension of phaseRange variable must to be (1,1)")
    end
    if varargin{3}<=0 || varargin{3}>90
        error("The value of phaseRange variable must to be between 0 to 90")
    end
%check list for variable at position 4 in varargin
     if numel(varargin{4})~=1
        error("The dimension of threshold variable must to be (1,1)")
    end
    if varargin{4}<0 || varargin{4}>1
        error("The value of threshold variable must to be between 0 to 1")
    end

%calculate min/max phase of the range
minPhase=-abs(phaseRange/2);
maxPhase=abs(phaseRange/2);

% calculate coherence, phase and frequency of the signals 
[wrsq,period,scale,coi,sig95]=wtc(signal1,signal2,'mcc',0);
[Wxy,Wperiod,Wscale,Wcoi,Wsig95]=xwt(signal1,signal2);


%choose the FOI range
%(example fs=7.8125 Hz).
period=period';
frequency=1./period;
rangeFreq=find(frequency>=lowFreq &frequency<=highFreq);
minfreq=rangeFreq(end);
maxFreq=rangeFreq(1);

%cut the Rsq and phase (complex number) according FOI
wcohFOI=wrsq(maxFreq:minfreq,:);
WxyFOI=Wxy(maxFreq:minfreq,:);

%use only data with Rsq>=threshold
sigCoh=wcohFOI>=threshold;
[R,C]=size(sigCoh);
for i=1:C
    tr=0;
    for j=1:R
        if sigCoh(j,i)==1
            tr=tr+1;
            temp(tr,1)=wcohFOI(j,i);
            temp2(tr,1)=WxyFOI(j,i);
        end
    end
    Avgcoh(:,i)=mean(temp,1);
    AvgComplex(:,i)=mean(temp2,1);
    temp=temp();
    temp2=temp2();
end


% Find Synchronize values
SynchLoc=find(atan2d(imag(AvgComplex),real(AvgComplex))>= minPhase & ...
              atan2d(imag(AvgComplex),real(AvgComplex))<=maxPhase);
if isempty(SynchLoc)
    cohSynch(1)=NaN;
    cohSynch(2)=NaN;
    cohSynch(3)=NaN;
    cohSynch(4)=NaN;
else
    cohSynch(1)=mean(Avgcoh(SynchLoc),'omitnan');
    cohSynch(2)=max(Avgcoh(SynchLoc));
    cohSynch(3)=median(Avgcoh(SynchLoc),'omitnan');
    cohSynch(4)=std(Avgcoh(SynchLoc),'omitnan');
end

% Find signal1 leader and signal2 follower 
signal1LeadLoc=find(atan2d(imag(AvgComplex),real(AvgComplex))>= (-90+minPhase) & ...
              atan2d(imag(AvgComplex),real(AvgComplex))<=(-90+maxPhase));
if isempty(signal1LeadLoc)
    coh1Lead(1)=NaN;
    coh1Lead(2)=NaN;
    coh1Lead(3)=NaN;
    coh1Lead(4)=NaN;
else          
    coh1Lead(1)=mean(Avgcoh(signal1LeadLoc),'omitnan');
    coh1Lead(2)=max(Avgcoh(signal1LeadLoc));
    coh1Lead(3)=median(Avgcoh(signal1LeadLoc),'omitnan');
    coh1Lead(4)=std(Avgcoh(signal1LeadLoc),'omitnan');
end

% Find signal2 leader and signal1 follower 
signal2LeadLoc=find(atan2d(imag(AvgComplex),real(AvgComplex))>= (90+minPhase) & ...
              atan2d(imag(AvgComplex),real(AvgComplex))<=(90+maxPhase));
if isempty(signal2LeadLoc)
    coh2Lead(1)=NaN;
    coh2Lead(2)=NaN;
    coh2Lead(3)=NaN;
    coh2Lead(4)=NaN;
else
    coh2Lead(1)=mean(Avgcoh(signal2LeadLoc),'omitnan');
    coh2Lead(2)=max(Avgcoh(signal2LeadLoc));
    coh2Lead(3)=median(Avgcoh(signal2LeadLoc),'omitnan');
    coh2Lead(4)=std(Avgcoh(signal2LeadLoc),'omitnan');
end

% Find Negative Synchronize values
NegSynchLoc=find((atan2d(imag(AvgComplex),real(AvgComplex))>= (-180) & ...
              atan2d(imag(AvgComplex),real(AvgComplex))<=(-180-minPhase))|...
              (atan2d(imag(AvgComplex),real(AvgComplex))>= (180-maxPhase) & ...
              atan2d(imag(AvgComplex),real(AvgComplex))<=(180)));
if isempty(NegSynchLoc)
    cohNegSynch(1)=NaN;
    cohNegSynch(2)=NaN;
    cohNegSynch(3)=NaN;
    cohNegSynch(4)=NaN;
else          
    cohNegSynch(1)=mean(Avgcoh(NegSynchLoc),'omitnan');
    cohNegSynch(2)=max(Avgcoh(NegSynchLoc));
    cohNegSynch(3)=median(Avgcoh(NegSynchLoc),'omitnan');
    cohNegSynch(4)=std(Avgcoh(NegSynchLoc),'omitnan');
end

% Prepair data for plots
VecSize(1)=numel(Avgcoh(SynchLoc));
VecSize(2)=numel(Avgcoh(signal1LeadLoc));
VecSize(3)=numel(Avgcoh(signal2LeadLoc));
VecSize(4)=numel(Avgcoh(NegSynchLoc));

% Calcute time percentage of any condition 
cohSynch(5)=VecSize(1)/C;
coh1Lead(5)=VecSize(2)/C;
coh2Lead(5)=VecSize(3)/C;
cohNegSynch(5)=VecSize(4)/C;

% Continue prepair data for plots
SynchCondition={'In-Phase';'Signal1Leader';'Signal2Leader';'Anti-Phase'};


ConditionVec(1:VecSize(1))=repmat(SynchCondition(1),1,VecSize(1));
ConditionVec((VecSize(1)+1):(VecSize(1)+VecSize(2)))=repmat(SynchCondition(2),1,VecSize(2));
ConditionVec((VecSize(1)+VecSize(2)+1):(VecSize(1)+VecSize(2)+VecSize(3)))=repmat(SynchCondition(3),1,VecSize(3));
ConditionVec((VecSize(1)+VecSize(2)+VecSize(3)+1):(VecSize(1)+VecSize(2)+VecSize(3)+VecSize(4)))=repmat(SynchCondition(4),1,VecSize(4));

CohForPlot=[Avgcoh(SynchLoc),Avgcoh(signal1LeadLoc),Avgcoh(signal2LeadLoc),...
    Avgcoh(NegSynchLoc)];
colors=[0 0 0;...
    0.35 0.35 0.35;...
    0.5 0.5 0.5;...
    0.4940 0.1840 0.5560];

CoherenceValues=[cohSynch;coh1Lead;coh2Lead;cohNegSynch];

datatable=table(cohSynch',coh1Lead',coh2Lead',cohNegSynch',...
    'VariableNames',{'In-Phase';'Signal1Leader';'Signal2Leader';'Anti-Phase'},...
    'RowNames',{'mean','max','median','std','percentage'});


% check if exist value to plot
ex=~isnan(CoherenceValues(:,1));
ro=0;
for k=1:size(ex)
    if ex(k)==1
        ro=ro+1;
        existValue(ro)=k;
    end
end
        
%print the values of coherence
datatable

%Box plot
%figure;
subplot(2,2,1)
boxplot(CohForPlot,ConditionVec,'Color',colors(existValue,:))
xlabel('Interaction')
ylabel('Coherence')
title('Coherence by Type  of Interaction')
text(0,0.8,'A','FontSize',14)
hold off

%Bar plot
%figure;
subplot(2,2,2)
BB=bar(categorical({'Mean','Max','Median'}),CoherenceValues(:,1:3)');
for ii=1:4
    BB(ii).FaceColor=colors(ii,:);
end
ylabel('Coherence')
title('Centeral Indices by Type of Interacition')
legend(SynchCondition)
hold off

%Scatter plot
LocVec=[SynchLoc,signal1LeadLoc,signal2LeadLoc,NegSynchLoc];
cl=zeros(C,3);
for i=1:VecSize(1)
    cl(SynchLoc(i),:)=colors(1,:);
end
for j=1:VecSize(2)
    cl(signal1LeadLoc(j),:)=colors(2,:);
end
for k=1:VecSize(3)
    cl(signal2LeadLoc(k),:)=colors(3,:);
end
for m=1:VecSize(4)
    cl(NegSynchLoc(m),:)=colors(4,:);
end
sampleLoc=1:C;

%figure;
subplot(2,2,3)
if ex(1)==1 
    scatter(sampleLoc(SynchLoc),Avgcoh(SynchLoc),9,cl(SynchLoc,:),"filled");
end
hold on
if ex(2)==1
    scatter(sampleLoc(signal1LeadLoc),Avgcoh(signal1LeadLoc),9,cl(signal1LeadLoc,:),"filled");
end
if ex(3)==1
    scatter(sampleLoc(signal2LeadLoc),Avgcoh(signal2LeadLoc),9,cl(signal2LeadLoc,:),"filled");
end
if ex(4)==1
    scatter(sampleLoc(NegSynchLoc),Avgcoh(NegSynchLoc),9,cl(NegSynchLoc,:),"filled");
end
xlabel('Time (samples)')
ylabel('Coherence')
title('Coherence Along Time')
legend(SynchCondition(existValue))
hold off

%pie plot
%figure;
subplot(2,2,4)
percen=table2array(datatable(5,:));
p=pie(percen);
if ex(1)==1
    p(1).FaceColor=colors(1,:);
end
if ex(2)==1 && ex(1)==1
    p(3).FaceColor=colors(2,:);
elseif ex(2)==1
    p(1).FaceColor=colors(2,:);
end
if ex(3)==1 && ex(2)==1 && ex(1)==1
    p(5).FaceColor=colors(3,:);
elseif ex(3)==1 && ((ex(2)==1 && ex(1)==0) || (ex(1)==1 && ex(2)==0))
    p(3).FaceColor=colors(3,:);
elseif ex(3)==1 && (ex(2)==0 && ex(1)==0)
    p(1).FaceColor=colors(3,:);
end
if ex(4)==1 && ex(3)==1 && ex(2)==1 && ex(1)==1
    p(7).FaceColor=colors(4,:);
elseif ex(4)==1 && ((ex(3)==1 && ex(2)==1 && ex(1)==0)||(ex(3)==1 && ex(2)==0 && ex(1)==1)...
        ||(ex(3)==0 && ex(2)==1 && ex(1)==1))
    p(5).FaceColor=colors(4,:);
elseif ex(4)==1 && ((ex(3)==1 && ex(2)==0 && ex(1)==0)||(ex(3)==0 && ex(2)==1 && ex(1)==0)...
        ||(ex(3)==0 && ex(2)==0 && ex(1)==1))
    p(3).FaceColor=colors(4,:);
elseif ex(4)==1 && ((ex(3)==0 && ex(2)==0 && ex(1)==0))
    p(1).FaceColor=colors(4,:);
end
title('Time Percentage')
legend(SynchCondition(existValue),'Location','northeastoutside','Orientation','vertical')
hold off

writetable(datatable,'datatable.xlsx','WriteRowNames',true);
end




