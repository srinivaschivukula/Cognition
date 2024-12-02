function Out=Dec_SubspaceProj(AllTrialData,Phase,unit,timeWindow,opts)
%%
% timeWindow=[-1.5 -.5];
% Task related Firing Rate
disp(timeWindow)
%%
LabelOrder={'Pinch'    'Press'    'Rub'    'Tap'};
cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase);

%%

% Conditions are basically {{Fc,Fs,Oc},{Fc,Fs,Os},{Fc,Oc,Os},{Fs,Oc,Os}}

Fc = {'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek'};
Fs = {'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder'};
Oc = {'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek'};
Os = {'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'};

ConditinListNames = {{'Fc','Fs','Oc'},...
    {'Fc','Fs','Os'},...
    {'Fc','Oc','Os'},...
    {'Fs','Oc','Os'}};

ConditinList={{Fc,Fs,Oc},...
    {Fc,Fs,Os},...
    {Fc,Oc,Os},...
    {Fs,Oc,Os}};

for WrapIDX=1:length(ConditinList)
    Conditions=ConditinList{WrapIDX};
    ConditionsName = ConditinListNames{WrapIDX};
    wrapName=cat(2,ConditionsName{:});
    AllMDL=[];
    
    %% Within Class Stats
    for i=1:length(Conditions)
        cCond=Conditions{i};
        cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond);
        
        
        L{i}=Analyze.returnFieldValues(cTrialData,'Action');
        FR{i}=squeeze(Analyze.getNeuralData(cTrialData, unit,timeWindow));
        
        [FRCell{i},LabelCell{i}]=DataConvert.list2cell(FR{i},L{i});
        [FRTrue{i},LabelTrue{i}]=DataConvert.cell2list(FRCell{i});
        
    end
    
    FormatAll=LabelTrue; for i=1:length(FormatAll); FormatAll{i}=FormatAll{i}*0+i; end
    
    FRBaseLine=squeeze(Analyze.getNeuralData(Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition','Null'), unit,[.5 2.5]));
    
    %%
    B=FRBaseLine(:,1)*0;
    
    XAll=[cat(1,FRTrue{:});FRBaseLine];
    YAll=[cat(1,LabelTrue{:});B];
    FAll=[cat(1,FormatAll{:});B];

    
    %%
    
    cmps={[1 2],[1 3],[2 3],[ 1 2 3],1,2,3};
    
    clear gAUC gAUCp gR2 gR2p ErrRedCoef AUC  R2  R2p  RefCond PairComp PairCompP
    pairIDX=1;
    ActClr=prism(4);
    tic
    for cmpIDX=1:length(cmps)
        FT=cmps{cmpIDX};
        
        for actIDX=1:4
            if length(FT)==1
                idx2test=1;
            else
                idx2test=1:4;
            end
            for shuffRep=idx2test
                disp([cmpIDX actIDX shuffRep])
                %
                
                
                Y=zeros(size(YAll));
                % activeIDX= any(YAll==AT,2) & any(FAll==FT,2);
                if length(FT)==3
                    if actIDX~=shuffRep
                        foo=setdiff(idx2test,[actIDX shuffRep]);
                        shuffRep2=foo(randi(length(foo)));
                    else
                        shuffRep2=shuffRep;
                    end
                    AT=[actIDX shuffRep shuffRep2];
                    GroupIDX=FAll==FT & YAll==AT;
                elseif length(FT)==1
                    AT=actIDX;
                    GroupIDX=FAll==FT & YAll==AT;
                else
                    AT=[actIDX shuffRep];
                    GroupIDX=FAll==FT & YAll==AT;
                end
                activeIDX=any(GroupIDX,2);
                
                Y(activeIDX)=1;
                X=XAll;
                
                
                Results=Analyze.FitPLSLM(X,Y,'GroupIDX',GroupIDX);
                
                %%
                
                % Plot what is being fit.
                if 1==1
                    clr={'r.','g.','b.'};
                    figure(16);clf; hold on
                    plot(Results.testTarg,'.','color',[.5 .5 .5],'markersize',15)
                    tmp=Results.testTarg;
                    for gidx=1:size(GroupIDX,2)
                        inds=find(GroupIDX(:,gidx));
                        plot(inds,tmp(inds),'.','markersize',20,'color',ActClr(AT(gidx),:))
                    end
                    axis tight
                    tm=ylim;
                    ylim([tm(1)-.1 tm(2)])
                    clr=lines(4);
                    for ii=1:4
                        xx=find(FAll==ii-1);
                        yy=tm(1)-.1;
                        plot([xx(2) xx(end-1)],[yy yy],'color',clr(ii,:),'linewidth',10)
                    end
                    
                    
                    for ii=1:4
                        xx=find(YAll==ii);
                        yy=tm(1)-.05;
                        %            plot([xx(2) xx(end-1)],[yy yy],'color',clr(ii,:),'linewidth',10)
                        plot(xx,yy,'+','color',ActClr(ii,:),'markersize',10)
                    end
                end
                %%
                [val,idx]=min(Results.gAUC);
                
                gAUC(cmpIDX,actIDX,shuffRep)=val;
                gAUCp(cmpIDX,actIDX,shuffRep)=Results.gAUCp(idx);
                
                [val,idx]=min(Results.gdP);
                gdP(cmpIDX,actIDX,shuffRep)=val;
                gdPp(cmpIDX,actIDX,shuffRep)=Results.gdPp(idx);
                
                
                
                [val,idx]=min(Results.gR2);
                gR2(cmpIDX,actIDX,shuffRep)=val;
                gR2p(cmpIDX,actIDX,shuffRep)=Results.gR2p(idx);
                
                
                ErrRedCoef(cmpIDX,actIDX,shuffRep)=Results.ErrRedCoef;
                AUC(cmpIDX,actIDX,shuffRep)=Results.AUC;
                R2(cmpIDX,actIDX,shuffRep)=Results.CVR2;
                R2p(cmpIDX,actIDX,shuffRep)=Results.CVR2p;
                
                %%
                testTarg=Results.testTarg;
                
                % FT is the vector of format indices included in this subspace
                % AT is the vector of actions each associated with the format from FT 
                for ftIDX=1:length(FT); 
                    S1=FAll==FT(ftIDX) & YAll==AT(ftIDX);
                    RefCond(pairIDX,ftIDX,:)=[FT(ftIDX) AT(ftIDX)];
                    for FT2=1:3
                        for Act2=1:4
                            ID2=Act2+((FT2-1)*4);
                            
                            S2=FAll==FT2 & YAll==Act2;
                            PairComp(pairIDX,ftIDX,ID2)=Analyze.TwoGroup.ShuffleCompare(testTarg(S1),testTarg(S2),'Type','dprime','NShuffles',0);
                            PairCompP(pairIDX,ftIDX,ID2)=ranksum(testTarg(S1),testTarg(S2));
                        end
                    end
                end
                
                %%
                  %% Get mean and std for each class
                    
      for refIDX=1:length(AT)
                       inGroup(refIDX)= AT(refIDX)+((FT(refIDX)-1)*4);
                       inGroupFA(refIDX,:)= [FT(refIDX) AT(refIDX) ]; 
                    end
 
                    testTarg=Results.testTarg;
                    for refFormatIDX=1:3
                        for refActionIDX=1:4
                            Sref=FAll==refFormatIDX & YAll==refActionIDX;
                            ID1=refActionIDX+((refFormatIDX-1)*4);
                            
                            MU(ID1)=mean(testTarg(Sref));
                            [tmp] = studentFit(testTarg(Sref),4);
                            tMU(ID1)=tmp.mu;
                            SIGMA(ID1)=std(testTarg(Sref))';
                        end
                    end
                    
                    % also baaseline
                    Sref=FAll==0 & YAll==0;                    
                    MU(ID1+1)=mean(testTarg(Sref));
                      [tmp] = studentFit(testTarg(Sref));
                            tMU(ID1+1)=tmp.mu;
                    SIGMA(ID1+1)=std(testTarg(Sref))';
                    
                    sigma=mean(SIGMA);
%                     figure(11);cla; plot(MU,'r.'); hold on; plot(tMU,'k.')
                    %%
                    % now compute pairwise distances.
                    for refIDX=1:12
                        for targIDX=1:length(MU)
                        PWC(refIDX,targIDX)=(MU(refIDX)-MU(targIDX))./sqrt((sigma.^2+sigma.^2)/2);
                        tPWC(refIDX,targIDX)=(tMU(refIDX)-tMU(targIDX))./sqrt((sigma.^2+sigma.^2)/2);
                        end
                    end
                    %%
               
                    
                    
                        %%
                                                        
                % FT is the vector of format indices included in this subspace
                % AT is the vector of actions each associated with the format from FT
                for refFormatIDX=1:3
                    for refActionIDX=1:4
                        Sref=FAll==refFormatIDX & YAll==refActionIDX;
                        ID1=refActionIDX+((refFormatIDX-1)*4);
                        
                        for targFormatIDX=1:3
                            for targActionIDX=1:4
                                Starg=FAll==targFormatIDX & YAll==targActionIDX;
                                ID2=targActionIDX+((targFormatIDX-1)*4);
                                
                                PWC2(ID1,ID2)=Analyze.TwoGroup.ShuffleCompare(testTarg(Sref),testTarg(Starg),'Type','dprime','NShuffles',0);
                                PWC2p(ID2,ID2)=ranksum(testTarg(Sref),testTarg(Starg));
                                
                            end
                        end
                    end
                end
                
                  for refFormatIDX=1:3
                    for refActionIDX=1:4
                        Sref=FAll==refFormatIDX & YAll==refActionIDX;
                        ID1=refActionIDX+((refFormatIDX-1)*4);
                        
                                Starg=FAll==0 & YAll==0;
                                ID2=16;
                                
                                PWC2(ID1,ID2)=Analyze.TwoGroup.ShuffleCompare(testTarg(Sref),testTarg(Starg),'Type','dprime','NShuffles',0);
                                PWC2p(ID2,ID2)=ranksum(testTarg(Sref),testTarg(Starg));
                        
                    end
                  end
                
                     PairWiseCom{pairIDX}=PWC;
                      PairWiseCom2{pairIDX}=PWC2;
                      PairWiseCom2p{pairIDX}=PWC2p;
                PairWiseComtDist{pairIDX}=tPWC;
                GroupSigma{pairIDX}=SIGMA;
                GrouptMeans{pairIDX}=tMU;
                GroupMeans{pairIDX}=MU;
                
                     IngroupIDX{pairIDX}=inGroup;
                       IngroupFormatAction{pairIDX}=inGroupFA;
                       FormatGroups{pairIDX}=cmpIDX;
                %%
                pairIDX=pairIDX+1;
                size(PairCompP)
            end
        end
    end
    toc


        Out.([wrapName '_' 'GroupMeans'])=GroupMeans;
        Out.([wrapName '_' 'IngroupIDX'])=IngroupIDX;
    Out.([wrapName '_' 'IngroupFormatAction'])=IngroupFormatAction;
    Out.([wrapName '_' 'FormatGroups'])=FormatGroups;
    
        Out.([wrapName '_' 'PairWiseCom'])=PairWiseCom;
    Out.([wrapName '_' 'PairWiseCom2'])=PairWiseCom2;
    Out.([wrapName '_' 'PairWiseCom2p'])=PairWiseCom2p;
        Out.([wrapName '_' 'PairWiseComtDist'])=PairWiseComtDist;
    Out.([wrapName '_' 'GroupSigma'])=GroupSigma;
    Out.([wrapName '_' 'GrouptMeans'])=GrouptMeans;
    
    
    
    
    Out.([wrapName '_' 'RefCond'])=RefCond;
    Out.([wrapName '_' 'PairComp'])=PairComp;
    Out.([wrapName '_' 'PairCompP'])=PairCompP;
    
    Out.([wrapName '_' 'gAUC'])=gAUC;
    Out.([wrapName '_' 'gAUCp'])=gAUCp;
    
    Out.([wrapName '_' 'gdP'])=gdP;
    Out.([wrapName '_' 'gdPp'])=gdPp;
    
    Out.([wrapName '_' 'gR2'])=gR2;
    Out.([wrapName '_' 'gR2p'])=gR2p;
    
    Out.([wrapName '_' 'ErrRedCoef'])=ErrRedCoef;
    Out.([wrapName '_' 'AUC'])=AUC;
    Out.([wrapName '_' 'R2'])=R2;
    Out.([wrapName '_' 'R2p'])=R2p;
end
%% Examples
if 1==1
    cmps={[1 2],[1 3],[2 3],[ 1 2 3],1,2,3};
    cmps={[2 3]};
    a=1;b=4;
    actCMPS={[a a],[a b],[b a],[b b]}
    
    cmps={[1 2],[1 3],[2 3],[ 1 2 3],1,2,3};
    cmps={[1 3]};
    a=3;b=4;
    actCMPS={[a a],[a b],[b a],[b b]}
    %     actCMPS={5};
    % cmps={[11 2],[1 3],};
    clear gAUC gAUCp gR2 gR2p
    clr= [255 102 204 ; 255 0 0; 153 0 204 ; 253 153 0 ; 0 102 255  ; 0 0 0 ]/255;
    
    ActClr=clr;
    tic
    for cmpIDX=1:length(cmps)
        plt.fig('units','inches','width',3,'height',10,'font','Arial','fontsize',9);
        pnl = panel();  pnl.margin=15; pnl.pack('v',4); pnl.fontname='arial';pnl.fontsize=12;
        
        for idx=1:length(actCMPS)
            %%
            FT=cmps{cmpIDX};
            
            
            Y=zeros(size(YAll));
            % activeIDX= any(YAll==AT,2) & any(FAll==FT,2);
            if length(FT)==3
                AT=actCMPS{idx};
                GroupIDX=FAll==FT & YAll==AT;
            elseif length(FT)==1
                AT=actCMPS{idx};
                GroupIDX=FAll==FT & YAll==AT;
            else
                AT=actCMPS{idx};
                GroupIDX=FAll==FT & YAll==AT;
            end
            activeIDX=any(GroupIDX,2);
            
            Y(activeIDX)=1;
            X=XAll;
            
            
            Results=Analyze.FitPLSLM(X,Y,'GroupIDX',GroupIDX);
            
            
            
            clr={'r.','g.','b.'};
            pnl(idx).select(); hold on
            plot(Results.testTarg,'.','color',[.5 .5 .5],'markersize',15)
            tmp=Results.testTarg;
            for gidx=1:size(GroupIDX,2)
                inds=find(GroupIDX(:,gidx));
                plot(inds,tmp(inds),'.','markersize',20,'color',ActClr(AT(gidx),:))
                text(mean(inds),mean(tmp(inds)), ...
                    sprintf('AUC/R2=%0.2f/%0.2f',Results.gAUC(gidx),Results.gR2(gidx)))
            end
            axis tight
            tm=ylim;
            ylim([tm(1)-.1 tm(2)])
            clr=lines(4);
            for ii=1:4
                xx=find(FAll==ii-1);
                yy=tm(1)-.1;
                plot([xx(2) xx(end-1)],[yy yy],'color',clr(ii,:),'linewidth',10)
            end
            
            
            for ii=1:4
                xx=find(YAll==ii);
                yy=tm(1)-.05;
                %            plot([xx(2) xx(end-1)],[yy yy],'color',clr(ii,:),'linewidth',10)
                plot(xx,yy,'+','color',ActClr(ii,:),'markersize',10)
            end
            
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            xlabel('Trial')
            ylabel('Subspace Response')
            
        end
    end
end
end
%%