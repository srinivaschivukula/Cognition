

%%

[varargin,rescale]   = Utilities.ProcVarargin(pltFCNArgs,'rescale',0);
[varargin,Sessions]   = Utilities.ProcVarargin(pltFCNArgs,'Sessions',[]);

% typeID={'DecodeTest_corrNM'};
%
% cvAccuracy=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,typeID{decodeType}));
%
% FigHandels(1)=plt.fig('units','inches','width',8,'height',8,'font','Arial','fontsize',12);
%
%  imagesc(mean(cvAccuracy,3))
%
% plt.vline([.5:5:20.5],{'k','linewidth',2}); plt.hline([.5:5:20.5],{'k','linewidth',2});
% colormap(jet)
% colorbar

tickLoc=7; StartInds=1;
%%

F=Analyze.returnFieldValues(PlotData,'DecodeTest_FormatCorrs');
FOut=[];
for sessionIDX=1:length(F)
    for format1=1:4
        for format2=1:4
            if sessionIDX==1
                FOut{format1,format2}=[];
            end
            FOut{format1,format2}=cat(3,FOut{format1,format2},F{sessionIDX}{format1,format2}(StartInds:end,StartInds:end,:));
        end
    end
end

timeWindow=opts.timeWindow(StartInds:end);
%%
for format1=1:4
    for format2=1:4
        if isempty(Sessions)
            FOutMu{format1,format2}=mean(FOut{format1,format2},3);
        else
            FOutMu{format1,format2}=mean(FOut{format1,format2}(:,:,Sessions),3);
        end
    end
end

%%
for format1=1:4
    for format2=1:4
        if isempty(Sessions)
            for sesID=1:size(FOut{format1,format2},3)
            FOutDiag{format1,format2}(:,sesID)=diag(FOut{format1,format2}(:,:,sesID));
            end
        else
            FOutDiag{format1,format2}=mean(FOut{format1,format2}(:,:,Sessions),3);
        end
    end
end
% for format1=1:4
%     for format2=1:4
%       FOutMu{format1,format2}=mean(FOut{format1,format2}(:,:,8:13),3)-mean(FOut{format1,format2}(:,:,1:7),3);
%     end
% end
%%
% plt.fig('units','inches','width',5,'height',4,'font','Arial','fontsize',12);
% hold on
% 
% F1=3;
% 
% textResp=FOutDiag{1,1}';v=mean(textResp,1); textResp=(textResp-v(1))/(max(v-v(1)));
% 
% 
% visResp=[FOutDiag{F1,F1}'];v=mean(visResp,1); visResp=(visResp-v(1))/(max(v-v(1)));
% 
% 
% 
% cross=[FOutDiag{1,F1}'];v=mean(cross,1); cross=(cross-v(1))/(max(v-v(1)));
% Analyze.plotEventRelatedAverage({textResp,visResp,cross},{'Text','Visual','CrossModal'},'TimeVec',timeWindow-.05,'legend')
% 
% xlim([-.1 1]);
% ylabel('Normalized Correlation')
% xlabel('Time')
% title('')

%%

plt.fig('units','inches','width',5,'height',4,'font','Arial','fontsize',12);
hold on

clear rescale;

allVals=cat(1,FOutMu{:}); mx=max(allVals(:)); mn=min(allVals(:));
pk95=prctile(allVals(:),[5 95]);

% F1=3;
tmp = [F1;F2;F3;F4;F5;F6];

F1=FOutDiag{2,1}';F1 = rescale(F1,pk95(1),pk95(2));
F2=FOutDiag{3,1}';F2 = rescale(F2,pk95(1),pk95(2));
F3=FOutDiag{4,1}';F3 = rescale(F3,pk95(1),pk95(2));
F4=FOutDiag{3,2}';F4 = rescale(F4,pk95(1),pk95(2));
F5=FOutDiag{4,2}';F5 = rescale(F5,pk95(1),pk95(2));
F6=FOutDiag{4,3}';F6 = rescale(F6,pk95(1),pk95(2));

F1=FOutDiag{2,1}';v=mean(F1,1); F1=(F1-v(1))/(max(v-v(1)));
F2=FOutDiag{3,1}';v=mean(F2,1); F2=(F2-v(1))/(max(v-v(1)));
F3=FOutDiag{4,1}';v=mean(F3,1); F3=(F3-v(1))/(max(v-v(1)));
F4=FOutDiag{3,2}';v=mean(F4,1); F4=(F4-v(1))/(max(v-v(1)));
F5=FOutDiag{4,2}';v=mean(F5,1); F5=(F5-v(1))/(max(v-v(1)));
F6=FOutDiag{4,3}';v=mean(F6,1); F6=(F6-v(1))/(max(v-v(1)));
% 
F1=FOutDiag{2,1}';v=mean(tmp,1); F1=(F1-v(1))/(max(v-v(1)));
F2=FOutDiag{3,1}';v=mean(tmp,1); F2=(F2-v(1))/(max(v-v(1)));
F3=FOutDiag{4,1}';v=mean(tmp,1); F3=(F3-v(1))/(max(v-v(1)));
F4=FOutDiag{3,2}';v=mean(tmp,1); F4=(F4-v(1))/(max(v-v(1)));
F5=FOutDiag{4,2}';v=mean(tmp,1); F5=(F5-v(1))/(max(v-v(1)));
F6=FOutDiag{4,3}';v=mean(tmp,1); F6=(F6-v(1))/(max(v-v(1)));

Analyze.plotEventRelatedAverage({F1,F2,F3,F4,F5,F6},{'Fc','Fs','Oc','Os','Oc','Os'},'TimeVec',timeWindow-.05,'legend')

legend('Fc-Fs','Fc-Oc','Fc-Os','Fs-Oc','Fs-Os','Oc-Os')

xlim([-.6 4]);
ylabel('Normalized Correlation')
xlabel('Time')
title('')
%%
plt.fig('units','inches','width',5,'height',4,'font','Arial','fontsize',12);
; hold on

F1=3;

textResp=FOutDiag{1,1}';v=mean(textResp,1); textResp=(textResp-v(1))/(max(v-v(1)));


visResp=[FOutDiag{2,2}';FOutDiag{3,3}';FOutDiag{4,4}'];
v=mean(visResp,1); visResp=(visResp-v(1))/(max(v-v(1)));



cross=[FOutDiag{1,2}';FOutDiag{1,3}';FOutDiag{1,4}'];

v=mean(cross,1); cross=(cross-v(1))/(max(v-v(1)));
Analyze.plotEventRelatedAverage({textResp,visResp,cross},{'Text','Visual','CrossModal'},'TimeVec',timeWindow-.05,'legend')

xlim([-.1 1]);
ylabel('Normalized Correlation')
xlabel('Time')
title('')

%%
% v=diag(FOutMu{1,1}); v=v-v(1);v=v/max(v);
%         plot(v,'k')
% v=diag(FOutMu{3,3});v=v-v(1); v=v/max(v);
%         plot(v,'b')
%             v=diag(FOutMu{1,3});v=v-v(1); v=v/max(v);
%         plot(v,'r')
