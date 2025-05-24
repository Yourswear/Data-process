clear all;clc

%Fidx{1}=dir('DATA_FIcurve_one-week*.mat');
Fidx{1}=dir('DATA_FIcurve_*.mat');
%Fidx{2}=dir('DATA_FIcurve_three-week*.mat');
GroupName={'DAN'};

Num=1;%选择大部分细胞都会达到的一个AP Num
thAP=1;
colorma=distinguishable_colors(length(Fidx));

%Binwidth=10;
Binwidth=50;
BinL=0:Binwidth:200-Binwidth;
BinR=BinL+Binwidth;
BinsFIcurveS=[BinL',BinR'];

Binwidth=50;
BinL=200:Binwidth:600-Binwidth;
BinR=BinL+Binwidth;
BinsFIcurveM=[BinL',BinR'];

%Binwidth=100;
Binwidth=50;
BinL=600:Binwidth:1500-Binwidth;
BinR=BinL+Binwidth;
BinsFIcurveL=[BinL',BinR'];

BinsFIcurve{1}=[BinsFIcurveS;BinsFIcurveM;BinsFIcurveL];
%%

VposiRecord=[];

for kk=1:length(Fidx)
    for k=1:length(Fidx{kk})
        filename=Fidx{kk}(k).name;
        load(filename)
        dt=tspanV(2)-tspanV(1);
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',num2str(kk),'/',num2str(length(Fidx)),'--',filename])
        FilenameRecord{kk}{k,1}=filename;
        FilenameRecord_row{kk}{1,k}=filename;
        %% fitting ReLU
        Fun=@(a,b,x) max([zeros(length(x),1),a.*x(:)+b],[],2);
        fo= fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0.001,-1000],...
            'Upper',[10,1000], ...
            'Startpoint',[0.1,-10]);
        
        uniqueStim=unique(StimAmp);
        StimAmp_act=zeros(size(uniqueStim))';
        spkFreq_act=zeros(size(uniqueStim))';
        for i=1:length(uniqueStim)
            idx=find(StimAmp==uniqueStim(i));
            StimAmp_act(i)=mean(StimAmp(idx));
            spkFreq_act(i)=mean(spkFreq(idx));
        end
        
        [tempx,idx]=sort(StimAmp_act);tempy=spkFreq_act(idx);
        idx=find(tempy>=20,1,'first');
        x=[0;tempx(1:idx)];
        y=[0;tempy(1:idx)];
        
%        [Model,~]=fit(x(:),y(:),fittype(Fun),fo);
%        xx=linspace(0,x(end),100);
%        yy=feval(Model,xx);
        
%        FIcurveParamRecord{kk}(k,:)=[-(Model.b/Model.a),Model.a];
        %% ------------- Bin FIcurve ---------------------------
        
        for i=1:size(BinsFIcurve{kk},1)
            idx=find(StimAmp>=BinsFIcurve{kk}(i,1)&StimAmp<BinsFIcurve{kk}(i,2));
            if ~isempty(idx)
                spkFreqBin{kk}(k,i)=mean(spkFreq(idx));
            else
                spkFreqBin{kk}(k,i)=nan;
            end
        end
        
        figure(1),clf
        plot(StimAmp_act,spkFreq_act,'o','markeredgecolor',colorma(kk,:),'markersize',6),hold on
%        plot(xx,yy,'m')
        drawnow
%         pause()        
        %% ------------- Bin spkshape,1~5 AP first spike ---------------------------
        
        Wantidx=find(ismember(spkFreq,Num*2)&(StimAmp(:)<600));
        
        if ~isempty(Wantidx)
            spkwaveform_temp=[];
            HWtemp=[];
            AHPtemp=[];
            Thresholdtemp=[];
            Amptemp=[];
            max_slopetemp=[];
            min_slopetemp=[];
            
            for jj=1:length(Wantidx)
                data=waveformV(:,Wantidx(jj));
                [spkwaveform,ttspanspk,HW, AHP, Threshold, Amp,max_slope,min_slope,...
                    HWidx,AHPidx, Thresholdidx,max_slopeidx,min_slopeidx] = AP_Statistic(data,spkTime{Wantidx(jj)},tspanV,1:min(Num));
                
                spkwaveform_temp(jj,:,:)=spkwaveform;
                HWtemp(:,jj)=HW;
                AHPtemp(:,jj)=AHP;
                Thresholdtemp(:,jj)=Threshold;
                Amptemp(:,jj)=Amp;
                max_slopetemp(:,jj)=max_slope;
                min_slopetemp(:,jj)=min_slope;
                
                figure(100),clf
                H=[];
                for j=1:size(spkwaveform,2)
                    H(j)=subplot(1,size(spkwaveform,2),j);
                    plot(ttspanspk,spkwaveform(:,j)),hold on
                    plot(ttspanspk(HWidx(j,:)),spkwaveform(HWidx(j,:),j),'ro')
                    plot(ttspanspk(AHPidx(j)),spkwaveform(AHPidx(j),j),'ko')
                    plot(ttspanspk(Thresholdidx(j)),spkwaveform(Thresholdidx(j),j),'co')
                    axis tight
                    drawnow
                end
                linkaxes(H,'x')
                 pause(0.1)
            end
            spkparamRecord{kk}{1}(:,k)=mean(HWtemp,2);
            spkparamRecord{kk}{2}(:,k)=mean(AHPtemp,2);
            spkparamRecord{kk}{3}(:,k)=mean(Thresholdtemp,2);
            spkparamRecord{kk}{4}(:,k)=mean(Amptemp,2);
            spkparamRecord{kk}{5}(:,k)=mean(max_slopetemp,2);
            spkparamRecord{kk}{6}(:,k)=mean(min_slopetemp,2);
            
            if size(spkwaveform_temp,2)==2600
                spkwaveformRecord{kk}(k,:,:)=squeeze(mean(spkwaveform_temp,1));
            else
                spkwaveformRecord{kk}(k,:,:)=nan(2600,min(Num));
            end
        else
            spkparamRecord{kk}{1}(:,k)=nan(min(Num),1);
            spkparamRecord{kk}{2}(:,k)=nan(min(Num),1);
            spkparamRecord{kk}{3}(:,k)=nan(min(Num),1);
            spkparamRecord{kk}{4}(:,k)=nan(min(Num),1);
            spkparamRecord{kk}{5}(:,k)=nan(min(Num),1);
            spkparamRecord{kk}{6}(:,k)=nan(min(Num),1);
            spkwaveformRecord{kk}(k,:,:)=nan(2600,min(Num));
        end
        
        %         %% Vposi,subthreshold
%         uniqueStimAmpSub=unique(StimAmp(~isnan(Vposi)));
%         if ~isempty(uniqueStimAmpSub)
%             for i=1:length(uniqueStimAmpSub)
%                 idx=StimAmp==uniqueStimAmpSub(i)&~isnan(Vposi);
%                 
%                 VposiRecord=[VposiRecord;kk,mean(StimAmp(idx)),mean(Vposi(idx))];
%                 
%             end
%         end
        %
%         Goodidx=[50,100,200,300];
%         for i=1:length(Goodidx)
%             idx=StimAmp==Goodidx(i)&~isnan(Vposi);
%             if ~isempty(idx)&size(waveformV(:,1),1)==45001;
%                 SubwaveformRecord{kk}{i}(:,k)=mean(waveformV(:,idx)-repmat(mean(waveformV(tspanV<0,idx)),length(tspanV),1),2);
%             else
%                 SubwaveformRecord{kk}{i}(:,k)=nan(45001,1);
%             end
%         end
        
    end
end

FIcurveParamName={'theta','slope'};
spkParamName={'HW', 'AHP', 'Threshold', 'Amp','max_slope','min_slope'};

save('Results_FIcurve.mat','Fidx','GroupName','colorma',...
    'tspanV',...
    'VposiRecord',...
    'spkwaveformRecord','spkparamRecord','ttspanspk',...
    'spkFreqBin','BinsFIcurve',...
    'Num','thAP',...
    'FIcurveParamName','spkParamName','GroupName')

for kk=1:length(Fidx)

    xlswrite('Results_FIcurve.xls',FilenameRecord_row{kk},['FIcurve_',GroupName{kk}],'B1')
    xlswrite('Results_FIcurve.xls',[BinsFIcurve{kk}(:,1),spkFreqBin{kk}'],['FIcurve_',GroupName{kk}],'A2')
    
    xlswrite('Results_FIcurve.xls',FilenameRecord{kk},['FIcurveParam_',GroupName{kk}],'A2')
    xlswrite('Results_FIcurve.xls',FIcurveParamName,['FIcurveParam_',GroupName{kk}],'B1')
    %xlswrite('Results_FIcurve.xls',FIcurveParamRecord{kk},['FIcurveParam_',GroupName{kk}],'B2')
    
    xlswrite('Results_FIcurve.xls',FilenameRecord{kk},['spkParam',GroupName{kk}],'A2')
    xlswrite('Results_FIcurve.xls',spkParamName,['spkParam_',GroupName{kk}],'B1')
    xlswrite('Results_FIcurve.xls',[spkparamRecord{kk}{1}(thAP,:)',...
                                             spkparamRecord{kk}{2}(thAP,:)',...
                                             spkparamRecord{kk}{3}(thAP,:)',...
                                             spkparamRecord{kk}{4}(thAP,:)',...
                                             spkparamRecord{kk}{5}(thAP,:)',...
                                             spkparamRecord{kk}{6}(thAP,:)'],['spkParam_',GroupName{kk}],'B2')
    
end

%%

clear all;clc
load('Results_FIcurve.mat')

% ==================== FIcurve ====================
figure(1),clf
for kk=1:length(Fidx)
    
    x=BinsFIcurve{kk}(:,1);
    y=nanmean(spkFreqBin{kk},1);
    N=sum(~isnan(spkFreqBin{kk}),1);
    sem=nanstd(spkFreqBin{kk},[],1)./sqrt(N);
    SD=nanstd(spkFreqBin{kk},[],1);
    VAR=var(spkFreqBin{kk},[],1);
    
    h3=subplot(313);
    plot(x,y./SD,'ro','color',colorma(kk,:),'markerfacecolor',colorma(kk,:),'markersize',8),hold on
    
    h2=subplot(312);
    plot(x,SD,'ro','color',colorma(kk,:),'markerfacecolor',colorma(kk,:),'markersize',8),hold on

    h1=subplot(311);
    errorbar(x,y,sem,'ro','color',colorma(kk,:),'markerfacecolor',colorma(kk,:),'markersize',8),hold on
    %     plot(x,y,'ro','color',colorma(kk,:),'markerfacecolor','w','markersize',8,'linestyle','none'),hold on

    idx=find(~isnan(y));
    x=[0;x(idx)];
    y=[0,y(idx)];
    xx=linspace(0,x(end),100);
    
    %     ====== log-normal distribution fit
    FunLogn=@(muA,muB,stdA,stdB,A,x) A.*logncdf(x,muA,stdA).*(1-logncdf(x,muB,stdB));
    foLogn= fitoptions('Method','NonlinearLeastSquares','MaxIter',10000,'TolX',1e-20,'TolFun',1e-20,...
        'Lower',[log(1),log(200),0,0,2],...
        'Upper',[log(2000),log(10000),log(1000),log(1000),1000], ...
        'Startpoint',[log(300),log(1300),log(100),log(1000),50]);

    [ModelLogn,~]=fit(x(:),y(:),fittype(FunLogn),foLogn)
    yy=feval(ModelLogn,xx);
    plot(xx,yy,'-','color',colorma(kk,:))
    
end
ylabel('Frequency (Hz)')
xlabel('Current injection (pA)')
linkaxes([h1,h2,h3],'x')
% 
% for i=1:size(spkFreqBin{1},2)
%     A=spkFreqBin{1}(:,i);A=A(~isnan(A));
%     B=spkFreqBin{2}(:,i);B=B(~isnan(B));
%     if ~isempty(A)&~isempty(B)
% %         [~,p] = ttest2(A,B);
%         p = ranksum(A,B);
% 
%         figure(1)
%         subplot(311);
%         if p<0.05&p>=0.01;
%             text(BinsFIcurve{1}(i,1),y(i),'*','fontsize',20);
%         elseif p<0.01&p>=0.001;
%             text(BinsFIcurve{1}(i,1),y(i),'**','fontsize',20);
%         elseif p<0.001
%             text(BinsFIcurve{1}(i,1),y(i),'***','fontsize',20);
%         end
%     end
% end

%

% figure(20),clf
% for kk=1:length(Fidx)
%     
%     x=mean(BinsFIcurve{kk},2);
%     y=nanmean(spkFreqBin{kk},1);
%     N=sum(~isnan(spkFreqBin{kk}),1);
%     sem=nanstd(spkFreqBin{kk},[],1)./sqrt(N);
%     SD=nanstd(spkFreqBin{kk},[],1);
%     
%     idx=find(~isnan(y));
%     x=x(idx);
%     y=y(idx);
%     N=N(idx);
%     sem=sem(idx);
%     SD=SD(idx);
%     
%     %     ====== log-normal distribution fit
%     FunLogn=@(muA,stdA,A,x) A.*logncdf(x,muA,stdA);
%     foLogn= fitoptions('Method','NonlinearLeastSquares','MaxIter',10000,'TolX',1e-20,'TolFun',1e-20,...
%         'Lower',[log(1),0,2],...
%         'Upper',[log(2000),log(1000),1000], ...
%         'Startpoint',[log(300),log(100),50]);
%     
%     Maxidx=find(y==max(y));
%     ModelLogn=[];
%     Record=[];
%     for j=1:5
%         x_act=x(1:Maxidx-j);
%         y_act=y(1:Maxidx-j);
%         [ModelLogn{j},GOF]=fit(x_act(:),y_act(:),fittype(FunLogn),foLogn);
%         Record(j)=GOF.rsquare;
%     end
%     idx=find(Record==max(Record),1,'first');
%      
%     x=x(1:Maxidx-idx);
%     y=y(1:Maxidx-idx);
%     N=N(1:Maxidx-idx);
%     sem=sem(1:Maxidx-idx);
%     SD=SD(1:Maxidx-idx);
%     xx=linspace(0,x(end),100);
%     
%     figure(20)
%     subplot(2,length(Fidx),kk);
% %     plot(x,spkFreqBin{kk}(:,1:Maxidx-idx),'k-o','markerfacecolor','w'),hold on
%     errorbar(x,y,sem,'ro','color',colorma(kk,:),'markerfacecolor',colorma(kk,:),'markersize',8),hold on
%     yy=feval(ModelLogn{idx},xx);
%     plot(xx,yy)
%     xlabel('pA')
%     ylabel('Hz')
%     title(GroupName{kk})
% 
%     subplot(2,length(Fidx),kk+length(Fidx));
%     plot(x,y./SD,'ro','color',colorma(kk,:),'markerfacecolor',colorma(kk,:),'markersize',8),hold on
%     x=x(SD>0);
%     y=y(SD>0);
%     SD=SD(SD>0);
%     
%     Modelpoly=polyfit(x',y./SD,1);
%     R=corrcoef(x',y./SD);
%     yy=polyval(Modelpoly,xx);
%     plot(xx,yy,'color',colorma(kk,:))
%     xlabel('pA')
%     ylabel('Mean/SD')
%     title({['k = '  num2str(Modelpoly(1))],['b = '  num2str(Modelpoly(2))],['R = '  num2str(R(2))]})
% 
% end

% ====================  statistic ====================
figure(4),clf
figure(400),clf
for k=1:length(FIcurveParamName)
    
    for kk=1:length(Fidx)
        data{kk}=FIcurveParamRecord{kk}(:,k);N=sum(~isnan(data{kk}));
        Mean=nanmean(data{kk});
        Sem=nanstd(data{kk})./sqrt(N);
        figure(4)
        subplot(1,length(FIcurveParamName),k)
        bar(kk,Mean,'facecolor',colorma(kk,:));hold on
        errorbar(kk,Mean,Sem,'color',colorma(kk,:))
        
        if k==1
        figure(400)
        subplot(length(FIcurveParamName),length(Fidx),kk)
        hist(data{kk},10)
        end
    end
    ylabel(FIcurveParamName{k})
    % ========== test
%     A=data{1};B=data{end};
%     if swtest(A,0.05)==0&swtest(B,0.05)==0
%         [~,p]=ttest2(A,B);
%         tmethod='two samples T-test';
%     else
%         [p,~]=ranksum(A,B);
%         tmethod='Wilcoxon rank sum ';% or Mann-Whitney U-test
%     end
%     title({['No ',num2str(nanmean(A)),'\pm',num2str(nanstd(A)./sqrt(length(A)))]
%         ['YES ',num2str(nanmean(B)),'\pm',num2str(nanstd(B)./sqrt(length(B)))]
%         [tmethod]
%         ['  p = ',num2str(p)]})
    % ==========
end

figure(5),clf
for k=1:length(spkParamName)
    
    subplot(1,length(spkParamName),k)
    for kk=1:length(Fidx)
        data{kk}=spkparamRecord{kk}{k}(thAP,:);N=sum(~isnan(data{kk}));
        Mean=nanmean(data{kk});
        Sem=nanstd(data{kk})./sqrt(N);
        bar(kk,Mean,'facecolor',colorma(kk,:));hold on
        errorbar(kk,Mean,Sem,'color',colorma(kk,:))
    end
    ylabel(spkParamName{k})
    % ========== test
%     A=data{1};B=data{end};
%     if swtest(A,0.05)==0&swtest(B,0.05)==0
%         [~,p]=ttest2(A,B);
%         tmethod='two samples T-test';
%     else
%         [p,~]=ranksum(A,B);
%         tmethod='Wilcoxon rank sum ';% or Mann-Whitney U-test
%     end
%     title({['No ',num2str(nanmean(A)),'\pm',num2str(nanstd(A)./sqrt(length(A)))]
%         ['YES ',num2str(nanmean(B)),'\pm',num2str(nanstd(B)./sqrt(length(B)))]
%         [tmethod]
%         ['  p = ',num2str(p)]})
    % ==========
end



% ====================  waveform  ====================

% figure(6),clf
% figure(7),clf

% xlim_max=0.008;
% xlim_min=-0.0015;
% 
% for kk=1:length(Fidx)
%     
%     for j=1:size(spkwaveformRecord{kk},3)
%         figure(6)
%         subplot(2,size(spkwaveformRecord{kk},3),j);
%         data=squeeze(nanmean(spkwaveformRecord{kk}(:,:,j),1));
%         plot(ttspanspk,data,'color',colorma(kk,:)),hold on
%         axis([xlim_min,xlim_max,-90,60])
%         
%         subplot(2,size(spkwaveformRecord{kk},3),j+size(spkwaveformRecord{kk},3))
%         tempidx=find(ttspanspk>-0.004&ttspanspk<0.0040);
%         plot(data(tempidx(2:end)),diff(data(tempidx))*100,'color',colorma(kk,:)),hold on
%         axis([-90,60,-400,1800])
%         
%     end
    
%     for j=1:length(Goodidx)
%         figure(7)
%         subplot(1,length(Goodidx),j)
%         data=squeeze(nanmean(SubwaveformRecord{kk}{j},2));
%         plot(tspanV,data,'color',colorma(kk,:)),hold on
%         axis tight
%     end
%     
% end

