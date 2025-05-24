clear all;clc
cd('E:\1.Scn2a_mDAN\L-DOPA\optogenetic-fiberphotometry\slow change with LB in TH-cre');
cwd=pwd;

Foldidx=dir('fiber-opto-*');
%% 
ipTime=1200
TimeStart=ipTime-300; %Train delay
TimeEnd=ipTime+1200;

TimeLaser=TimeStart:Interval:TimeEnd;

for kk=1:length(Foldidx)
    cd(Foldidx(kk).name)
    Filename=Foldidx(kk).name;

    fidx=[dir('*405*.mat')
    dir('*580*.mat')];

    for k=1:length(fidx)
        load(fidx(k).name)
        tspan=data(1:end-1,1);
        data=data(2:end,3);
        
        dataAll{k}.tspan=tspan;
        dataAll{k}.data=data;
    end

 %%
    Fun=@(Y0,Plateau,KFast,KSlow,PercentFast,x) Plateau + (Y0-Plateau)*PercentFast*.01*exp(-KFast*x)...
        + (Y0-Plateau)*(100-PercentFast)*.01*exp(-KSlow*x);
    fo= fitoptions('Method','NonlinearLeastSquares','MaxIter',10000,...
        'Lower',[5.0, 0,   1/10000,   1/10000,0],...
        'Upper',[200,1000, 1/100, 1/1000,100,], ...
        'Startpoint',[50,20,1/1000,1/10000,50]);
    foNew= fitoptions('Method','NonlinearLeastSquares','MaxIter',10000,...
        'Lower',[5.0, 0],...
        'Upper',[200,1000], ...
        'Startpoint',[50,20]);

    
    yy={};
    for k=1:length(dataAll)
        nanidx=[];
        pretime=10;
        postime=30;
        for i=1:length(TimeLaser)
            tspan=dataAll{k}.tspan;
            idx=find(tspan>TimeLaser(i)-pretime&tspan<TimeLaser(i)+postime);
            nanidx=[nanidx;idx(:)];
        end
        if k==1
            fitidx=tspan>59.5; %可调整：设置曲线参数区段
        else
            fitidx=(tspan>59.5&tspan<TimeLaser(3)|tspan>5609.5); %可调整：设置曲线参数对应区段
        end
        x=tspan(fitidx);x=x-x(1);

        y=dataAll{k}.data;
        y(nanidx)=nan;
        y=naninterp(y);
        y=y(fitidx);
        if k==1
            [model,~]=fit(x(:),y(:),Fun,fo);
        else
            param=coeffvalues(model);param=param(3:end);
            FunNew=@(Y0,Plateau,x) Plateau + (Y0-Plateau)*param(3)*.01*exp(-param(1)*x)...
                + (Y0-Plateau)*(100-param(3))*.01*exp(-param(2)*x);
            [model,~]=fit(x(:),y(:),FunNew,foNew);
        end
        
        yy{k}.data=feval(model,tspan);
        yy{k}.tspan=tspan;
   

    end

    dataSM=dataAll{2}.data-yy{2}.data;   %+yy{2}.data(1)
    dataSM(nanidx)=nan;
    dataSM=naninterp(dataSM);
    dataSM=dataSM(659.5*40:6120); %可调整：AUC与amplitude分析范围
    peakidx=find(dataSM==max(dataSM));
    Amplitude=mean(dataSM(peakidx-10:peakidx+10));
    Area=trapz(dataSM);


%    baseline=median(dataSM(60:720)); %可调整：baseline参考范围 
%     EdgePoint=[];
%     Threshold=[];
%     for j=1:2
%         Threshold(j,1)=Amplitude*(0.20*(j==1)+0.8*(j==2))+baseline;
%         %temp=eventedgefinder(dataSM,Threshold(j),1,0)*dt;
%         temp=eventedgefinder(dataSM,Threshold(j),1,0)*0.025;
%         EdgePoint(j)=temp(end,2);
%     end
%     HW=EdgePoint(1)-EdgePoint(2);
    
    paramWhole=[Amplitude,Area]; %,HW
    paramWholeName={'Amplitude','Area'}; %,'20%~80% decay'
    FilenameRecord{kk,1}=Filename;
    paramWholeRecord(kk,:)=paramWhole;

%%
    figure(1),clf
    Dist=40;
    for k=1:length(dataAll)
          plot(dataAll{k}.tspan,dataAll{k}.data+(k-1)*Dist),hold on  %raw data
          plot(yy{k}.tspan,yy{k}.data+(k-1)*Dist,'k','linewidth',3),hold on
          plot(yy{k}.tspan,dataAll{k}.data-yy{k}.data+(k-1)*Dist+yy{k}.data(1),'r') %correction base
    end
    plot(TimeLaser,zeros(size(TimeLaser)),'ro')
    
    
    %%
    cd(cwd)
    
end   

xlswrite('Results.xlsx',FilenameRecord,'paramWhole','A2')
xlswrite('Results.xlsx',paramWholeName,'paramWhole','B1')
xlswrite('Results.xlsx',paramWholeRecord,'paramWhole','B2')
