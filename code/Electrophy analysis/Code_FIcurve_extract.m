clear all;clc
Fidx{1}=dir('FI*.mat');
%Fidx{2}=dir('LL_2021050813_Thy1-Cas9 P65 E22_Slice2C8_GFP+ Ctrl_15.5-16.1_2nd recording.mat');
%Fidx{3}=dir('three-week*.mat');

%%
for kk=1:length(Fidx)
    for k=1:length(Fidx{kk})
        %% loading data
        filename=Fidx{kk}(k).name;
        load(filename)
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',num2str(kk),'/',num2str(length(Fidx)),'--',filename])
        %% ------------- reading data ---------------------------
        clear(('*Ch31*'));
        varia=who('*Ch*');
        
        for i=1:length(varia)
            if isfield(eval(varia{i}),'units')
                str=eval([varia{i},'.units']);
                
                if ismember('pA',str)
                    Im=eval([varia{i} '.values']);
                    Im=Im-mean(Im(1:10));
                    dtI=eval([varia{i} '.interval']);
                end
                
                if ismember('mV',str)
                    Vm=eval([varia{i} '.values']);
                    dtV=eval([varia{i} '.interval']);
                end
                
            end
        end
        
        clear (varia{:})
        
        %%  finding key time point
        TimeStim = eventedgefinder(Im,5,1/dtI,0.4,1,1,0).*dtI;%(Im,5,1/dtI,0.2,0.5,1,0)其中第4个参数可控制刺激长度
        Dura=round(range(TimeStim,2)*100)./100;
        TimeStim=TimeStim(:,1);
        disp(length(TimeStim))
        %% waveform extract
        pretime=0.1;
        postime=0.8;
        
        waveformV=[];
        waveformI=[];
        Vrest=[];
        for i=1:length(TimeStim)
            dataidxV=round(TimeStim(i)/dtV-pretime/dtV:TimeStim(i)/dtV+postime/dtV);
            baseidxV=round(TimeStim(i)/dtV-pretime/dtV:TimeStim(i)/dtV);
            dataidxI=round(TimeStim(i)/dtI-pretime/dtI:TimeStim(i)/dtI+postime/dtI);
            baseidxI=round(TimeStim(i)/dtI-pretime/dtI:TimeStim(i)/dtI);
            
            Vrest=[Vrest,mean(Vm(baseidxV))];
            waveformV=[waveformV,Vm(dataidxV)];
            waveformI=[waveformI,Im(dataidxI)-mean(Im(baseidxI))];
        end
        tspanV=linspace(-pretime,postime,size(waveformV,1));
        tspanI=linspace(-pretime,postime,size(waveformI,1));
        
        wrongidx=[];
%         temp=mean(waveformV(tspanV<0|(tspanV>0.6&tspanV<1),:),1);
%         wrongidx=[wrongidx;find(temp>-45)'];
%         
%         temp=min(waveformV(tspanV<0|tspanV>0.5,:),[],1);
%         wrongidx=[wrongidx;find(temp<-100)'];
%         
%         temp=mean(waveformV(tspanV<0,:),1);
%         wrongidx=[wrongidx;find(temp>-55)'];
        
        
        waveformI(:,wrongidx)=[];
        waveformV(:,wrongidx)=[];
        Vrest(wrongidx)=[];
        TimeStim(wrongidx)=[];
        Dura(wrongidx)=[];
        
        
        %% make same necessary Tag
        % spk shape
        spkFreq=[];
        spkTime=[];
        
        for i=1:size(waveformV,2)
            peaktime=peakfinder(waveformV(:,i)-Vrest(i),30,0, 1, 0).*dtV;
            worngidx=[];
            for j=1:length(peaktime)
                idx=round(peaktime(j)/dtV-0.01/dtV:peaktime(j)/dtV);
                temp=min(waveformV(idx,i));
                if temp>waveformV(idx(end),i)-10;
                    worngidx=[worngidx;j];
                end
                
                if waveformV(idx(end),i)<-20&peaktime(j)>0.5+pretime;
                    worngidx=[worngidx;j];
                end
                
            end
            peaktime(worngidx)=[];
            
            if ~isempty(peaktime);peaktime(peaktime>0.51+pretime|peaktime<pretime+0.0005)=[];end
            %             idx=find(diff(peaktime)>0.001);
            %             idx=[1;idx(:)+1];
            %             peaktime=peaktime(idx);
            
            spkFreq=[spkFreq;length(peaktime)/0.5];
            spkTime{i}=peaktime;
            
            figure(1),clf
            plot(tspanV,waveformV(:,i)),hold on
%             plot(tspanI,waveformI(:,i)),hold on
            plot(tspanV(round(peaktime/dtV)),waveformV(round(peaktime/dtV),i),'ro')
            drawnow
            %pause()
        end
        
        temp=waveformV-repmat(Vrest,size(waveformV,1),1);
        VAHP=trapz(temp(tspanV>0.51&tspanV<0.8,:),1);
% %         
%         Vposi=mean(waveformV(tspanV>0.4&tspanV<Dura,:))-mean(waveformV(tspanV<0,:));
%         Vposi(spkFreq~=0)=nan;
        
        % StimAmp
        StimAmp=mean(waveformI(tspanI>0.1&tspanI<0.4,:))-mean(waveformI(tspanI<0,:));

        % PC PV protocol
        StimAmp(StimAmp<200)=round(StimAmp(StimAmp<200)/10)*10;
        StimAmp(StimAmp>=200)=round(StimAmp(StimAmp>=200)/50)*50;
 
        % SST protocol
%         StimAmp(StimAmp<200)=round(StimAmp(StimAmp<200)/10)*10;
%         StimAmp(StimAmp>=200)=round(StimAmp(StimAmp>=200)/100)*100;
        
        % VIP protocol
%         StimAmp(StimAmp<100)=round(StimAmp(StimAmp<100)/10)*10;
%         StimAmp(StimAmp>=100)=round(StimAmp(StimAmp>=100)/50)*50;
%        
        save(['DATA_FIcurve_',filename],'tspanV','tspanI','waveformV','waveformI','TimeStim','Dura',...
            'Vrest','VAHP','spkFreq','spkTime','StimAmp')
        
          
        figure(1),clf
        plot(StimAmp,VAHP,'ro')
        drawnow
    end
    
end
