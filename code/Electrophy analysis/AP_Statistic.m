function [spkwaveform,ttspan,...
    HW, AHP, Threshold, Amp,AIS_maxslope,Soma_maxslope...
    HWidx,AHPidx, Thresholdidx,min_slope3idx] = AP_Statistic(data,peaktime,tspan,Goodidx)
% peaktime=spkTime{Wantidx(jj)};
% tspan=tspanV;
% Goodidx=1:min(Num);
%INPUT
% data :   a vector
% tspan:  a vector,have the same elememt with data
% NumSpk: successive number of spike to be analized
%%
dt=tspan(end)-tspan(end-1);
spktime=peaktime(Goodidx);
pretime=0.005;
postime=0.05;

spkwaveform=[];
for i=1:length(spktime)
    idx=round(spktime(i)/dt-pretime/dt:spktime(i)/dt+postime/dt);
    ttspan=linspace(-pretime,postime,numel(idx));
    spkwaveform=[spkwaveform,data(idx)];
end
%%
for i=1:length(spktime)
    data=spkwaveform(:,i);
    slope=[nan;diff(data)/(dt*1000)];
    slope2=[nan;diff(slope)/(dt*1000)];
    slope2=smooth(slope2,5);
    slope3=[nan;diff(slope2)/(dt*1000)];
    slope3=smooth(slope3,5);
    % find slope3的两个trough，分别对应AIS和soma相的max slope
    timewindow=[1:pretime/dt];
    min_slope3idx(i,:)=peakfinder(slope3(timewindow), 5, -4000, -1, 0);
    if length(min_slope3idx(i,:))>1    
        AIS_maxslope(i)=slope3(min_slope3idx(i,1));
        Soma_maxslope(i)=slope3(min_slope3idx(i,2));
    else
        AIS_maxslope(i)=nan;
        Soma_maxslope(i)=slope3(min_slope3idx(i,1));
    end
    
    idxA=find(slope(:)>20);
    idxB=find(data(:)>-60);
    idxC=intersect(idxA(:),idxB(:));
    
     %  if length(spktime)>1
    if length(spktime)>1&i<length(peaktime) % change
        endtime=(peaktime(i+1)-peaktime(i));
    else
        endtime=postime*0.3;
    end  
    
    if ~isempty(idxC)
        Thresholdidx(i)=idxC(1);
        Threshold(i)=data(Thresholdidx(i));
        AHPidx(i)=find(data(ttspan>0&ttspan<endtime)==min(data(ttspan>0&ttspan<endtime)),1,'first')+sum(ttspan<0);
        AHP(i)=data(Thresholdidx(i))-data(AHPidx(i));
        Amp(i)=max(data)-data(Thresholdidx(i));
        HWidx(i,:) = eventedgefinder(data,data(Thresholdidx(i))+Amp(i)/2,1/dt,0.0001,0.0001,1,0);
        HW(i)=(HWidx(i,2)-HWidx(i,1))*dt*1000;
    else
        Thresholdidx(i)=nan;
        Threshold(i)=nan;
        AHPidx(i)=nan;
        AHP(i)=nan;
        Amp(i)=nan;
        HWidx(i,:)=nan;
        HW(i)=nan;
    end    
end

end
