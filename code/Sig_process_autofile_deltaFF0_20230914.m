clc;
clear;
close all;
tic
%%
    fs=10;%Sampling frequency when collecting images(Hz)
    r_um=50;%The radius of the circle used for framing,(um)
    pixsize=0.425;%Pixel size，（um/pix）
    thresh=0.3;%The maximum white value in colormap corresponds to ΔF/F0
    thresh_select=0.12;%Fluorescence region selection threshold
    bp1=1.5;%Breakpoint 1 for calculating photobleaching
    bp2=5;%Breakpoint 2 for calculating photobleaching
%%%%%%%%%%%%%%%%The above are adjustable parameters%%%%%%%%%%%%%%%%%%
folderpath=uigetdir;
namelist = dir([folderpath,'\*.tif']);
r_pix=r_um/pixsize;
dt=1/fs;
dtt=1/fs/30;
%%
for n=1:length(namelist)
    %Image Loading
    Info=imfinfo([namelist(n).folder,'\',namelist(n).name]);
    format=Info.Format;
    Slice=size(Info,1);                                          %%Get the number of frames in the z direction of the image
    Width=Info.Width;
    Height=Info.Height;
    t=0:dt:(Slice-1)/fs;
    disp('Information');
    toc
    Image = double(tiffreadVolume([namelist(n).folder,'\',namelist(n).name]));
    disp('Read Image ready');
    toc
    
    %%
    % Maximum Fluorescence Change Display - High Precision Version
    t_img_bp=bp1/dt:bp2/dt-dt;
    bas=zeros(Height,Width,length(t_img_bp));
    bas_full=zeros(Height,Width,Slice);
    bas=Image(:,:,bp1/dt+1:bp2/dt);
    for i=1:Height
        for j=1:Width
            dif=(bas(i,j,end)-bas(i,j,1))./((bp2-bp1)/dt-1);
            bas_start=bas(i,j,1)-bp1/dt*dif;
            bas_end=bas(i,j,end)+(Slice-bp2/dt-1)*dif;
            temp=zeros(1,1,Slice);
            temp(1,1,:)=linspace(bas_start,bas_end,Slice);  %Extend the length of the bas from the size of the distance between two bp points to the size of the entire image frame size 
            bas_full(i,j,:)=temp;
        end
    end
    clear temp;
    delta_Image=(Image-bas_full)./bas_full;  %delta_Image is the stack of changes in fluorescence
    maximg=max(delta_Image(:,:,5:14/dt),[],3); %The maximum amount of fluorescence change within 5~14s (i.e. the maximum degree of stimulation)
    maximg(maximg<thresh_select)=0;
    mask=im2bw(maximg,0);
    pixnum=length(find(mask>0));
    
    % Maximum fluorescence change display-Rapid Edition
    % aver_50=mean(Image(:,:,10:50),3);%Take the average of the first 50 sheets as the background for subtracting
    % subbackimg=Image-aver_50;%The stack after subtracting the background
    % maximg=max(subbackimg,[],3);
    % maximg=maximg/mean(aver_50,'all');
    % maximg(maximg<thresh_select)=0;
    % mask=im2bw(maximg,0);
    % pixnum=length(find(mask>0));
    
    maximg=maximg./thresh;
    imwrite(maximg,[namelist(n).folder,'\max_',namelist(n).name]);
    figure(666);
    set(figure(666),'position',[100,200,Height,Width]);
    imshow(maximg);
    mymap=[0         0         0
        0.1137         0    0.1765
        0.2275         0    0.3569
        0.3412         0    0.5333
        0.2275         0    0.6902
        0.1137         0    0.8431
             0         0    1.0000
        0.0118    0.2353    0.6863
        0.0196    0.4667    0.3725
        0.0314    0.7020    0.0588
        0.5176    0.8510    0.0314
        1.0000    1.0000         0
        1.0000    0.6667         0
        1.0000    0.3333         0
        1.0000         0         0
        0.9804    0.1059    0.3137
        0.9608    0.2078    0.6275
        0.9412    0.3137    0.9412
        0.9608    0.5412    0.9608
        0.9804    0.7725    0.9804
        1.0000    1.0000    1.0000];
    colormap(mymap);
    labeldis=linspace(0,1,11);
    label=linspace(0,thresh,11);
    bar=colorbar('Ticks',labeldis,'TickLabels',label);
    bar.Label.String='ΔF / F0';
    disp('Image background extraction');
    toc
    
    %%
    %Specific data processing, including signal extraction, curve fitting and plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%The signal is extracted according to the selected area of the box%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Signal=zeros(1,Slice);
    Signal=squeeze(sum(Image,[1 2]));%Signal extraction
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Fitting curve%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fitob=fit(t',Signal,'smoothingspline','SmoothingParam',0.99);%The fitting curve increases the amount of data and the accuracy of calculation results. fitob is the curve model after fitting
    max_t0=find(Signal==max(Signal));%The peak position of the original signal
    fitob2=fit(t(max_t0:end)',Signal(max_t0:end),'smoothingspline','SmoothingParam',0.5);
    tt=0:dtt:Slice*dt-dtt;
    tt_2=bp1:dtt:bp2-dtt;
    lent=length(tt);
    yy=feval(fitob,tt');
    max_t=find(yy==max(yy));%Peak position of the signal after fitting
    yy(max_t:end)=feval(fitob2,tt(max_t:end));%Signals with lower peaks are fitted with smoother curves
    %The baseline is fitted with the signal between the two breakpoints bp1 and bp2
    bp_Sig=yy(bp1/dtt+1:bp2/dtt);
    Cor_Sig=detrend(bp_Sig);
    bias=yy(bp1/dtt+1:bp2/dtt)-Cor_Sig; x
    fitob3=fit(tt_2',bias,'smoothingspline','SmoothingParam',1);
    bias=feval(fitob3,tt');The fitted attenuation curve bas
    %Restore the corrected signal：
    Cor_Sig=100.*(yy-bias)./bias;
    maxsig=max(Cor_Sig(5/dtt:20/dtt));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1);
    set(figure(1),'position',[380,270,length(t)*6,450]);%figure Window position and size adjustment
    subplot(1,2,2);  
    sigmin=min(Cor_Sig,[],'all');
    sigmax=max(Cor_Sig,[],'all');
    ylim([sigmin sigmax]);
    xlabel('Time/s');
    ylabel('ΔF/F0 %');
    plot(tt,Cor_Sig,'Color','r','LineWidth',1.5);hold on;
    text(0.5*tt(end),0.9*max(Cor_Sig),['max ΔF/F0 = ',num2str(maxsig),'%']);
    legend('ΔF/F0');  hold off;
    subplot(1,2,1); 
    plot(t,Signal,'Color',"b");hold on;
    xlabel('Time/s');
    ylabel('Pixvalue');
    legend('raw curve');  hold off;
    res=['max ΔF/F0 = ',num2str(maxsig),'%'];
    toc
    disp(['percentage of pix exceeding the threshold is: ',num2str(100*(pixnum/Height/Width)),'%']);
    disp(res);
    result(n).name=namelist(n).name;
    result(n).Percentage=100*(pixnum/Height/Width);
    result(n).maxdff=maxsig;
end
save([folderpath,'\result.mat'],"result");
