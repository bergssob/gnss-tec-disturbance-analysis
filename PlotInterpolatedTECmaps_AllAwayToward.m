%TEC map plotting routines.
% 1. TEC map using all LOSs
% 2. TEC maps, 3 panels using: All LOSs, LOSs Away, LOSs Toward

%Required input: NEXRAD data and TEC data in S-RAID format fulltimedata.mat
%--------------------------------------------------------------------------


%Date = [2025 6 9]; %year, month, day


%Loading NEXRAD data
NexRad_str = ['D:\MATLAB\NEXRAD_data/nexrad_' num2str(Date(1)) '_' num2str(Date(2), '%02d') '_' num2str(Date(3), '%02d') '.mat'];
load(NexRad_str)

nexradlat = lat(1:4:end);
nexradlon = lon(1:4:end);
prectime = timefull;

%Loading vTEC data FOR LOCAL
vTEC_str = ['D:\MATLAB\fulltimedata_for_LoS_analysis\fulltimedata_' num2str(Date(2), '%02d') '_' num2str(Date(3), '%02d') '_' num2str(Date(1)) '-GRE.mat'];
load(vTEC_str)

%Continental US REGION
minlat=25;
maxlat=50;
minlon=-125;
maxlon=-75;
llat=200; %200 is 0.125°
llon=400; %400 is 0.125°
elevcut = 20; % Elevation cut-off angle

% % %Europe
% minlat=35;
% maxlat=65;
% minlon=-15;
% maxlon=30;
% llat=200; %200 is 0.15°
% llon=300; %300 is 0.15°
% elevcut = 20; % Elevation cut-off angle


latgrid=linspace(minlat,maxlat,llat);
longrid=linspace(minlon,maxlon,llon);
dlat=latgrid(2)-latgrid(1);
dlon=longrid(2)-longrid(1);
[LON,LAT]=meshgrid(longrid,latgrid);


%defining assistance vector to help plot the NEXRAD data, which helps to match the
%correct 10 minute NEXRAD frame with each 2 minute TEC frame
v = 1:length(timefull);
%u = repelem(v,5);
u = repelem(v,20);%To match NEXRAD with TEC every 30 secs
%u = repelem(v,4);%To match NEXRAD data which has 2 min resolution with TEC every 30 secs


FS = 18; %Fontsize


%filterindx = 7; %10 min
filterindx = 8; %5-40 min
% filterindx = 9; %60 min
% filterindx = 10; %120 min

if filterindx == 7
    filtertext = '10 min high-pass';
    cbar = [-0.05 0.05];
elseif filterindx == 8
    filtertext = '5-40 min band-pass';
    cbar = [-0.1 0.1];
elseif filterindx == 9
    filtertext = '60 min high-pass';
    cbar = [-0.15 0.15];
elseif filterindx == 10
    filtertext = '120 min high-pass';
    cbar = [-1 1];
else
    filtertext = '';
    cbar = [];
end

%Start by making all vTEC values above 1 a nan (because those values are the unrealistic fluctuations at the end values of time intervals).
for tt = 1:length(fulltimedata)-1 %for every 30 secs
    index = abs(fulltimedata{1,tt}(:,filterindx)) > 1;    
    fulltimedata{1,tt}(index,filterindx) = nan;
end



%% 1. TEC map using all LOSs

clear movieoutput
movieoutput(length(UTT)) = struct('cdata', [], 'colormap', []); %To use parfor
parpool('local',10)
parfor t = 1*120:(24*120);

t
tic

currentdata = fulltimedata{1,t}; %t to get every 30 secs.

azimuths = currentdata(:,3);
%Cleaning up nan values
indxx = find(isnan(currentdata(:,filterindx)));
currentdata(indxx,:)=[];
azimuths(indxx)=[];
%Because azimuth runs from -180 to 180 in the data, we want to make it from 0-360.
xxx = find(azimuths<0);
azimuths(xxx)=360+azimuths(xxx);


vTECgrid2=[];
          

 for ilat=1:llat
     for ilon=1:llon
           lattmpg=LAT(ilat,ilon);
           lontmpg=LON(ilat,ilon);                            
                
           inds=find(currentdata(:,4)>lattmpg-dlat/2 & currentdata(:,4)<=lattmpg+dlat/2 ...
                    & currentdata(:,5)>lontmpg-dlon/2 & currentdata(:,5)<=lontmpg+dlon/2 ...
                    & currentdata(:,2)>elevcut); 
                
            if(isempty(inds))
                vTECgrid2(ilat,ilon)=NaN;
            else
                if(all(currentdata(inds,filterindx)>=0))
                    vTECgrid2(ilat,ilon)=max(currentdata(inds,filterindx));
                elseif(all(currentdata(inds,filterindx)<0))
                    vTECgrid2(ilat,ilon)=min(currentdata(inds,filterindx));
                else
                    vTECgrid2(ilat,ilon)=mean(currentdata(inds,filterindx),'omitnan');
                end
            end



                
     end
 end

     
mapaverage = find(abs(vTECgrid2)>1);
if(~isempty(mapaverage))
    vTECgrid2(mapaverage) = NaN;
end

f = vTECgrid2;
f(1,:)=0;
f(end,:)=0;
f(:,1)=0;
f(:,end)=0;
vTECgrid2 = inpaint_nans(f,4);
vTECgrid2 = smooth2a(vTECgrid2,1);

FS = 18;
 hg = figure();
 set(hg,'Position',[10,10,1000,800]);
 set(gcf,'color','white','DefaultFigureVisible','off');
 setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 0, 0]);
        set(gcf,'color','white','visible','off');

% %----------------Plotting NEXRAD overlaid on vTEC----------------------------
ax6 = subplot_tight(1,1,1,[0.15 0.15]); 
%ax6 = subplot_tight(3,2,[1 4],[0.1 0.1]);  %for 3x2
worldmap([minlat maxlat],[minlon maxlon]);
setm(gca,'FontSize',16)
set(gcf, 'CurrentAxes', ax6);
    h=pcolorm(latgrid,longrid,vTECgrid2);
    shading interp
    set(gca,'FontSize',FS);
    set(h,'alphadata',~isnan(vTECgrid2));
    hold on;

    caxis(cbar)
    c=colorbar;
    %c.FontSize = 9;
    colormap(ax6, flipud(lbmap(512,'Blue')))
    title([filtertext ' filtered dvTEC (TECu)']);

    hold off
geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',1);
geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',0.5);
set(gca,'LineWidth',1)
set(gca,'FontSize',FS)

pause(0.5) %For some reasons maps were not aligned unless pausing between them

ax8 = axes('Position', get(ax6, 'Position'), 'Color', 'none');
worldmap([minlat maxlat],[minlon maxlon]);
setm(gca,'FontSize',16)
set(gcf, 'CurrentAxes', ax8);
    h2 = pcolorm(nexradlat,nexradlon,squeeze(precipfull(u(t),:,:))');
    %shading interp
    cmap8 = jet(200);
    cmap8 = (cmap8(101:200,:));
    cmap8(1,:) = [1,1,1];
    %colormap(h8.Parent,cmap8)
    colormap(ax8,cmap8)
    caxis([0 7])
alphaData = 0.6 * double(squeeze(precipfull(u(t),:,:))' ~= 0); % Creating a transparency map, 0 for zero values and 0.6 for non-zero values.
set(h2, 'FaceAlpha', 'flat', 'AlphaData', alphaData, 'AlphaDataMapping', 'none');
setm(ax8, 'Frame', 'off');

%  %Annotation - timestamp:
    dim = [0.36 0.24 0.1 0.6];
    a = annotation('textbox',dim,'String',[datestr(UTT(t)),' UT'],'FitBoxToText','on','BackgroundColor','white');
    a.FontSize = 16;

% %For Europe
%     dim = [0.34 0.34 0.1 0.6];
%     a = annotation('textbox',dim,'String',[datestr(UTT(t)),' UT'],'FitBoxToText','on','BackgroundColor','white');
%     a.FontSize = 16;



movieoutput(t) = getframe(gcf);

close all

toc
end
empty_frames = arrayfun(@(f) isempty(f.cdata), movieoutput);
movieoutput(empty_frames) = [];


save_video_string = ['TEC_Overlaid_NEXRAD_'  num2str(Date(1)) '_' num2str(Date(2), '%02d') '_' num2str(Date(3), '%02d')];

%v2 = VideoWriter(save_video_string2,'Motion JPEG AVI'); %for VEGA
v = VideoWriter(save_video_string,'MPEG-4'); %for local
v.FrameRate = 20; %Frames per second
open(v)
writeVideo(v,movieoutput(:));
close(v);






%% 2. TEC maps, 3 panels using: All LOSs, LOSs Away, LOSs Toward. Interpolated maps

t_ind = 0; %index for the assistent vector for a moving epicenter
%correct 10 minute NEXRAD frame with each 30 sec TEC frame
v = 1:length(t_start:t_end);
shortened_v = v(1:20:end);
new_v = 1:1:length(shortened_v);
v_assistant = repelem(new_v,20);

elevroof = 90;
clear movieoutput
movieoutput(length(UTT)) = struct('cdata', [], 'colormap', []); %To use parfor
parpool('local',10)
%parfor t = 21:6*120+1
parfor t = (2*120+1):(4*120+1)
    
t
tic

currentdata = fulltimedata{1,t}; %t to get every 30 secs.

azimuths = currentdata(:,3);
%Cleaning up nan values
indxx = find(isnan(currentdata(:,filterindx)));
currentdata(indxx,:)=[];
azimuths(indxx)=[];
%Because azimuth runs from -180 to 180 in the data, we want to make it from 0-360.
xxx = find(azimuths<0);
azimuths(xxx)=360+azimuths(xxx);


%epicenter_at_timestep = [26 -100]; %RF01 21st May 2025
%epicenter_at_timestep = [32 -101]; %RF10 7th June 2025
%epicenter_at_timestep = [30 -102]; %RF09 6th June 2025
%epicenter_at_timestep = [33 -100]; %RF11 9th June 2025
%epicenter_at_timestep = [30 -83]; %Hurricane Helene
%epicenter_at_timestep = [33 -98]; %11th March 2026
epicenter_at_timestep = [33.4 -97.3]; %26th May 2024 stationary epicenter

% t_ind = t_ind + 1;
% epicenter_at_timestep = epicenter(v_assistant(t_ind),:); %For moving epicenter


vTECgrid_All = [];
vTECgrid_Away = [];
vTECgrid_Toward = [];
azi_range = 45; %the width of the azimuth range of LOSs towards the center

   for ilat=1:llat
       for ilon=1:llon



                lattmpg=LAT(ilat,ilon);
                lontmpg=LON(ilat,ilon);

                inds=find(currentdata(:,4)>lattmpg-dlat/2 & currentdata(:,4)<=lattmpg+dlat/2 ...
                    & currentdata(:,5)>lontmpg-dlon/2 & currentdata(:,5)<=lontmpg+dlon/2 ...
                    & (currentdata(:,2) > elevcut & currentdata(:,2) < elevroof) );

                    if(isempty(inds))
                        vTECgrid_All(ilat,ilon)=NaN;
                    else
                            if(all(currentdata(inds,filterindx)>=0))
                        vTECgrid_All(ilat,ilon)=max(currentdata(inds,filterindx));
                            elseif(all(currentdata(inds,filterindx)<0))
                        vTECgrid_All(ilat,ilon)=min(currentdata(inds,filterindx));
                            else
                        vTECgrid_All(ilat,ilon)=mean(currentdata(inds,filterindx),'omitnan');
                            end
                    end

%Away from
    azimuth_from_center = azimuth(epicenter_at_timestep(1),epicenter_at_timestep(2),lattmpg,lontmpg); %Now the azimuth is calculated from the center to the grid cell
                        if azimuth_from_center > 360-azi_range
                                inds=find(currentdata(:,4)>lattmpg-dlat/2 & currentdata(:,4)<=lattmpg+dlat/2 ...
                                & currentdata(:,5)>lontmpg-dlon/2 & currentdata(:,5)<=lontmpg+dlon/2 ...
                                & (currentdata(:,2) > elevcut & currentdata(:,2) < elevroof) &  (azimuths > azimuth_from_center-azi_range | azimuths < azi_range - (360-azimuth_from_center) ) );
                        elseif azimuth_from_center < azi_range
                                inds=find(currentdata(:,4)>lattmpg-dlat/2 & currentdata(:,4)<=lattmpg+dlat/2 ...
                                & currentdata(:,5)>lontmpg-dlon/2 & currentdata(:,5)<=lontmpg+dlon/2 ...
                                & (currentdata(:,2) > elevcut & currentdata(:,2) < elevroof) &  (azimuths > 360 - (azi_range - azimuth_from_center) | azimuths < azimuth_from_center+azi_range ) );
                        else
                                inds=find(currentdata(:,4)>lattmpg-dlat/2 & currentdata(:,4)<=lattmpg+dlat/2 ...
                                & currentdata(:,5)>lontmpg-dlon/2 & currentdata(:,5)<=lontmpg+dlon/2 ...
                                & (currentdata(:,2) > elevcut & currentdata(:,2) < elevroof) & (azimuths > azimuth_from_center-azi_range & azimuths < azimuth_from_center+azi_range) );
                        end


                    if(isempty(inds))
                        vTECgrid_Away(ilat,ilon)=NaN;
                    else
                            if(all(currentdata(inds,filterindx)>=0))
                        vTECgrid_Away(ilat,ilon)=max(currentdata(inds,filterindx));
                            elseif(all(currentdata(inds,filterindx)<0))
                        vTECgrid_Away(ilat,ilon)=min(currentdata(inds,filterindx));
                            else
                        vTECgrid_Away(ilat,ilon)=mean(currentdata(inds,filterindx),'omitnan');
                            end
                    end

%Toward
  azimuth_to_center = azimuth(lattmpg,lontmpg,epicenter_at_timestep(1),epicenter_at_timestep(2));
                        if azimuth_to_center > 360-azi_range
                                inds=find(currentdata(:,4)>lattmpg-dlat/2 & currentdata(:,4)<=lattmpg+dlat/2 ...
                                & currentdata(:,5)>lontmpg-dlon/2 & currentdata(:,5)<=lontmpg+dlon/2 ...
                                & (currentdata(:,2) > elevcut & currentdata(:,2) < elevroof) &  (azimuths > azimuth_to_center-azi_range | azimuths < azi_range - (360-azimuth_to_center) ) );
                        elseif azimuth_to_center < azi_range
                                inds=find(currentdata(:,4)>lattmpg-dlat/2 & currentdata(:,4)<=lattmpg+dlat/2 ...
                                & currentdata(:,5)>lontmpg-dlon/2 & currentdata(:,5)<=lontmpg+dlon/2 ...
                                & (currentdata(:,2) > elevcut & currentdata(:,2) < elevroof) &  (azimuths > 360 - (azi_range - azimuth_to_center) | azimuths < azimuth_to_center+azi_range ) );
                        else
                                inds=find(currentdata(:,4)>lattmpg-dlat/2 & currentdata(:,4)<=lattmpg+dlat/2 ...
                                & currentdata(:,5)>lontmpg-dlon/2 & currentdata(:,5)<=lontmpg+dlon/2 ...
                                & (currentdata(:,2) > elevcut & currentdata(:,2) < elevroof) & (azimuths > azimuth_to_center-azi_range & azimuths < azimuth_to_center+azi_range) );
                        end


                    if(isempty(inds))
                        vTECgrid_Toward(ilat,ilon)=NaN;
                    else
                            if(all(currentdata(inds,filterindx)>=0))
                        vTECgrid_Toward(ilat,ilon)=max(currentdata(inds,filterindx));
                            elseif(all(currentdata(inds,filterindx)<0))
                        vTECgrid_Toward(ilat,ilon)=min(currentdata(inds,filterindx));
                            else
                        vTECgrid_Toward(ilat,ilon)=mean(currentdata(inds,filterindx),'omitnan');
                            end
                    end


       end
   end

     
mapaverage = find(abs(vTECgrid_All)>1);
if(~isempty(mapaverage))
    vTECgrid_All(mapaverage) = NaN;
end

mapaverage = find(abs(vTECgrid_Away)>1);
if(~isempty(mapaverage))
    vTECgrid_Away(mapaverage) = NaN;
end

mapaverage = find(abs(vTECgrid_Toward)>1);
if(~isempty(mapaverage))
    vTECgrid_Toward(mapaverage) = NaN;
end


f = vTECgrid_All;
f(1,:)=0;
f(end,:)=0;
f(:,1)=0;
f(:,end)=0;
vTECgrid_All = inpaint_nans(f,4);
vTECgrid_All = smooth2a(vTECgrid_All,1);

f = vTECgrid_Away;
f(1,:)=0;
f(end,:)=0;
f(:,1)=0;
f(:,end)=0;
vTECgrid_Away = inpaint_nans(f,4);
vTECgrid_Away = smooth2a(vTECgrid_Away,1);

f = vTECgrid_Toward;
f(1,:)=0;
f(end,:)=0;
f(:,1)=0;
f(:,end)=0;
vTECgrid_Toward = inpaint_nans(f,4);
vTECgrid_Toward = smooth2a(vTECgrid_Toward,1);

subplot_lims = [0.04 0.04];
FS = 18;
 hg = figure();
 %set(hg,'Position',[10,10,1000,800]);
 set(hg,'Position',[10,10,2000,600]);
 set(gcf,'color','white','DefaultFigureVisible','off');
 setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 0, 0]);
        set(gcf,'color','white','visible','off');

% %----------------Plotting NEXRAD overlaid on vTEC----------------------------
ax6 = subplot_tight(1,3,1,subplot_lims); 
%ax6 = subplot_tight(3,2,[1 4],[0.1 0.1]);  %for 3x2
worldmap([minlat maxlat],[minlon maxlon]);
setm(gca,'FontSize',16)
set(gcf, 'CurrentAxes', ax6);
    h=pcolorm(latgrid,longrid,vTECgrid_All);
    shading interp
    set(gca,'FontSize',FS);
    set(h,'alphadata',~isnan(vTECgrid_All));
    hold on;
    caxis(cbar)
    c=colorbar;
    %ylabel(c, '[TECu]')
    %c.FontSize = 9;
    colormap(ax6, flipud(lbmap(512,'Blue')))
    title({[filtertext ' filtered dvTEC'],'All LOSs'});
hold off
geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',1);
geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',0.5);
set(gca,'LineWidth',1)
set(gca,'FontSize',FS)

pause(0.2) %For some reasons maps were not aligned unless pausing between them

ax8 = axes('Position', get(ax6, 'Position'), 'Color', 'none');
worldmap([minlat maxlat],[minlon maxlon]);
setm(gca,'FontSize',16)
set(gcf, 'CurrentAxes', ax8);
    h2 = pcolorm(nexradlat,nexradlon,squeeze(precipfull(u(t),:,:))');
    %shading interp
    cmap8 = jet(200);
    cmap8 = (cmap8(101:200,:));
    cmap8(1,:) = [1,1,1];
    %colormap(h8.Parent,cmap8)
    colormap(ax8,cmap8)
    caxis([0 7])
alphaData = 0.6 * double(squeeze(precipfull(u(t),:,:))' ~= 0); % Creating a transparency map, 0 for zero values and 0.6 for non-zero values.
set(h2, 'FaceAlpha', 'flat', 'AlphaData', alphaData, 'AlphaDataMapping', 'none');
setm(ax8, 'Frame', 'off');


ax2 = subplot_tight(1,3,2,subplot_lims); 
%ax6 = subplot_tight(3,2,[1 4],[0.1 0.1]);  %for 3x2
worldmap([minlat maxlat],[minlon maxlon]);
setm(gca,'FontSize',16)
set(gcf, 'CurrentAxes', ax2);
    h=pcolorm(latgrid,longrid,vTECgrid_Away);
    shading interp
    set(gca,'FontSize',FS);
    set(h,'alphadata',~isnan(vTECgrid_Away));
    hold on;
    caxis(cbar)
    c=colorbar;
    ylabel(c, 'dvTEC [TECu]')
    %c.FontSize = 9;
    colormap(ax2, flipud(lbmap(512,'Blue')))
    title(['LOSs away from source (Texas)'],'color','b');
hold off
geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',1);
geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',0.5);
set(gca,'LineWidth',1)
set(gca,'FontSize',FS)

ax3 = subplot_tight(1,3,3,subplot_lims); 
%ax6 = subplot_tight(3,2,[1 4],[0.1 0.1]);  %for 3x2
worldmap([minlat maxlat],[minlon maxlon]);
setm(gca,'FontSize',16)
set(gcf, 'CurrentAxes', ax3);
    h=pcolorm(latgrid,longrid,vTECgrid_Toward);
    shading interp
    set(gca,'FontSize',FS);
    set(h,'alphadata',~isnan(vTECgrid_Toward));
    hold on;
    caxis(cbar)
    c=colorbar;
    ylabel(c, '[TECu]')
    %c.FontSize = 9;
    colormap(ax3, flipud(lbmap(512,'Blue')))
    title(['LOSs toward source (Texas)'],'color','r');
hold off
geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',1);
geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',0.5);
set(gca,'LineWidth',1)
set(gca,'FontSize',FS)



    dim = [0.42 0.30 0.1 0.6];
    a = annotation('textbox',dim,'String',[datestr(UTT(t)),' UT'],'FitBoxToText','on','BackgroundColor','white');
    a.FontSize = 16;


movieoutput(t) = getframe(gcf);

close all

toc
end
empty_frames = arrayfun(@(f) isempty(f.cdata), movieoutput);
movieoutput(empty_frames) = [];


save_video_string = ['TEC_Overlaid_NEXRAD_LOSsAllAwayToward_'  num2str(Date(1)) '_' num2str(Date(2), '%02d') '_' num2str(Date(3), '%02d')];
%v2 = VideoWriter(save_video_string2,'Motion JPEG AVI'); %for VEGA
v = VideoWriter(save_video_string,'MPEG-4'); %for local
v.FrameRate = 20; %Frames per second
open(v)
writeVideo(v,movieoutput(:));
close(v);

