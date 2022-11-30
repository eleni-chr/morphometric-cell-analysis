function plotMorphologyParamDyn
%% PREPARE DATA
load('dataAnalysis.mat','data');

%SEPARATE DATA ACCORDING TO CONDITION
%Initialise variables.
LPS0_Dyn0=[];
LPS0_Dyn25=[];
LPS10_Dyn0=[];
LPS10_Dyn25=[];

for ii=1:length(data)
    if ~isempty(strfind(data{ii,1},'_10ugLPS+DMSO'))
        LPS10_Dyn0=[LPS10_Dyn0;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'_10ugLPS+Dynarrestin'))
        LPS10_Dyn25=[LPS10_Dyn25;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'_0ugLPS+DMSO'))
        LPS0_Dyn0=[LPS0_Dyn0;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'_0ugLPS+Dyanarrestin'))
        LPS0_Dyn25=[LPS0_Dyn25;data{ii,3}];
    end
end

%EXTRACT PARAMETERS FOR EACH CONDITION
%Initialise variables.
LPS0_Dyn0_Area=zeros(length(LPS0_Dyn0),1);
LPS0_Dyn25_Area=zeros(length(LPS0_Dyn25),1);
LPS10_Dyn0_Area=zeros(length(LPS10_Dyn0),1);
LPS10_Dyn25_Area=zeros(length(LPS10_Dyn25),1);

LPS0_Dyn0_MajorAxisLength=zeros(length(LPS0_Dyn0),1);
LPS0_Dyn25_MajorAxisLength=zeros(length(LPS0_Dyn25),1);
LPS10_Dyn0_MajorAxisLength=zeros(length(LPS10_Dyn0),1);
LPS10_Dyn25_MajorAxisLength=zeros(length(LPS10_Dyn25),1);

LPS0_Dyn0_MinorAxisLength=zeros(length(LPS0_Dyn0),1);
LPS0_Dyn25_MinorAxisLength=zeros(length(LPS0_Dyn25),1);
LPS10_Dyn0_MinorAxisLength=zeros(length(LPS10_Dyn0),1);
LPS10_Dyn25_MinorAxisLength=zeros(length(LPS10_Dyn25),1);

LPS0_Dyn0_Eccentricity=zeros(length(LPS0_Dyn0),1);
LPS0_Dyn25_Eccentricity=zeros(length(LPS0_Dyn25),1);
LPS10_Dyn0_Eccentricity=zeros(length(LPS10_Dyn0),1);
LPS10_Dyn25_Eccentricity=zeros(length(LPS10_Dyn25),1);

LPS0_Dyn0_Circularity=zeros(length(LPS0_Dyn0),1);
LPS0_Dyn25_Circularity=zeros(length(LPS0_Dyn25),1);
LPS10_Dyn0_Circularity=zeros(length(LPS10_Dyn0),1);
LPS10_Dyn25_Circularity=zeros(length(LPS10_Dyn25),1);

LPS0_Dyn0_Extent=zeros(length(LPS0_Dyn0),1);
LPS0_Dyn25_Extent=zeros(length(LPS0_Dyn25),1);
LPS10_Dyn0_Extent=zeros(length(LPS10_Dyn0),1);
LPS10_Dyn25_Extent=zeros(length(LPS10_Dyn25),1);

LPS0_Dyn0_Perimeter=zeros(length(LPS0_Dyn0),1);
LPS0_Dyn25_Perimeter=zeros(length(LPS0_Dyn25),1);
LPS10_Dyn0_Perimeter=zeros(length(LPS10_Dyn0),1);
LPS10_Dyn25_Perimeter=zeros(length(LPS10_Dyn25),1);

LPS0_Dyn0_Elongation=zeros(length(LPS0_Dyn0),1);
LPS0_Dyn25_Elongation=zeros(length(LPS0_Dyn25),1);
LPS10_Dyn0_Elongation=zeros(length(LPS10_Dyn0),1);
LPS10_Dyn25_Elongation=zeros(length(LPS10_Dyn25),1);

%CONDITION 1 PARAMETERS
for ii=1:length(LPS0_Dyn0)
    LPS0_Dyn0_Area(ii,1)=LPS0_Dyn0(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS0_Dyn0_MajorAxisLength(ii,1)=LPS0_Dyn0(ii).MajorAxisLength*0.18;
    LPS0_Dyn0_MinorAxisLength(ii,1)=LPS0_Dyn0(ii).MinorAxisLength*0.18;
    LPS0_Dyn0_Eccentricity(ii,1)=LPS0_Dyn0(ii).Eccentricity;
    LPS0_Dyn0_Circularity(ii,1)=LPS0_Dyn0(ii).Circularity;
    LPS0_Dyn0_Extent(ii,1)=LPS0_Dyn0(ii).Extent;
    LPS0_Dyn0_Perimeter(ii,1)=LPS0_Dyn0(ii).Perimeter*0.18;
    LPS0_Dyn0_Elongation(ii,1)=LPS0_Dyn0(ii).MajorAxisLength*0.18/LPS0_Dyn0(ii).MinorAxisLength*0.18; %calculation additional parameter.
end

%CONDITION 2 PARAMETERS
for ii=1:length(LPS0_Dyn25)
    LPS0_Dyn25_Area(ii,1)=LPS0_Dyn25(ii).Area*0.18;
    LPS0_Dyn25_MajorAxisLength(ii,1)=LPS0_Dyn25(ii).MajorAxisLength*0.18;
    LPS0_Dyn25_MinorAxisLength(ii,1)=LPS0_Dyn25(ii).MinorAxisLength*0.18;
    LPS0_Dyn25_Eccentricity(ii,1)=LPS0_Dyn25(ii).Eccentricity;
    LPS0_Dyn25_Circularity(ii,1)=LPS0_Dyn25(ii).Circularity;
    LPS0_Dyn25_Extent(ii,1)=LPS0_Dyn25(ii).Extent;
    LPS0_Dyn25_Perimeter(ii,1)=LPS0_Dyn25(ii).Perimeter*0.18;
    LPS0_Dyn25_Elongation(ii,1)=LPS0_Dyn25(ii).MajorAxisLength*0.18/LPS0_Dyn25(ii).MinorAxisLength*0.18; %calculation additional parameter.
end

%CONDITION 3 PARAMETERS
for ii=1:length(LPS10_Dyn0)
    LPS10_Dyn0_Area(ii,1)=LPS10_Dyn0(ii).Area*0.18;
    LPS10_Dyn0_MajorAxisLength(ii,1)=LPS10_Dyn0(ii).MajorAxisLength*0.18;
    LPS10_Dyn0_MinorAxisLength(ii,1)=LPS10_Dyn0(ii).MinorAxisLength*0.18;
    LPS10_Dyn0_Eccentricity(ii,1)=LPS10_Dyn0(ii).Eccentricity;
    LPS10_Dyn0_Circularity(ii,1)=LPS10_Dyn0(ii).Circularity;
    LPS10_Dyn0_Extent(ii,1)=LPS10_Dyn0(ii).Extent;
    LPS10_Dyn0_Perimeter(ii,1)=LPS10_Dyn0(ii).Perimeter*0.18;
    LPS10_Dyn0_Elongation(ii,1)=LPS10_Dyn0(ii).MajorAxisLength*0.18/LPS10_Dyn0(ii).MinorAxisLength*0.18; %calculation additional parameter.
end

%CONDITION 4 PARAMETERS
for ii=1:length(LPS10_Dyn25)
    LPS10_Dyn25_Area(ii,1)=LPS10_Dyn25(ii).Area*0.18;
    LPS10_Dyn25_MajorAxisLength(ii,1)=LPS10_Dyn25(ii).MajorAxisLength*0.18;
    LPS10_Dyn25_MinorAxisLength(ii,1)=LPS10_Dyn25(ii).MinorAxisLength*0.18;
    LPS10_Dyn25_Eccentricity(ii,1)=LPS10_Dyn25(ii).Eccentricity;
    LPS10_Dyn25_Circularity(ii,1)=LPS10_Dyn25(ii).Circularity;
    LPS10_Dyn25_Extent(ii,1)=LPS10_Dyn25(ii).Extent;
    LPS10_Dyn25_Perimeter(ii,1)=LPS10_Dyn25(ii).Perimeter*0.18;
    LPS10_Dyn25_Elongation(ii,1)=LPS10_Dyn25(ii).MajorAxisLength*0.18/LPS10_Dyn25(ii).MinorAxisLength*0.18; %calculation additional parameter.
end

%% REMOVE OUTLIERS

%CONDITION 1
for ii=1:length(LPS0_Dyn0)
    LPS0_Dyn0_Area=rmoutliers(LPS0_Dyn0_Area);
    LPS0_Dyn0_MajorAxisLength=rmoutliers(LPS0_Dyn0_MajorAxisLength);
    LPS0_Dyn0_MinorAxisLength=rmoutliers(LPS0_Dyn0_MinorAxisLength);
    LPS0_Dyn0_Eccentricity=rmoutliers(LPS0_Dyn0_Eccentricity);
    LPS0_Dyn0_Circularity=rmoutliers(LPS0_Dyn0_Circularity);
    LPS0_Dyn0_Extent=rmoutliers(LPS0_Dyn0_Extent);
    LPS0_Dyn0_Perimeter=rmoutliers(LPS0_Dyn0_Perimeter);
    LPS0_Dyn0_Elongation=rmoutliers(LPS0_Dyn0_Elongation);
end

%CONDITION 2
for ii=1:length(LPS0_Dyn25)
    LPS0_Dyn25_Area=rmoutliers(LPS0_Dyn25_Area);
    LPS0_Dyn25_MajorAxisLength=rmoutliers(LPS0_Dyn25_MajorAxisLength);
    LPS0_Dyn25_MinorAxisLength=rmoutliers(LPS0_Dyn25_MinorAxisLength);
    LPS0_Dyn25_Eccentricity=rmoutliers(LPS0_Dyn25_Eccentricity);
    LPS0_Dyn25_Circularity=rmoutliers(LPS0_Dyn25_Circularity);
    LPS0_Dyn25_Extent=rmoutliers(LPS0_Dyn25_Extent);
    LPS0_Dyn25_Perimeter=rmoutliers(LPS0_Dyn25_Perimeter);
    LPS0_Dyn25_Elongation=rmoutliers(LPS0_Dyn25_Elongation);
end

%CONDITION 3
for ii=1:length(LPS10_Dyn0)
    LPS10_Dyn0_Area=rmoutliers(LPS10_Dyn0_Area);
    LPS10_Dyn0_MajorAxisLength=rmoutliers(LPS10_Dyn0_MajorAxisLength);
    LPS10_Dyn0_MinorAxisLength=rmoutliers(LPS10_Dyn0_MinorAxisLength);
    LPS10_Dyn0_Eccentricity=rmoutliers(LPS10_Dyn0_Eccentricity);
    LPS10_Dyn0_Circularity=rmoutliers(LPS10_Dyn0_Circularity);
    LPS10_Dyn0_Extent=rmoutliers(LPS10_Dyn0_Extent);
    LPS10_Dyn0_Perimeter=rmoutliers(LPS10_Dyn0_Perimeter);
    LPS10_Dyn0_Elongation=rmoutliers(LPS10_Dyn0_Elongation);
end

%CONDITION 4
for ii=1:length(LPS10_Dyn25)
    LPS10_Dyn25_Area=rmoutliers(LPS10_Dyn25_Area);
    LPS10_Dyn25_MajorAxisLength=rmoutliers(LPS10_Dyn25_MajorAxisLength);
    LPS10_Dyn25_MinorAxisLength=rmoutliers(LPS10_Dyn25_MinorAxisLength);
    LPS10_Dyn25_Eccentricity=rmoutliers(LPS10_Dyn25_Eccentricity);
    LPS10_Dyn25_Circularity=rmoutliers(LPS10_Dyn25_Circularity);
    LPS10_Dyn25_Extent=rmoutliers(LPS10_Dyn25_Extent);
    LPS10_Dyn25_Perimeter=rmoutliers(LPS10_Dyn25_Perimeter);
    LPS10_Dyn25_Elongation=rmoutliers(LPS10_Dyn25_Elongation);
end

%% PERFORM TWO-WAY ANOVAS FOR EACH PARAMETER

%Prepare data for Kruskall Wallis tests.
g1=ones(length(LPS0_Dyn0_Area),1); %condition 1.
g2=2*ones(length(LPS0_Dyn25_Area),1); %condition 2.
g3=3*ones(length(LPS10_Dyn0_Area),1); %condition 3.
g4=4*ones(length(LPS10_Dyn25_Area),1); %condition 4.
group_Area=[g1;g2;g3;g4]; %merge conditions in order.

g1=ones(length(LPS0_Dyn0_MajorAxisLength),1); %condition 1.
g2=2*ones(length(LPS0_Dyn25_MajorAxisLength),1); %condition 2.
g3=3*ones(length(LPS10_Dyn0_MajorAxisLength),1); %condition 3.
g4=4*ones(length(LPS10_Dyn25_MajorAxisLength),1); %condition 4.
group_MajorAxisLength=[g1;g2;g3;g4]; %merge conditions in order.

g1=ones(length(LPS0_Dyn0_MinorAxisLength),1); %condition 1.
g2=2*ones(length(LPS0_Dyn25_MinorAxisLength),1); %condition 2.
g3=3*ones(length(LPS10_Dyn0_MinorAxisLength),1); %condition 3.
g4=4*ones(length(LPS10_Dyn25_MinorAxisLength),1); %condition 4.
group_MinorAxisLength=[g1;g2;g3;g4]; %merge conditions in order.

g1=ones(length(LPS0_Dyn0_Eccentricity),1); %condition 1.
g2=2*ones(length(LPS0_Dyn25_Eccentricity),1); %condition 2.
g3=3*ones(length(LPS10_Dyn0_Eccentricity),1); %condition 3.
g4=4*ones(length(LPS10_Dyn25_Eccentricity),1); %condition 4.
group_Eccentricity=[g1;g2;g3;g4]; %merge conditions in order.

g1=ones(length(LPS0_Dyn0_Circularity),1); %condition 1.
g2=2*ones(length(LPS0_Dyn25_Circularity),1); %condition 2.
g3=3*ones(length(LPS10_Dyn0_Circularity),1); %condition 3.
g4=4*ones(length(LPS10_Dyn25_Circularity),1); %condition 4.
group_Circularity=[g1;g2;g3;g4]; %merge conditions in order.

g1=ones(length(LPS0_Dyn0_Extent),1); %condition 1.
g2=2*ones(length(LPS0_Dyn25_Extent),1); %condition 2.
g3=3*ones(length(LPS10_Dyn0_Extent),1); %condition 3.
g4=4*ones(length(LPS10_Dyn25_Extent),1); %condition 4.
group_Extent=[g1;g2;g3;g4]; %merge conditions in order.

g1=ones(length(LPS0_Dyn0_Perimeter),1); %condition 1.
g2=2*ones(length(LPS0_Dyn25_Perimeter),1); %condition 2.
g3=3*ones(length(LPS10_Dyn0_Perimeter),1); %condition 3.
g4=4*ones(length(LPS10_Dyn25_Perimeter),1); %condition 4.
group_Perimeter=[g1;g2;g3;g4]; %merge conditions in order.

g1=ones(length(LPS0_Dyn0_Elongation),1); %condition 1.
g2=2*ones(length(LPS0_Dyn25_Elongation),1); %condition 2.
g3=3*ones(length(LPS10_Dyn0_Elongation),1); %condition 3.
g4=4*ones(length(LPS10_Dyn25_Elongation),1); %condition 4.
group_Elongation=[g1;g2;g3;g4]; %merge conditions in order.

%Merge data in order.
area=[LPS0_Dyn0_Area;LPS0_Dyn25_Area;LPS10_Dyn0_Area;LPS10_Dyn25_Area];
majorAxisLength=[LPS0_Dyn0_MajorAxisLength;LPS0_Dyn25_MajorAxisLength;LPS10_Dyn0_MajorAxisLength;LPS10_Dyn25_MajorAxisLength];
minorAxisLength=[LPS0_Dyn0_MinorAxisLength;LPS0_Dyn25_MinorAxisLength;LPS10_Dyn0_MinorAxisLength;LPS10_Dyn25_MinorAxisLength];
eccentricity=[LPS0_Dyn0_Eccentricity;LPS0_Dyn25_Eccentricity;LPS10_Dyn0_Eccentricity;LPS10_Dyn25_Eccentricity];
circularity=[LPS0_Dyn0_Circularity;LPS0_Dyn25_Circularity;LPS10_Dyn0_Circularity;LPS10_Dyn25_Circularity];
extent=[LPS0_Dyn0_Extent;LPS0_Dyn25_Extent;LPS10_Dyn0_Extent;LPS10_Dyn25_Extent];
perimeter=[LPS0_Dyn0_Perimeter;LPS0_Dyn25_Perimeter;LPS10_Dyn0_Perimeter;LPS10_Dyn25_Perimeter];
elongation=[LPS0_Dyn0_Elongation;LPS0_Dyn25_Elongation;LPS10_Dyn0_Elongation;LPS10_Dyn25_Elongation];

%Perform Kruskal-Wallis test.
[p_area,~,stats_area]=kruskalwallis(area,group_Area);
[p_majorAxisLength,~,stats_majorAxisLength]=kruskalwallis(majorAxisLength,group_MajorAxisLength);
[p_minorAxisLength,~,stats_minorAxisLength]=kruskalwallis(minorAxisLength,group_MinorAxisLength);
[p_eccentricity,~,stats_eccentricity]=kruskalwallis(eccentricity,group_Eccentricity);
[p_circularity,~,stats_circularity]=kruskalwallis(circularity,group_Circularity);
[p_extent,~,stats_extent]=kruskalwallis(extent,group_Extent);
[p_perimeter,~,stats_perimeter]=kruskalwallis(perimeter,group_Perimeter);
[p_elongation,~,stats_elongation]=kruskalwallis(elongation,group_Elongation);
close all

%Perform multiple comparison tests for all parameters.
c_area=multcompare(stats_area,'CType','hsd','Display','off'); %first two columns are the groups being compared. Last column is the p-value.
c_majorAxisLength=multcompare(stats_majorAxisLength,'CType','hsd','Display','off');
c_minorAxisLength=multcompare(stats_minorAxisLength,'CType','hsd','Display','off');
c_eccentricity=multcompare(stats_eccentricity,'CType','hsd','Display','off');
c_circularity=multcompare(stats_circularity,'CType','hsd','Display','off');
c_extent=multcompare(stats_extent,'CType','hsd','Display','off');
c_perimeter=multcompare(stats_perimeter,'CType','hsd','Display','off');
c_elongation=multcompare(stats_elongation,'CType','hsd','Display','off');

save('KruskalWallisResults.mat','p_area','p_majorAxisLength','p_minorAxisLength','p_eccentricity','p_circularity','p_extent','p_perimeter','p_elongation');
save('multCompResults.mat','c_area','c_majorAxisLength','c_minorAxisLength','c_eccentricity','c_circularity','c_extent','c_perimeter','c_elongation');

%% PLOT AREA FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Dyn0_Area));
x2=5*ones(1,length(LPS0_Dyn25_Area));
x3=8*ones(1,length(LPS10_Dyn0_Area));
x4=11*ones(1,length(LPS10_Dyn25_Area));

plot(x1,LPS0_Dyn0_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Dyn0_Area),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0_Dyn25_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0_Dyn25_Area),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS10_Dyn0_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS10_Dyn0_Area),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS10_Dyn25_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS10_Dyn25_Area),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 12]);
xticks(0:12);
ylabel('Cell area (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','Vehicle','','','Dynarrestin','','','LPS','','','LPS + Dynarrestin',''});

%Indicate significance on plot, if any.
sig=c_area(:,6)<0.05; %extract only significant results.
sig_area=c_area(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_area,1)
    if sig_area(ii,6)<0.05
        if sig_area(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_area(ii,1)==2
            x1=5;
        elseif sig_area(ii,1)==3
            x1=8;
        else
            x1=11;
        end
        if sig_area(ii,2)==2
            x2=5;
        elseif sig_area(ii,2)==3
            x2=8;
        else
            x2=11;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_area(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_area(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_area(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'CellArea'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'CellArea.svg'); %save figure as an SVG file.
saveas(fig,'CellArea.jpg'); %save figure as an JPEG file.
close

%% PLOT MAJOR AXIS LENGTH FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Dyn0_MajorAxisLength));
x2=5*ones(1,length(LPS0_Dyn25_MajorAxisLength));
x3=8*ones(1,length(LPS10_Dyn0_MajorAxisLength));
x4=11*ones(1,length(LPS10_Dyn25_MajorAxisLength));

plot(x1,LPS0_Dyn0_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Dyn0_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0_Dyn25_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0_Dyn25_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS10_Dyn0_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS10_Dyn0_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS10_Dyn25_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS10_Dyn25_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 12]);
xticks(0:12);
ylabel('Major axis length (μm)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','Vehicle','','','Dynarrestin','','','LPS','','','LPS + Dynarrestin',''});

%Indicate significance on plot, if any.
sig=c_majorAxisLength(:,6)<0.05; %extract only significant results.
sig_majorAxisLength=c_majorAxisLength(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_majorAxisLength,1)
    if sig_majorAxisLength(ii,6)<0.05
        if sig_majorAxisLength(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_majorAxisLength(ii,1)==2
            x1=5;
        elseif sig_majorAxisLength(ii,1)==3
            x1=8;
        else
            x1=11;
        end
        if sig_majorAxisLength(ii,2)==2
            x2=5;
        elseif sig_majorAxisLength(ii,2)==3
            x2=8;
        else
            x2=11;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_majorAxisLength(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_majorAxisLength(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_majorAxisLength(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'MajorAxisLength'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'MajorAxisLength.svg'); %save figure as an SVG file.
saveas(fig,'MajorAxisLength.jpg'); %save figure as an JPEG file.
close

%% PLOT MINOR AXIS LENGTH FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Dyn0_MinorAxisLength));
x2=5*ones(1,length(LPS0_Dyn25_MinorAxisLength));
x3=8*ones(1,length(LPS10_Dyn0_MinorAxisLength));
x4=11*ones(1,length(LPS10_Dyn25_MinorAxisLength));

plot(x1,LPS0_Dyn0_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Dyn0_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0_Dyn25_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0_Dyn25_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS10_Dyn0_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS10_Dyn0_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS10_Dyn25_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS10_Dyn25_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 12]);
xticks(0:12);
ylabel('Minor axis length (μm)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','Vehicle','','','Dynarrestin','','','LPS','','','LPS + Dynarrestin',''});

%Indicate significance on plot, if any.
sig=c_minorAxisLength(:,6)<0.05; %extract only significant results.
sig_minorAxisLength=c_minorAxisLength(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_minorAxisLength,1)
    if sig_minorAxisLength(ii,6)<0.05
        if sig_minorAxisLength(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_minorAxisLength(ii,1)==2
            x1=5;
        elseif sig_minorAxisLength(ii,1)==3
            x1=8;
        else
            x1=11;
        end
        if sig_minorAxisLength(ii,2)==2
            x2=5;
        elseif sig_minorAxisLength(ii,2)==3
            x2=8;
        else
            x2=11;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_minorAxisLength(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_minorAxisLength(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_minorAxisLength(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'MinorAxisLength'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'MinorAxisLength.svg'); %save figure as an SVG file.
saveas(fig,'MinorAxisLength.jpg'); %save figure as an JPEG file.
close

%% PLOT ECCENTRICITY FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Dyn0_Eccentricity));
x2=5*ones(1,length(LPS0_Dyn25_Eccentricity));
x3=8*ones(1,length(LPS10_Dyn0_Eccentricity));
x4=11*ones(1,length(LPS10_Dyn25_Eccentricity));

plot(x1,LPS0_Dyn0_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Dyn0_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0_Dyn25_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0_Dyn25_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS10_Dyn0_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS10_Dyn0_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS10_Dyn25_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS10_Dyn25_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 12]);
xticks(0:12);
ylabel('Eccentricity','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','Vehicle','','','Dynarrestin','','','LPS','','','LPS + Dynarrestin',''});

%Indicate significance on plot, if any.
sig=c_eccentricity(:,6)<0.05; %extract only significant results.
sig_eccentricity=c_eccentricity(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_eccentricity,1)
    if sig_eccentricity(ii,6)<0.05
        if sig_eccentricity(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_eccentricity(ii,1)==2
            x1=5;
        elseif sig_eccentricity(ii,1)==3
            x1=8;
        else
            x1=11;
        end
        if sig_eccentricity(ii,2)==2
            x2=5;
        elseif sig_eccentricity(ii,2)==3
            x2=8;
        else
            x2=11;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_eccentricity(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_eccentricity(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_eccentricity(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'Eccentricity'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'Eccentricity.svg'); %save figure as an SVG file.
saveas(fig,'Eccentricity.jpg'); %save figure as an JPEG file.
close

%% PLOT CIRCULARITY FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Dyn0_Circularity));
x2=5*ones(1,length(LPS0_Dyn25_Circularity));
x3=8*ones(1,length(LPS10_Dyn0_Circularity));
x4=11*ones(1,length(LPS10_Dyn25_Circularity));

plot(x1,LPS0_Dyn0_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Dyn0_Circularity),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0_Dyn25_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0_Dyn25_Circularity),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS10_Dyn0_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS10_Dyn0_Circularity),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS10_Dyn25_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS10_Dyn25_Circularity),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 12]);
xticks(0:12);
ylabel('Circularity','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','Vehicle','','','Dynarrestin','','','LPS','','','LPS + Dynarrestin',''});

%Indicate significance on plot, if any.
sig=c_circularity(:,6)<0.05; %extract only significant results.
sig_circularity=c_circularity(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_circularity,1)
    if sig_circularity(ii,6)<0.05
        if sig_circularity(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_circularity(ii,1)==2
            x1=5;
        elseif sig_circularity(ii,1)==3
            x1=8;
        else
            x1=11;
        end
        if sig_circularity(ii,2)==2
            x2=5;
        elseif sig_circularity(ii,2)==3
            x2=8;
        else
            x2=11;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_circularity(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_circularity(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_circularity(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'Circularity'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'Circularity.svg'); %save figure as an SVG file.
saveas(fig,'Circularity.jpg'); %save figure as an JPEG file.
close

%% PLOT EXTENT FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Dyn0_Extent));
x2=5*ones(1,length(LPS0_Dyn25_Extent));
x3=8*ones(1,length(LPS10_Dyn0_Extent));
x4=11*ones(1,length(LPS10_Dyn25_Extent));

plot(x1,LPS0_Dyn0_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Dyn0_Extent),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0_Dyn25_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0_Dyn25_Extent),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS10_Dyn0_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS10_Dyn0_Extent),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS10_Dyn25_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS10_Dyn25_Extent),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 12]);
xticks(0:12);
ylabel('Extent','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','Vehicle','','','Dynarrestin','','','LPS','','','LPS + Dynarrestin',''});

%Indicate significance on plot, if any.
sig=c_extent(:,6)<0.05; %extract only significant results.
sig_extent=c_extent(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_extent,1)
    if sig_extent(ii,6)<0.05
        if sig_extent(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_extent(ii,1)==2
            x1=5;
        elseif sig_extent(ii,1)==3
            x1=8;
        else
            x1=11;
        end
        if sig_extent(ii,2)==2
            x2=5;
        elseif sig_extent(ii,2)==3
            x2=8;
        else
            x2=11;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_extent(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_extent(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_extent(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'Extent'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'Extent.svg'); %save figure as an SVG file.
saveas(fig,'Extent.jpg'); %save figure as an JPEG file.
close

%% PLOT PERIMETER FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Dyn0_Perimeter));
x2=5*ones(1,length(LPS0_Dyn25_Perimeter));
x3=8*ones(1,length(LPS10_Dyn0_Perimeter));
x4=11*ones(1,length(LPS10_Dyn25_Perimeter));

plot(x1,LPS0_Dyn0_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Dyn0_Perimeter),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0_Dyn25_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0_Dyn25_Perimeter),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS10_Dyn0_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS10_Dyn0_Perimeter),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS10_Dyn25_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS10_Dyn25_Perimeter),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 12]);
xticks(0:12);
ylabel('Perimeter (μm)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','Vehicle','','','Dynarrestin','','','LPS','','','LPS + Dynarrestin',''});

%Indicate significance on plot, if any.
sig=c_perimeter(:,6)<0.05; %extract only significant results.
sig_perimeter=c_perimeter(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_perimeter,1)
    if sig_perimeter(ii,6)<0.05
        if sig_perimeter(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_perimeter(ii,1)==2
            x1=5;
        elseif sig_perimeter(ii,1)==3
            x1=8;
        else
            x1=11;
        end
        if sig_perimeter(ii,2)==2
            x2=5;
        elseif sig_perimeter(ii,2)==3
            x2=8;
        else
            x2=11;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_perimeter(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_perimeter(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_perimeter(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'Perimeter'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'Perimeter.svg'); %save figure as an SVG file.
saveas(fig,'Perimeter.jpg'); %save figure as an JPEG file.
close

%% PLOT ELONGATION FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Dyn0_Elongation));
x2=5*ones(1,length(LPS0_Dyn25_Elongation));
x3=8*ones(1,length(LPS10_Dyn0_Elongation));
x4=11*ones(1,length(LPS10_Dyn25_Elongation));

plot(x1,LPS0_Dyn0_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Dyn0_Elongation),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0_Dyn25_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0_Dyn25_Elongation),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS10_Dyn0_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS10_Dyn0_Elongation),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS10_Dyn25_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS10_Dyn25_Elongation),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 12]);
xticks(0:12);
ylabel('Elongation (μm)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','Vehicle','','','Dynarrestin','','','LPS','','','LPS + Dynarrestin',''});

%Indicate significance on plot, if any.
sig=c_elongation(:,6)<0.05; %extract only significant results.
sig_elongation=c_elongation(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_elongation,1)
    if sig_elongation(ii,6)<0.05
        if sig_elongation(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_elongation(ii,1)==2
            x1=5;
        elseif sig_elongation(ii,1)==3
            x1=8;
        else
            x1=11;
        end
        if sig_elongation(ii,2)==2
            x2=5;
        elseif sig_elongation(ii,2)==3
            x2=8;
        else
            x2=11;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_elongation(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_elongation(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_elongation(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'Elongation'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'Elongation.svg'); %save figure as an SVG file.
saveas(fig,'Elongation.jpg'); %save figure as an JPEG file.
close

end