function plotMorphologyParamPM_incTimes(data)
%% PREPARE DATA

%load('dataAnalysis.mat','data');

%SEPARATE DATA ACCORDING TO CONDITION
%Initialise variables.
LPS0h=[];
LPS3h=[];
LPS6h=[];
Vehicle6h=[];

for ii=1:length(data)
    if ~isempty(strfind(data{ii,1},'6h_0'))
        Vehicle6h=[Vehicle6h;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'6h_250'))
        LPS6h=[LPS6h;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'3h'))
        LPS3h=[LPS3h;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'0h'))
        LPS0h=[LPS0h;data{ii,3}];
    end
end

%EXTRACT PARAMETERS FOR EACH CONDITION
%Initialise variables.
LPS0h_Area=zeros(length(LPS0h),1);
LPS3h_Area=zeros(length(LPS3h),1);
LPS6h_Area=zeros(length(LPS6h),1);
Vehicle6h_Area=zeros(length(Vehicle6h),1);

LPS0h_MajorAxisLength=zeros(length(LPS0h),1);
LPS3h_MajorAxisLength=zeros(length(LPS3h),1);
LPS6h_MajorAxisLength=zeros(length(LPS6h),1);
Vehicle6h_MajorAxisLength=zeros(length(Vehicle6h),1);

LPS0h_MinorAxisLength=zeros(length(LPS0h),1);
LPS3h_MinorAxisLength=zeros(length(LPS3h),1);
LPS6h_MinorAxisLength=zeros(length(LPS6h),1);
Vehicle6h_MinorAxisLength=zeros(length(Vehicle6h),1);

LPS0h_Eccentricity=zeros(length(LPS0h),1);
LPS3h_Eccentricity=zeros(length(LPS3h),1);
LPS6h_Eccentricity=zeros(length(LPS6h),1);
Vehicle6h_Eccentricity=zeros(length(Vehicle6h),1);

LPS0h_Circularity=zeros(length(LPS0h),1);
LPS3h_Circularity=zeros(length(LPS3h),1);
LPS6h_Circularity=zeros(length(LPS6h),1);
Vehicle6h_Circularity=zeros(length(Vehicle6h),1);

LPS0h_Extent=zeros(length(LPS0h),1);
LPS3h_Extent=zeros(length(LPS3h),1);
LPS6h_Extent=zeros(length(LPS6h),1);
Vehicle6h_Extent=zeros(length(Vehicle6h),1);

LPS0h_Perimeter=zeros(length(LPS0h),1);
LPS3h_Perimeter=zeros(length(LPS3h),1);
LPS6h_Perimeter=zeros(length(LPS6h),1);
Vehicle6h_Perimeter=zeros(length(Vehicle6h),1);

LPS0h_Elongation=zeros(length(LPS0h),1);
LPS3h_Elongation=zeros(length(LPS3h),1);
LPS6h_Elongation=zeros(length(LPS6h),1);
Vehicle6h_Elongation=zeros(length(Vehicle6h),1);

%CONDITION 1 PARAMETERS
for ii=1:length(LPS0h)
    LPS0h_Area(ii,1)=LPS0h(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS0h_MajorAxisLength(ii,1)=LPS0h(ii).MajorAxisLength*0.18;
    LPS0h_MinorAxisLength(ii,1)=LPS0h(ii).MinorAxisLength*0.18;
    LPS0h_Eccentricity(ii,1)=LPS0h(ii).Eccentricity;
    LPS0h_Circularity(ii,1)=LPS0h(ii).Circularity;
    LPS0h_Extent(ii,1)=LPS0h(ii).Extent;
    LPS0h_Perimeter(ii,1)=LPS0h(ii).Perimeter*0.18;
    LPS0h_Elongation(ii,1)=(LPS0h(ii).MajorAxisLength*0.18)/(LPS0h(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%CONDITION 2 PARAMETERS
for ii=1:length(LPS3h)
    LPS3h_Area(ii,1)=LPS3h(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS3h_MajorAxisLength(ii,1)=LPS3h(ii).MajorAxisLength*0.18;
    LPS3h_MinorAxisLength(ii,1)=LPS3h(ii).MinorAxisLength*0.18;
    LPS3h_Eccentricity(ii,1)=LPS3h(ii).Eccentricity;
    LPS3h_Circularity(ii,1)=LPS3h(ii).Circularity;
    LPS3h_Extent(ii,1)=LPS3h(ii).Extent;
    LPS3h_Perimeter(ii,1)=LPS3h(ii).Perimeter*0.18;
    LPS3h_Elongation(ii,1)=(LPS3h(ii).MajorAxisLength*0.18)/(LPS3h(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%CONDITION 3 PARAMETERS
for ii=1:length(LPS6h)
    LPS6h_Area(ii,1)=LPS6h(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS6h_MajorAxisLength(ii,1)=LPS6h(ii).MajorAxisLength*0.18;
    LPS6h_MinorAxisLength(ii,1)=LPS6h(ii).MinorAxisLength*0.18;
    LPS6h_Eccentricity(ii,1)=LPS6h(ii).Eccentricity;
    LPS6h_Circularity(ii,1)=LPS6h(ii).Circularity;
    LPS6h_Extent(ii,1)=LPS6h(ii).Extent;
    LPS6h_Perimeter(ii,1)=LPS6h(ii).Perimeter*0.18;
    LPS6h_Elongation(ii,1)=(LPS6h(ii).MajorAxisLength*0.18)/(LPS6h(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%CONDITION 4 PARAMETERS
for ii=1:length(Vehicle6h)
    Vehicle6h_Area(ii,1)=Vehicle6h(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    Vehicle6h_MajorAxisLength(ii,1)=Vehicle6h(ii).MajorAxisLength*0.18;
    Vehicle6h_MinorAxisLength(ii,1)=Vehicle6h(ii).MinorAxisLength*0.18;
    Vehicle6h_Eccentricity(ii,1)=Vehicle6h(ii).Eccentricity;
    Vehicle6h_Circularity(ii,1)=Vehicle6h(ii).Circularity;
    Vehicle6h_Extent(ii,1)=Vehicle6h(ii).Extent;
    Vehicle6h_Perimeter(ii,1)=Vehicle6h(ii).Perimeter*0.18;
    Vehicle6h_Elongation(ii,1)=(Vehicle6h(ii).MajorAxisLength*0.18)/(Vehicle6h(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%% REMOVE OUTLIERS

%CONDITION 1
for ii=1:length(LPS0h)
    LPS0h_Area=rmoutliers(LPS0h_Area);
    LPS0h_MajorAxisLength=rmoutliers(LPS0h_MajorAxisLength);
    LPS0h_MinorAxisLength=rmoutliers(LPS0h_MinorAxisLength);
    LPS0h_Eccentricity=rmoutliers(LPS0h_Eccentricity);
    LPS0h_Circularity=rmoutliers(LPS0h_Circularity);
    LPS0h_Extent=rmoutliers(LPS0h_Extent);
    LPS0h_Perimeter=rmoutliers(LPS0h_Perimeter);
    LPS0h_Elongation=rmoutliers(LPS0h_Elongation);
end

%CONDITION 2
for ii=1:length(LPS3h)
    LPS3h_Area=rmoutliers(LPS3h_Area);
    LPS3h_MajorAxisLength=rmoutliers(LPS3h_MajorAxisLength);
    LPS3h_MinorAxisLength=rmoutliers(LPS3h_MinorAxisLength);
    LPS3h_Eccentricity=rmoutliers(LPS3h_Eccentricity);
    LPS3h_Circularity=rmoutliers(LPS3h_Circularity);
    LPS3h_Extent=rmoutliers(LPS3h_Extent);
    LPS3h_Perimeter=rmoutliers(LPS3h_Perimeter);
    LPS3h_Elongation=rmoutliers(LPS3h_Elongation);
end

%CONDITION 3
for ii=1:length(LPS6h)
    LPS6h_Area=rmoutliers(LPS6h_Area);
    LPS6h_MajorAxisLength=rmoutliers(LPS6h_MajorAxisLength);
    LPS6h_MinorAxisLength=rmoutliers(LPS6h_MinorAxisLength);
    LPS6h_Eccentricity=rmoutliers(LPS6h_Eccentricity);
    LPS6h_Circularity=rmoutliers(LPS6h_Circularity);
    LPS6h_Extent=rmoutliers(LPS6h_Extent);
    LPS6h_Perimeter=rmoutliers(LPS6h_Perimeter);
    LPS6h_Elongation=rmoutliers(LPS6h_Elongation);
end

%CONDITION 4
for ii=1:length(Vehicle6h)
    Vehicle6h_Area=rmoutliers(Vehicle6h_Area);
    Vehicle6h_MajorAxisLength=rmoutliers(Vehicle6h_MajorAxisLength);
    Vehicle6h_MinorAxisLength=rmoutliers(Vehicle6h_MinorAxisLength);
    Vehicle6h_Eccentricity=rmoutliers(Vehicle6h_Eccentricity);
    Vehicle6h_Circularity=rmoutliers(Vehicle6h_Circularity);
    Vehicle6h_Extent=rmoutliers(Vehicle6h_Extent);
    Vehicle6h_Perimeter=rmoutliers(Vehicle6h_Perimeter);
    Vehicle6h_Elongation=rmoutliers(Vehicle6h_Elongation);
end

%% PERFORM STATISTICS FOR EACH PARAMETER

%Prepare data for Kruskall Wallis tests.
g1=ones(length(Vehicle6h_Area),1); %condition 1.
g2=2*ones(length(LPS0h_Area),1); %condition 2.
g3=3*ones(length(LPS3h_Area),1); %condition 3.
g4=4*ones(length(LPS6h_Area),1); %condition 4.
group_area=[g1;g2;g3;g4]; %merge conditions in order.

g1=ones(length(Vehicle6h_MajorAxisLength),1); %condition 1.
g2=2*ones(length(LPS0h_MajorAxisLength),1); %condition 2.
g3=3*ones(length(LPS3h_MajorAxisLength),1); %condition 3.
g4=4*ones(length(LPS6h_MajorAxisLength),1); %condition 4.
group_majorAxisLength=[g1;g2;g3;g4]; %merge conditions in order.

g1=ones(length(Vehicle6h_MinorAxisLength),1); %condition 1.
g2=2*ones(length(LPS0h_MinorAxisLength),1); %condition 2.
g3=3*ones(length(LPS3h_MinorAxisLength),1); %condition 3.
g4=4*ones(length(LPS6h_MinorAxisLength),1); %condition 4.
group_minorAxisLength=[g1;g2;g3;g4]; %merge conditions in order.

g1=ones(length(Vehicle6h_Eccentricity),1); %condition 1.
g2=2*ones(length(LPS0h_Eccentricity),1); %condition 2.
g3=3*ones(length(LPS3h_Eccentricity),1); %condition 3.
g4=4*ones(length(LPS6h_Eccentricity),1); %condition 4.
group_eccentricity=[g1;g2;g3;g4]; %merge conditions in order.

g1=ones(length(Vehicle6h_Circularity),1); %condition 1.
g2=2*ones(length(LPS0h_Circularity),1); %condition 2.
g3=3*ones(length(LPS3h_Circularity),1); %condition 3.
g4=4*ones(length(LPS6h_Circularity),1); %condition 4.
group_circularity=[g1;g2;g3;g4]; %merge conditions in order.

g1=ones(length(Vehicle6h_Extent),1); %condition 1.
g2=2*ones(length(LPS0h_Extent),1); %condition 2.
g3=3*ones(length(LPS3h_Extent),1); %condition 3.
g4=4*ones(length(LPS6h_Extent),1); %condition 4.
group_extent=[g1;g2;g3;g4]; %merge conditions in order.

g1=ones(length(Vehicle6h_Perimeter),1); %condition 1.
g2=2*ones(length(LPS0h_Perimeter),1); %condition 2.
g3=3*ones(length(LPS3h_Perimeter),1); %condition 3.
g4=4*ones(length(LPS6h_Perimeter),1); %condition 4.
group_perimeter=[g1;g2;g3;g4]; %merge conditions in order.

g1=ones(length(Vehicle6h_Elongation),1); %condition 1.
g2=2*ones(length(LPS0h_Elongation),1); %condition 2.
g3=3*ones(length(LPS3h_Elongation),1); %condition 3.
g4=4*ones(length(LPS6h_Elongation),1); %condition 4.
group_elongation=[g1;g2;g3;g4]; %merge conditions in order.

%Merge data in order.
area=[Vehicle6h_Area;LPS0h_Area;LPS3h_Area;LPS6h_Area];
majorAxisLength=[Vehicle6h_MajorAxisLength;LPS0h_MajorAxisLength;LPS3h_MajorAxisLength;LPS6h_MajorAxisLength];
minorAxisLength=[Vehicle6h_MinorAxisLength;LPS0h_MinorAxisLength;LPS3h_MinorAxisLength;LPS6h_MinorAxisLength];
eccentricity=[Vehicle6h_Eccentricity;LPS0h_Eccentricity;LPS3h_Eccentricity;LPS6h_Eccentricity];
circularity=[Vehicle6h_Circularity;LPS0h_Circularity;LPS3h_Circularity;LPS6h_Circularity];
extent=[Vehicle6h_Extent;LPS0h_Extent;LPS3h_Extent;LPS6h_Extent];
perimeter=[Vehicle6h_Perimeter;LPS0h_Perimeter;LPS3h_Perimeter;LPS6h_Perimeter];
elongation=[Vehicle6h_Elongation;LPS0h_Elongation;LPS3h_Elongation;LPS6h_Elongation];

%Perform Kruskal-Wallis test.
[p_area,~,stats_area]=kruskalwallis(area,group_area);
[p_majorAxisLength,~,stats_majorAxisLength]=kruskalwallis(majorAxisLength,group_majorAxisLength);
[p_minorAxisLength,~,stats_minorAxisLength]=kruskalwallis(minorAxisLength,group_minorAxisLength);
[p_eccentricity,~,stats_eccentricity]=kruskalwallis(eccentricity,group_eccentricity);
[p_circularity,~,stats_circularity]=kruskalwallis(circularity,group_circularity);
[p_extent,~,stats_extent]=kruskalwallis(extent,group_extent);
[p_perimeter,~,stats_perimeter]=kruskalwallis(perimeter,group_perimeter);
[p_elongation,~,stats_elongation]=kruskalwallis(elongation,group_elongation);
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

x1=2*ones(1,length(Vehicle6h_Area));
x2=4*ones(1,length(LPS0h_Area));
x3=6*ones(1,length(LPS3h_Area));
x4=8*ones(1,length(LPS6h_Area));

plot(x1,Vehicle6h_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(Vehicle6h_Area),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0h_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0h_Area),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS3h_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS3h_Area),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS6h_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS6h_Area),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 9]); 
xticks(1:9);
ylabel('Cell area (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','6hVehicle','','0hLPS','','3hLPS','','6hLPS',''});

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
            x1=4;
        elseif sig_area(ii,1)==3
            x1=6;
        else
            x1=8;
        end
        if sig_area(ii,2)==2
            x2=4;
        elseif sig_area(ii,2)==3
            x2=6;
        else
            x2=8;
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

x1=2*ones(1,length(Vehicle6h_MajorAxisLength));
x2=4*ones(1,length(LPS0h_MajorAxisLength));
x3=6*ones(1,length(LPS3h_MajorAxisLength));
x4=8*ones(1,length(LPS6h_MajorAxisLength));

plot(x1,Vehicle6h_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(Vehicle6h_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0h_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0h_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS3h_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS3h_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS6h_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS6h_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 9]); 
xticks(1:9);
ylabel('Cell MajorAxisLength (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','6hVehicle','','0hLPS','','3hLPS','','6hLPS',''});

%Indicate significance on plot, if any.
sig=c_majorAxisLength(:,6)<0.05; %extract only significant results.
sig_MajorAxisLength=c_majorAxisLength(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_MajorAxisLength,1)
    if sig_MajorAxisLength(ii,6)<0.05
        if sig_MajorAxisLength(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_MajorAxisLength(ii,1)==2
            x1=4;
        elseif sig_MajorAxisLength(ii,1)==3
            x1=6;
        else
            x1=8;
        end
        if sig_MajorAxisLength(ii,2)==2
            x2=4;
        elseif sig_MajorAxisLength(ii,2)==3
            x2=6;
        else
            x2=8;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_MajorAxisLength(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_MajorAxisLength(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_MajorAxisLength(ii,6)<0.01
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

x1=2*ones(1,length(Vehicle6h_MinorAxisLength));
x2=4*ones(1,length(LPS0h_MinorAxisLength));
x3=6*ones(1,length(LPS3h_MinorAxisLength));
x4=8*ones(1,length(LPS6h_MinorAxisLength));

plot(x1,Vehicle6h_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(Vehicle6h_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0h_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0h_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS3h_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS3h_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS6h_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS6h_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 9]); 
xticks(1:9);
ylabel('Cell MinorAxisLength (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','6hVehicle','','0hLPS','','3hLPS','','6hLPS',''});

%Indicate significance on plot, if any.
sig=c_minorAxisLength(:,6)<0.05; %extract only significant results.
sig_MinorAxisLength=c_minorAxisLength(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_MinorAxisLength,1)
    if sig_MinorAxisLength(ii,6)<0.05
        if sig_MinorAxisLength(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_MinorAxisLength(ii,1)==2
            x1=4;
        elseif sig_MinorAxisLength(ii,1)==3
            x1=6;
        else
            x1=8;
        end
        if sig_MinorAxisLength(ii,2)==2
            x2=4;
        elseif sig_MinorAxisLength(ii,2)==3
            x2=6;
        else
            x2=8;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_MinorAxisLength(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_MinorAxisLength(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_MinorAxisLength(ii,6)<0.01
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

x1=2*ones(1,length(Vehicle6h_Eccentricity));
x2=4*ones(1,length(LPS0h_Eccentricity));
x3=6*ones(1,length(LPS3h_Eccentricity));
x4=8*ones(1,length(LPS6h_Eccentricity));

plot(x1,Vehicle6h_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(Vehicle6h_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0h_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0h_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS3h_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS3h_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS6h_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS6h_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 9]); 
xticks(1:9);
ylabel('Cell Eccentricity (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','6hVehicle','','0hLPS','','3hLPS','','6hLPS',''});

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
            x1=4;
        elseif sig_eccentricity(ii,1)==3
            x1=6;
        else
            x1=8;
        end
        if sig_eccentricity(ii,2)==2
            x2=4;
        elseif sig_eccentricity(ii,2)==3
            x2=6;
        else
            x2=8;
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

x1=2*ones(1,length(Vehicle6h_Circularity));
x2=4*ones(1,length(LPS0h_Circularity));
x3=6*ones(1,length(LPS3h_Circularity));
x4=8*ones(1,length(LPS6h_Circularity));

plot(x1,Vehicle6h_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(Vehicle6h_Circularity),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0h_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0h_Circularity),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS3h_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS3h_Circularity),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS6h_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS6h_Circularity),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 9]); 
xticks(1:9);
ylabel('Cell Circularity (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','6hVehicle','','0hLPS','','3hLPS','','6hLPS',''});

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
            x1=4;
        elseif sig_circularity(ii,1)==3
            x1=6;
        else
            x1=8;
        end
        if sig_circularity(ii,2)==2
            x2=4;
        elseif sig_circularity(ii,2)==3
            x2=6;
        else
            x2=8;
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

x1=2*ones(1,length(Vehicle6h_Extent));
x2=4*ones(1,length(LPS0h_Extent));
x3=6*ones(1,length(LPS3h_Extent));
x4=8*ones(1,length(LPS6h_Extent));

plot(x1,Vehicle6h_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(Vehicle6h_Extent),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0h_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0h_Extent),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS3h_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS3h_Extent),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS6h_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS6h_Extent),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 9]); 
xticks(1:9);
ylabel('Cell Extent (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','6hVehicle','','0hLPS','','3hLPS','','6hLPS',''});

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
            x1=4;
        elseif sig_extent(ii,1)==3
            x1=6;
        else
            x1=8;
        end
        if sig_extent(ii,2)==2
            x2=4;
        elseif sig_extent(ii,2)==3
            x2=6;
        else
            x2=8;
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

x1=2*ones(1,length(Vehicle6h_Perimeter));
x2=4*ones(1,length(LPS0h_Perimeter));
x3=6*ones(1,length(LPS3h_Perimeter));
x4=8*ones(1,length(LPS6h_Perimeter));

plot(x1,Vehicle6h_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(Vehicle6h_Perimeter),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0h_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0h_Perimeter),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS3h_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS3h_Perimeter),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS6h_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS6h_Perimeter),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 9]); 
xticks(1:9);
ylabel('Cell Perimeter (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','6hVehicle','','0hLPS','','3hLPS','','6hLPS',''});

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
            x1=4;
        elseif sig_perimeter(ii,1)==3
            x1=6;
        else
            x1=8;
        end
        if sig_perimeter(ii,2)==2
            x2=4;
        elseif sig_perimeter(ii,2)==3
            x2=6;
        else
            x2=8;
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

x1=2*ones(1,length(Vehicle6h_Elongation));
x2=4*ones(1,length(LPS0h_Elongation));
x3=6*ones(1,length(LPS3h_Elongation));
x4=8*ones(1,length(LPS6h_Elongation));

plot(x1,Vehicle6h_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(Vehicle6h_Elongation),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS0h_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS0h_Elongation),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS3h_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS3h_Elongation),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS6h_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS6h_Elongation),'kx','MarkerSize',18); %plot mean for condition 4.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 9]); 
xticks(1:9);
ylabel('Cell Elongation (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','6hVehicle','','0hLPS','','3hLPS','','6hLPS',''});

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
            x1=4;
        elseif sig_elongation(ii,1)==3
            x1=6;
        else
            x1=8;
        end
        if sig_elongation(ii,2)==2
            x2=4;
        elseif sig_elongation(ii,2)==3
            x2=6;
        else
            x2=8;
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