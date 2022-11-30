function morphology(b)
%% Function written by Eleni Christoforidou in MATLAB R2020b.

%This function performs image processing of individual cells previously
%extracted from confocal images in ImageJ. Then uses the processed images
%to calculte several morphology parameters for each cell.

%Run this function from inside the folder containing the subfolders with
%the images to be analysed. The images must begin with 'cell' in the
%filename, and must be in TIFF format. They must also be a single frame
%(i.e., not a z-stack).

%INPUT ARGUMENTS:
%b     A vector of integers, each representing the maximum background
%      pixel value of each image to be analysed. The values must be in the
%      same order as the order that the images are analysed in. If this
%      information is in an Excel file, use something like 
%      b=xlsread('20201105.xlsx','D2:D122'); to get the data into MATLAB.

%For EACH image analysed, the following are saved in the same directory as
%the image:
%(1) A FIG file showing the image processing performed.
%(2) A MAT file with the same filename as the image, containing the 'props'
%    variable, which contains all the calculations performed for that cell.

%For ALL images analysed, the following are saved in the parent directory:
%(1) A MAT file called 'dataAnalysis', containing the combined data from
%    ALL images analysed.

%% FIND ALL IMAGES TO ANALYSE
folders=dir('*/cell*.tif');
idx=[];
for i=1:length(folders)
    if isfolder(folders(i).name)
        idx=[idx;i]; %Ignore items that are not folders
    end
end
folders=folders(idx);
wd=cd;
back=0;

for ii=1:length(folders)
    cd(folders(ii).folder);
    fn=folders(ii).name;
    back=back+1;
    thresh=b(back);

    %LOAD IMAGE OF CELL AND PLOT IT
    I=imread(fn); %read original image.
    subplot(2,3,1)
    imshow(I)
    title('Original image');

    set(gcf,'WindowStyle','docked');

    %REMOVE BACKGROUND PIXELS
    I(I<=thresh)=0; %set pixels below the background threshold to 0.
    subplot(2,3,2)
    imshow(I)
    title('Background removed');

    %BINARISE THE IMAGE AND PLOT IT
    BW=imbinarize(I,'adaptive','Sensitivity',1);
    subplot(2,3,3)
    imshow(BW)
    title('Binary image');

    %FILL ANY HOLES IN THE CELL AFTER IMAGE IS BINARISED
    BWfilled=imfill(BW,'holes');
    subplot(2,3,4)
    imshow(BWfilled)
    title('Holes filled');

    %EXTRACT LARGEST OBJECT IN IMAGE (I.E., THE MAIN PART OF THE CELL)
    largestBlob=bwareafilt(BWfilled,1);
    subplot(2,3,5)
    imshow(largestBlob)
    title('Largest object');
    
    %EXTRACT ALL OBJECTS APART FROM LARGEST IN IMAGE (I.E., ALL THE PARTS OF
    %THE CELL THAT ARE NOT USED FOR THE CALCULATIONS BELOW)
    boundaries=bwboundaries(BWfilled); %cell array where each cell contains the row & column coordinates for an object in the image.
    numberOfBoundaries=length(boundaries); %total number of objects found.
    smallestBlobs=bwareafilt(BWfilled,numberOfBoundaries-1,'smallest'); %extract all objects apart from the largest one.
    subplot(2,3,6)
    imshow(smallestBlobs)
    title('All other objects');

    %SAVE FIGURE
    filename=fn(1:end-4);
    savefig(filename);
    close

    %CALUCLATE CELL PROPERTIES (using the largest object in the image only)
    props=regionprops(largestBlob,'Area','Circularity','Eccentricity','Extent','MajorAxisLength','MinorAxisLength','Perimeter');
    %Fields in props variable are as follows:
       %Area of the cell.
       %Circularity of the cell.
       %Eccentricity=0 is a circle, while eccentricity=1 is a line. Shows how much of an elliptical shape the cell has.
       %Extent is the ratio of pixels in the object to pixels in the total bounding box (i.e., how far the cell extents).
       %MajorAxisLength is the length (in pixels) of the major axis of the cell.
       %MinorAxisLength is the length (in pixels) of the minor axis of the cell.
       %Perimeter of the cell.
        
    %Calculate additional parameters and append to props variable.
    [N,R]=boxcount(largestBlob);
    FD=log(N)/log(R); %calculate fractal dimension.
    
    CH=bwconvhull(BWfilled); %convex hull.
    CHprops=regionprops(CH,'Area','Perimeter','MajorAxisLength','MinorAxisLength','Circularity');
    
    props.CHA=CHprops.Area; %CHA = convex hull area.
    props.CHP=CHprops.Perimeter; %CHP = convex hull perimeter.
    props.Density=props.Area/CHprops.Area;
    props.Roughness=props.Perimeter/CHprops.Perimeter;
    props.CHspan=props.MajorAxisLength/props.MinorAxisLength; %CHspan = convex hull span ratio.
    props.CHC=CHprops.Circularity; %CHC = convex hull circularity.
    props.FD=FD; %FD = fractal dimension.
    
    %SAVE CELL PROPERTIES
    save(strcat(filename,'.mat'),'props');    
    cd(wd);
end
clearvars -except wd;

%% COMBINE DATA FROM ALL IMAGES TO ONE VARIABLE
len=length(dir('*/cell*.mat'));
data=cell(len,3);
r=1;
fols=dir;

for fol=3:length(fols)
    folname=fols(fol).name;
    cd(fols(fol).name);
    mats=dir('cell*.mat');

    for jj=1:length(mats)
        fname=mats(jj).name;
        load(fname,'props');
        fname=fname(1:end-4);

        %APPEND DATA TO MASTER VARIABLE
        data{r,1}=folname; %column 1 is the sample name.
        data{r,2}=fname; %column 2 is the cell ID.
        data{r,3}=props; %column 3 is the 'props' variable.
        r=r+1;
    end
    cd(wd);
end

%SAVE MASTER VARIABLE
save(strcat('dataAnalysis','.mat'),'data');
clear
end