%FIND NUMBER OF OBJECTS IN FILLED BINARY IMAGE
numberOfBoundaries=inf;
while numberOfBoundaries>1

    boundaries=bwboundaries(BWfilled); %cell array where each cell contains the row & column coordinates for an object in the image.
    numberOfBoundaries=length(boundaries); %total number of objects found.
    
    if numberOfBoundaries>1 %if there are multiple objects in the image.
        largestBoundary=bwboundaries(largestBlob); %find coordinates of largest object.
        largeBoundary=largestBoundary{1}; %coordinates of largest object.
        xLarge=largeBoundary(:,2); %column coordinates of largest object.
        yLarge=largeBoundary(:,1); %row coordinates of largest object.
        cent=regionprops(largestBlob,'centroid'); %get centroid of largest object.
        refs=cent.Centroid; %store centroid coordinates.   
        hold on; plot(refs(:,1),refs(:,2),'b*') %plot centroid on image.

        BWobjects=bwareafilt(BWfilled,numberOfBoundaries-1,'smallest'); %remove largest object from image.
        boundaries=bwboundaries(BWobjects); %remove the coordinates of the largest object from the list.
    
        overallMinDistance=inf;
        for b=1:numberOfBoundaries %loop through each small object.

            %Find the object closest to the largest object (i.e., closest to the centroid of the largest object).
            x=boundaries{b}(:,1);
            y=boundaries{b}(:,2);
            distances=sqrt((x - refs(1)).^2 + (y - refs(2)).^2); %Pythagorean theorem to find the least distance between the current object and the largest object.
            [minDistance,indexOfMinDistance]=min(distances);
            if minDistance<overallMinDistance %this object is the closest one to the largest object.
                index1=b;
                index2=indexOfMinDistance;
                
                subplot(3,3,7)
                imshow(largestBlob)
                title('All objects connected');
                h=findobj(gcf,'type','axes');
                f=get(h(1),'children'); %get current axes of subplot with largest object.

                hLine=drawline('Position',[[x(index1), xLarge(index2)]; [y(index1), yLarge(index2)]]); %define the line that connects this object to the largest object.
                singleLineBinaryImage=hLine.createMask(); %create a binary image ("mask") from the object.
                largestBlob(singleLineBinaryImage)=255; %burn line into image by setting the values to 255 wherever the mask is true.
            end
        end
    end
    
    subplot(3,3,8);
    imshow(largestBlob,[]); %display the image with all objects connected.
    title('BW continuous');
    
    numberOfBoundaries=numberOfBoundaries-1;
end