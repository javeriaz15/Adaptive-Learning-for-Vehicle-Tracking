
%==========================================================================
% This script provide object detection, segmentation and tracking using adaptive learning
% The approach involve frame subtraction, morphological operations, canny
% edge detector and using computation geometry for instance segmentation of
% object(s) and plotting centroids for object tracking
% Outline to display segmented object
% Mask to display segmented object
% -------------------------------------------------------------------------
% Prepared by James F. Peters & Juwairiah Zia
%==========================================================================
clear all; close all;
videoReader = vision.VideoFileReader('sample.mp4');  % Read sample matlab video
writerObj = VideoWriter('output.mp4','MPEG-4');        % Create video writer object
writerObj.FrameRate = 24;                                   % Set the fps (Optional)
open(writerObj);                                            % Open the video writer
framen = 350;
feature = zeros(framen,4);
for count = 1:framen                       % The number of frames to process
    feature(count,1)= count;            % Saving frame number in a feature vector
    disp(count)                         % Display the current frame number
    frame = step(videoReader);          % Read the next video frame
    
    frameGray = rgb2gray(frame);        % Convert the frame to grayscale
    if count ==1
        bkFrame = frameGray;            % Get the first frame as background frame
    end
    
        frameDif = abs(bkFrame - frameGray);    % The difference between the two frames
        frameBW1 = imbinarize(frameDif , 0.15); % Covert the frame to binary
        frameBW = bwareaopen(frameBW1,50);      % Remove blobs smaller than 50 (Turn dark foreground to white)
        se = strel('square',15);                % Use a square of radius 10 to dialte and erode object shape
        frameDil = imdilate(frameBW,se);        % Dialate the object shape
        frameEr = imerode(frameDil,se);         % Erode the object shape
        thisFrame = imfill(frameEr,'holes');    % holes filled CC
        
        frameEdge1 = edge(thisFrame,'Canny');   % apply canny to detect edges
        frameEdge = imdilate(frameEdge1,se);    % dialate the edges
    
        figure('visible','off');
        subplot(2,3,1);imshow(frameDif);title('frame subtraction');
        subplot(2,3,2);imshow(thisFrame);title('Morph oprt');
        subplot(2,3,3);imshow(frameEdge);title('Edges:Canny');
        subplot(2,3,4);imshow(frame);title('Corners');
        subplot(2,3,5);imshow(frame);title('Segmented:Mask');
        subplot(2,3,6);imshow(frame);title('Segmented:Triangulation');

        B=[];
        [B] = bwboundaries(thisFrame);         % Save the object boundaries

        corners =[];
        corners = corner(frameEdge, 30);       % Detect corners on edges 
        subplot(2,3,4);imshow(frame);title('Corners');hold on
        plot(corners(:,1),corners(:,2),'mo','MarkerSize',3,'MarkerEdgeColor','w','MarkerFaceColor','m');
        

        %% Object Segmentation: Mask and triangulation
        % Create Delaunay triangulation
        if size(corners,1)>3
            segObj = labeloverlay(frame,thisFrame);
            subplot(2,3,5);imshow(segObj);title('Segmented:Mask');
            
%           % triangulation of segmented object
            s = regionprops(frameEdge,frameGray,'Centroid');   % Get the stats of the blobs
            blbCentroids = cat(1,s.Centroid);                  % Exatract the centroids
            feature(count,2)= size(blbCentroids,1);            % Saving number of centroid in a feature vector
            blbCentroid = blbCentroids(1,:);
            feature(count,[3,4])= blbCentroid;                 % Saving centroid in a feature vector
            seedPoints = [blbCentroid;corners];
            delaunay_tri=delaunay(seedPoints(:,1),seedPoints(:,2)); % Generating delaunay triangles

            % Extract the triangle vertices
            x_tri=[]; % Clear the variable (Not to redraw triangles from earlier frame)
            y_tri=[];
            for i=1:size(delaunay_tri,1)
                x_tri(:,i)=seedPoints(delaunay_tri(i,:),1)';
                y_tri(:,i)=seedPoints(delaunay_tri(i,:),2)';
            end
            
            % Plotting the triangles and the vertices on the RGB frame
            subplot(2,3,6);imshow(frame);title('Segmented:Triangulation');hold on
            faceclr='yellow';           % False color the triangles
            faceopac=0.5;               % Triangle color opacity
            vrtxmrkr='o';               % Triangle vertex marker type
            vrtxmrkrsize=3;             % Triangle vertex marker size
            vrtxmrkrclr='g';            % Triangle vertex marker fill color (Boundary points)
            mrkredgeclr=[0 0 0];        % Tringle vertex edge color
            edgeclr=[0.3 0.3 0.3];      % Triangle edge color
            edgewidth=1;                % Triangle edge width

            patch(x_tri,y_tri,faceclr,'FaceAlpha',faceopac,'Marker',vrtxmrkr,'MarkerFaceColor',vrtxmrkrclr ...
                ,'MarkerSize',vrtxmrkrsize,'MarkerEdgeColor',mrkredgeclr,'EdgeColor',edgeclr,'LineWidth',edgewidth);
            hold on
            
            % Plotting centroid of segmented object
            plot(blbCentroid(:,1),blbCentroid(:,2),'r+');
            hold off
        end
    
    writeFrame = getframe(gcf);             % Capture the current displayed frame
    close();                                % Close the figure
    writeVideo(writerObj,writeFrame.cdata); % Write the captured frame to the video
end

close(writerObj); % Close the video writer object

%%
% xcoord = feature(:,2);
% ycoord = feature(:,3);
% zcoord = feature(:,1);
% figure
% subplot(1,2,1)
% plot3(xcoord,ycoord,zcoord,'-o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
% % axis equal
% xlabel('x of cntr');ylabel('y of of cntr');zlabel('frame');title('Using plot3');
% 
% subplot(1,2,2)
% scatter3(xcoord,ycoord,zcoord,'filled')
% xlabel('x of cntr');ylabel('y of of cntr');zlabel('frame');title('Using scatter3');
% view(40,35)