% _________________________________________________________________________
% This script provide the answers for the Problem 7.12 
% A pair of cell complexes K, K' covering a pair of different 
% triangulated video frames
%--------------------------------------------------------------------------
% Prepared by Juwairiah Zia
% Student Number: 7895845
% June 9th, 2020
% ECE 8300 Computer Vision
% University of Manitoba
%==========================================================================

clc; close all; clear all;  % Erase all existing variables.
detri = [];csz = [];vf = [];po = [];triindt = []; triind = []; vf1 = []; 
brcnt= []; brcntx= []; brcnty= [];brcnt3= []; brcntx1= []; brcnty1= []; 
prox =[]; bo=[]; b1=[]; b2=[]; p1= [];p0=[];
% dm3=[];rowno=[];rowall=[];rowallf=[];
videoReader = vision.VideoFileReader('test1.mp4');  % Read matlab video for input
writerObj = VideoWriter('p4.7_1.mp4','MPEG-4');  % Create video writer object
writerObj.FrameRate = 25;   % Set the fps
open(writerObj);            % Open the video writer
fn = 250;
prox=zeros(fn,7);   
for count = 1:fn    % The number of frames to process
    disp(count)                         % Display the current frame number
    frame = step(videoReader);          % Read the next video frame
    
    frameGray = rgb2gray(frame);        % Convert the frame to grayscale
    if count ==1
        bkFrame = frameGray;            % Get the first frame as background frame
    end
    
    frameDif = abs(bkFrame - frameGray);    % The difference between the two frames
    frameBW = imbinarize( frameDif , 0.15); % Covert the frame to binary
    frameBW2 = bwareaopen(frameBW,50);      % Remove blobs smaller than 50 (Turn dark foreground to white)
    se = strel('disk',10);                  % Use a disk of radius 10 to dialte and erode object shape
    frameBW3 = imdilate(frameBW2,se);       % Dialate the object shape
    frameBW3 = imerode(frameBW3,se);        % Erode the object shape
    
    s = regionprops(frameBW2,'centroid');   % Get the stats of the blobs
    blbCentroids = cat(1,s.Centroid);       % Exatract the centroids
    
    [B] = bwboundaries(frameBW3);           % Save the object boundaries
    corners = corner(frameBW3);             % Detect corners on edges
    figure('visible','off');                % To run the script faster
    imshow(frame)                           % Display the current frame
    seedPoints = [blbCentroids;corners];
    hold on 

% %   fC=flip(frameCentroids);
%     xp=(frameCentroids(:,1));
%     yp=(frameCentroids(:,2));
    

    %% Create Delaunay triangulation
    if size(seedPoints,1)>3
    detri=delaunayTriangulation(seedPoints);
    %Find the MNC on the triangulated shape
    itr=size(detri.Points,1);%determine number of iteration by counting vertices
    idtr=size(detri.ConnectivityList,1);
    %finding and concatenating all the triangles as per vertices
    for db = 1:itr      %go to ith row of detri.Points and extract dth vertice
        da = vertexAttachments(detri,db); %finds the triangles attached to the dth vertice.
        opt1(db,:) = da;%concatenating da
        %d = ID, da returns the IDs of the triangles
        %e.g. if ID= [3,4] then get IDs from 3 and 4 row in connectivity list. 
        %daa=cell2mat(da);
    end
    
    ab=db;%a way to remove previous values in cell
        for aa=1:ab
            abc=opt1(aa,:);
            vvv=cell2mat(abc);
            csz1(aa,:)=length(vvv);    
        end
           csz=csz1(1:aa,:);%a way to get current length in the frame
           %Find maximum triangles connected together FOR MAXIMUM NEUCLEUS CLUSTER
           mx1 = max(csz);%highest # of row = max of the triangles connected to each other through a vertice
           [row,col] = find(csz==mx1); %extract row number(s) in mx1; saved in row
           triindt = row;%saving row number(s), gives number of MNC(s)
           %for first MNC
           triindt = row;%saving row number(s), gives number of MNC(s)
           nmncs = size(triindt,1);
           for nmncr = 1:nmncs
           %for first MNC
           triind = triindt(nmncr,:);%first MNC/NrvH 
           %%%%%%%
           p1(nmncr,:)=detri.Points(triind,:);%Nucleus of MNC/NrvH
           p=p1(nmncr,:);
           %%%%%%%
           nmnc=size(triind,1);%number of MNC (may have multiple)%%%%%%%%%%%
           %IDs of triangles in MNC
%           if nmnc == 1 
              Vbr = vertexAttachments(detri,triind);
              mn4 = Vbr{:};
              %e.g. mn4=[1,2,3,4...] means cluster contain triangle 1, triangle 2 etc
              %on the connectivity list
              %mn4 contain triangle IDs that form MNC using dab
              %fprintf('tri IDs in MNC= %d \n', mn4);   
              %IDS of every coordinate of triangle in MNC 
              vf=detri(mn4,:);%connectivity list of MNC each row contain IDs of each triangle in mn4
              %to extract coordinates of MNC 
              ntri=size(vf,1);
              for mn=1:ntri
                  trix=vf(mn,:);
                  tricrd=detri.Points(trix,:);%cord of mnth triangle
                  xmnctri = (tricrd(:,1)); %x point list of MNC of all triangles in MNC
                  ymnctri = (tricrd(:,2)); %y point list of MNC of all triangles in MNC
  
                  %barycenter
                  brcnt1(mn,:) = mean(tricrd); %barrycenter = mean of every column containing x and y coordinates of triangles
              end
                  %barycenters of MNC (NrvE)
                  %a way to get current barycenters in the frame
                  brcnt=brcnt1(1:mn,:);%point list of cycE
                  brcntx=brcnt(:,1);
                  brcnty=brcnt(:,2);
                           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %NRVE'
                  %triangles along the border of MNC and finding their barycenters
                  %triangles with points sans neucleus
                  mnc2 = size(vf,1);
                  mnc3 = zeros(mnc2,2);
                  for mncc=1:mnc2
                      mnc1=vf(mncc,:);
                      mnc2=mnc1;
                      mnc2(ismember(mnc2,triind)) = [];
                      mnc3(mncc,:)=mnc2;
                  end 
                      vf4 = mnc3(1:mncc,:);%current MNC points ID without neucleus
                      %Use edge attachments
                      %finding barycenters of triangles opposite to the edges of MNC
                      %finding triangles opposite to the edges of MNC
                      mnnc1=size(vf4,1);
                      for mmnnc = 1:mnnc1
                          ids = vf4(mmnnc,1);ide = vf4(mmnnc,2);%edges IDs from vf4 i.e. MNC without neucleus
                          ID = edgeAttachments(detri,ids,ide);%gives the ID(row#) from connectivity list of detri
                          idmat = cell2mat(ID);
                          idmats = size(idmat,2);%triangle edge should make two triangles
                          if idmats == 1
                             mnnc6(mmnnc,:)= detri.ConnectivityList(idmat,:); 
                             elseif idmats > 1
                          mnnc4 = detri.ConnectivityList(idmat,:);%gives the ID connected to the edge
                          %one triangle will be from MNC, compare it with vf 
                          vf5 = ismember(mnnc4,vf,'rows');
                          %eliminate the one in vf5
                          mnnc5 = mnnc4;
                          [row2,col2] = find(vf5==1);%row numbers
                          %removing MNC triangles
                          mnnc5(row2,:)=[];
                          mnnc6(mmnnc,:)= mnnc5;
                          end
                      end
%connectivity list (like vf) of triangles along the edges of MNC (taking out current)
                          nrv2=mnnc6(1:mmnnc,:);
                          %finding barycenters
                          ntri2=size(nrv2,1);
                          for mnnc7=1:ntri2
                              trix2=nrv2(mnnc7,:);
                              tricrd2=detri.Points(trix2,:);%cord of mnth triangle
                              xnrvtri = (tricrd2(:,1)); %x point list of MNC of all triangles in MNC
                              ynrvtri = (tricrd2(:,2)); %y point list of MNC of all triangles in MNC
                              %barycenter
                              brnrv2(mnnc7,:) = mean(tricrd2); %barrycenter = mean of every column containing x and y coordinates of triangles
                          end
                          %barycenters of triangles opp. to edges of MNC
                          %a way to get current barycenters in the frame
                          brnrv=brnrv2(1:mnnc7,:);%point list of cycE
                          brnrvx=brnrv(:,1);
                          brnrvy=brnrv(:,2);
                          areabr1 = polyarea(brnrvx,brnrvy); 
                          areabr (nmncr,:) =  areabr1;
                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PLOT
%%%%
%triangulation
% hold on
num_tri=size(detri.ConnectivityList,1);%determine number of triangles = number of rows on connectivitylist
for i=1:num_tri      %go to ith row of connectivitylist and take out indices
    coord_tri=detri.Points(detri.ConnectivityList(i,:),:);
    %separating x and y coordinates
    x_tri = (coord_tri(:,1)); %contain x coordinate of all triangles
    y_tri = (coord_tri(:,2)); %contain y coordinate of all triangles
    patch(x_tri,y_tri,'yellow','FaceAlpha',0.005,'Marker','o', ...
    'MarkerFaceColor','w','MarkerSize',2,'MarkerEdgeColor','green', ...
    'EdgeColor','blue','LineWidth',0.5);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%false color bary centers of NrvE' (barycenter along the
%triangles of MNC
     scatter(brnrv(:,1),brnrv(:,2),'s','yellow');%NrvE' points
%      hold on 
     %plot Nrv;
     patch(brnrvx,brnrvy,'g','FaceAlpha',.15, 'EdgeColor','magenta','LineWidth',0.8);%NrvE'   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      hold on 
%false color bary centers of NrvE (barycenters of MNC
     scatter(brcnt(:,1),brcnt(:,2),'s','r');%NrvE points
%      hold on 
     %plot NrvH;
     patch(brcntx,brcnty,'y','FaceAlpha',.4, 'EdgeColor','red','LineWidth',0.8);%NrvE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cusp filament between NrvE and NrvE'
brcntn = size(brcnt,1);
for brcntc = 1:brcntn
    brcntc1 = brcnt(brcntc,:);
    brcntc2 = brnrv(brcntc,:);
    cuspfl=[brcntc1;brcntc2];
    plot(cuspfl(:,1),cuspfl(:,2),'--w');
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hold on
%nucleus of NrvE
plot(p(:,1),p(:,2),'ro');
%%plotting MNC
%hold on 
%triplot(detri(mn4,:),xp,yp,'Color','m')%plot MNC(s)


% %plotting MNC without neucleus
% hold on 
       
       for mnnc = 1:mnnc1
           mnnc2 = vf4(mnnc,:);
           mnnc3 = detri.Points(mnnc2,:);
%            plot(mnnc3(:,1),mnnc3(:,2),'k--+');
       end
%        
%Bo filaments
bfilbrc(nmncr,:) = size(brcnt,1); 
           end
%     hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Feature selection
%B(skcyclicNrvE), holeCount)
%betty number = the number of nerve ribbon filaments (B0(rbE)) 
%+ nerve ribbon cycles (B1(rbE)) + nerve ribbon holes (B2(rbE)). 
%B0
b0 = sum(bfilbrc);
prox(count,1) = b0;
%B1
b1 = nmncs*2;
prox(count,2) = b1;
%B2
b2 = nmncr;
prox(count,3) = b2;
%betti number of nrvH
b012 = b0+b1+b2;
prox(count,4) = b012;
%hole count NrvE
holee = nmncr * mnnc1;
prox(count,5) = holee;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mass of shape
shar = sum(areabr);
prox(count,6) = shar;
%Energy of shape
%velocity = f-f'/t-t'
tc=(count*35)/800;%current time
prox(count,7)=tc;
%mass of shape E
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
hold off

writeFrame = getframe(gcf);              % Capture the current displayed frame
close();                                 % Close the figure
writeVideo(writerObj,writeFrame.cdata);  % Write the captured frame to the video
%     end
end
close(writerObj);   % Close the video writer object
%%%%%%%Comparison step
% fprintf('\phi(skcyclicNrvE)'
%select threshole
%%%%%NrvE
k1 = 'Select frame for skcyclicNrvE=';
e1 = input(k1);%6
%%%%%%NrvE'
k2 = 'Select another frame for skcyclicNrvE=';
e2 = input(k2);%7
%compare if the betti number and hole count of e1 and e2 are less than th
%betti number for e1(skcyclicNrvE)
bet1 = prox(e1,4);
%betti number for e2(skcyclicNrvE')
bet2 = prox(e2,4);
betti = abs(bet1-bet2);
%hole count for e1(skcyclicNrvE)
eh1 = prox(e1,5); 
%hole count for e1(skcyclicNrvE')
eh2 = prox(e2,5);
ehole = abs(eh1-eh2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Energy of shape
%velocity = f-f'/t-t'
sht1=prox(e1,7);%time at K
sht2=prox(e2,7);%time at K'
sht=abs(sht1-sht2);%difference in time
shv1 = abs((e1-e2)/sht);%velocity
shv = shv1^2;%v^2
%kinetic energy for E
shm1 = prox(e1,6);
shke1 = 0.5*shm1*shv;
%kinetic energy for E'
shm2 = prox(e2,6);
shke2 = 0.5*shm2*shv;
shkee = abs(shke1-shke2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('1.(skcyclicNrvE)=(%d, %d)\n',bet1,eh1);
% fprintf('2.(skcyclicNrvE)=(%d, %d)\n',bet2,eh2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('1.(skcyclicNrvE)=(%d, %d)\n',shke1,eh1);
fprintf('2.(skcyclicNrvE)=(%d, %d)\n',shke2,eh2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for threshold
thr1 = 'Set the threshold=';
thr = input(thr1);%5
if shkee<thr || ehole<thr
   des = 'Yes';
   fprintf('%s, the norm of difference between the features is less than threshold %d. \n skcyclicNrvE on frame %d and skcyclicNrvE on frame %d are close to each other approximately descriptively.\n', des,thr,e1,e2);
else
    des = 'No';
    fprintf('%s, the norm of difference between the features is not less than threshold %d. \n skcyclicNrvE on frame %d and skcyclicNrvE on frame %d are not close to each other approximately descriptively.\n', des,thr,e1,e2); 
end