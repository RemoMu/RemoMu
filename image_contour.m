clc 
clear
close all
%% Instructions
calibrationfactor=0.009406;
ratio=0.189;         %ratio of wall thickness to outer diameter (t/D) 
filenum=[0:2:204];   %according segmentation data
image_spacing=0.4;   %has not to be the actual one, choose for better alphashape
threshold=8;         %alpharadius threshold (run lines 288-289 for check alphashape)
                     %create folders Data and Images, add data and change
                     %paths in the section below
                     
                     %add parseXML and Polyhedron_Plane to folder
                     
                     %code is only capable for 4 calc contours, an error
                     %message should occure when this number is exceeded 
                     %and the code is finished (variable error=1)
                     
                     %A centerline txt file is generated through the lumen
                     %centroid
                     
                     %by Remo Muri (remomuri@msn.com; remo.muri@unibe.ch)
%% Set the different folder PATH

% required to locate Polyhedron_Plan (because of multiple "cd" to data
% folders during execution
addpath('C:\Users\remom\Studium\Master_Thesis\Matlab\PT_12'); 

% path the the data and image folders
path_data_root    = 'C:\Users\remom\Studium\Master_Thesis\Matlab\PT_12\Data';
path_data_folder  = 'C:\Users\remom\Studium\Master_Thesis\Matlab\PT_12\Data';
path_image_folder = 'C:\Users\remom\Studium\Master_Thesis\Matlab\PT_12\Images';

%%
firstcalc=1;
error=0;
for i=1:length(filenum)
cd(path_data_folder)
if (filenum(i)>=10 && filenum(i)<100)   %lumen contour name
    lumcon=sprintf('trans00%d.lumConAut',filenum(i));
elseif (filenum(i)>=0 && filenum(i)<10)  
    lumcon=sprintf('trans000%d.lumConAut',filenum(i));
else
    lumcon=sprintf('trans0%d.lumConAut',filenum(i));
end
lum=importdata(lumcon);         %import lumen contour
nlum=lum(1);
lum(1)=[];                      %delete first
lum(2*nlum+1)=[];               %delete last
lum=reshape(lum,2,nlum)';       %n by 2 matrix
lum(nlum+1,1)=lum(1,1);         %close lumen
lum(nlum+1,2)=lum(1,2);         %close lumen
lum=lum*calibrationfactor;      %convert to [mm] (Calibration factor)

     lumin = polyshape({lum(:,1)'},{lum(:,2)'}, 'simplify', false);
     [x_lum,y_lum] = centroid(lumin);
     
     cd(path_data_folder)
     fid1=fopen('Centerline.txt','a+');     %plots celterline
     fprintf(fid1,'%d, %d\n',x_lum,y_lum);
     fclose(fid1);
     cd(path_data_folder)
     
     lum_x=lum(:,1);
     lum_y=lum(:,2);
     lum_z=linspace(image_spacing*i-image_spacing,image_spacing*i-image_spacing,length(lum));
     if i==1
     lum_x_new=[];
     lum_y_new=[];
     lum_z_new=[];
     end
     lum_x_new=[lum_x_new;(lum_x)];   %move back to original
     lum_y_new=[lum_y_new;(lum_y)];   %move back to original
     lum_z_new=[lum_z_new;(lum_z)'];

if (filenum(i)<100)             %calc contour name
    calccon=sprintf('trans00%d.strutPoints',filenum(i));
else
    calccon=sprintf('trans0%d.strutPoints',filenum(i));
end 
if isfile(calccon)              %import calc contour if file exist
     cd(path_data_folder)
     c=[];                      %delete previous
     calc=parseXML(calccon);
     F=getfield(calc,{1,1},'Children');
     G=getfield(F,{1,2},'Children');
     H=getfield(G,{1,2},'Children');
     
     length_calc=(size(H,2)-1)/2;
     for j=1:length_calc        %load all coordinates
         I=getfield(H,{1,2*j},'Attributes');
         X=getfield(I,{1,1},'Value');
         Y=getfield(I,{1,2},'Value');
         c(j,:)=[str2num(X),str2num(Y)];
     end
     c(length_calc+1,:)=c(1,:); %close contour
     c=c*calibrationfactor;
     polyin = polyshape({lum(:,1)',c(:,1)'},{lum(:,2)',c(:,2)'}, 'simplify', false);
     [x,y] = centroid(polyin);
       
     calc_x=c(:,1);
     calc_y=c(:,2);     
     calc_z=linspace(image_spacing*i-image_spacing,image_spacing*i-image_spacing,length(calc_x));
     if firstcalc==1
     calc_x_new=[];
     calc_y_new=[];
     calc_z_new=[];
     calc_x_new2=[];             %for 2nd calc contour
     calc_y_new2=[];
     calc_z_new2=[];
     calc_x_new3=[];             %for 2nd calc contour
     calc_y_new3=[];
     calc_z_new3=[];
     calc_x_new4=[];             %for 2nd calc contour
     calc_y_new4=[];
     calc_z_new4=[];
     firstcalc=2;
     end
     calc_x_new=[calc_x_new;(calc_x)];   %move back to original
     calc_y_new=[calc_y_new;(calc_y)];   %move back to original
     calc_z_new=[calc_z_new;(calc_z)'];
     
     A_c=polyarea(c(:,1)',c(:,2)');
     A_lum=polyarea(lum(:,1)',lum(:,2)');
     
%2nd calc contour
     if size(G,2)>3
     H2=getfield(G,{1,4},'Children');
          length_calc2=(size(H2,2)-1)/2;
     else
          length_calc2=0;
     end
     if length_calc2>10
         c2=[];
     for j=1:length_calc2          %load all coordinates
         I2=getfield(H2,{1,2*j},'Attributes');
         X2=getfield(I2,{1,1},'Value');
         Y2=getfield(I2,{1,2},'Value');
         c2(j,:)=[str2num(X2),str2num(Y2)];
     end
     c2(length_calc2+1,:)=c2(1,:); %close contour
     c2=c2*calibrationfactor;
     polyin2 = polyshape({lum(:,1)',c(:,1)',c2(:,1)'},{lum(:,2)',c(:,2)',c2(:,2)'}, 'simplify', false);
     [x,y] = centroid(polyin2);
       
     calc_x2=c2(:,1);
     calc_y2=c2(:,2);     
     calc_z2=linspace(image_spacing*i-image_spacing,image_spacing*i-image_spacing,length(calc_x2));

     calc_x_new2=[calc_x_new2;(calc_x2)];   %move back to original
     calc_y_new2=[calc_y_new2;(calc_y2)];   %move back to original
     calc_z_new2=[calc_z_new2;(calc_z2)'];
     
     A_c2=polyarea(c2(:,1)',c2(:,2)');
     else
     c2=[];
     A_c2=0;
     end
     
%3rd calc contour
     if size(G,2)>5
     H3=getfield(G,{1,6},'Children');
          length_calc3=(size(H3,2)-1)/2;
     else
          length_calc3=0;
     end
     if length_calc3>10
         c3=[];
     for j=1:length_calc3         %load all coordinates
         I3=getfield(H3,{1,2*j},'Attributes');
         X3=getfield(I3,{1,1},'Value');
         Y3=getfield(I3,{1,2},'Value');
         c3(j,:)=[str2num(X3),str2num(Y3)];
     end
     c3(length_calc3+1,:)=c3(1,:); %close contour
     c3=c3*calibrationfactor;
     polyin2 = polyshape({lum(:,1)',c(:,1)',c2(:,1)',c3(:,1)'},{lum(:,2)',c(:,2)',c2(:,2)',c3(:,2)'}, 'simplify', false);
     [x,y] = centroid(polyin2);
       
     calc_x3=c3(:,1);
     calc_y3=c3(:,2);     
     calc_z3=linspace(image_spacing*i-image_spacing,image_spacing*i-image_spacing,length(calc_x3));

     calc_x_new3=[calc_x_new3;(calc_x3)];   %move back to original
     calc_y_new3=[calc_y_new3;(calc_y3)];   %move back to original
     calc_z_new3=[calc_z_new3;(calc_z3)'];
     
     A_c3=polyarea(c3(:,1)',c3(:,2)');
     else
     c3=[];
     A_c3=0;
     end
     
%4th calc contour
     if size(G,2)>7
     H4=getfield(G,{1,8},'Children');
          length_calc4=(size(H4,2)-1)/2;
     else
          length_calc4=0;
     end
     if length_calc4>10
         c4=[];
     for j=1:length_calc4         %load all coordinates
         I4=getfield(H4,{1,2*j},'Attributes');
         X4=getfield(I4,{1,1},'Value');
         Y4=getfield(I4,{1,2},'Value');
         c4(j,:)=[str2num(X4),str2num(Y4)];
     end
     c4(length_calc4+1,:)=c4(1,:); %close contour
     c4=c4*calibrationfactor;
     polyin2 = polyshape({lum(:,1)',c(:,1)',c2(:,1)',c3(:,1)',c4(:,1)'},{lum(:,2)',c(:,2)',c2(:,2)',c3(:,2)',c4(:,2)'}, 'simplify', false);
     [x,y] = centroid(polyin2);
       
     calc_x4=c4(:,1);
     calc_y4=c4(:,2);     
     calc_z4=linspace(image_spacing*i-image_spacing,image_spacing*i-image_spacing,length(calc_x4));

     calc_x_new4=[calc_x_new4;(calc_x4)];   %move back to original
     calc_y_new4=[calc_y_new4;(calc_y4)];   %move back to original
     calc_z_new4=[calc_z_new4;(calc_z4)'];
     
     A_c4=polyarea(c4(:,1)',c4(:,2)');
     else
     c4=[];
     A_c4=0;
     end

%5th calc contour (error)
     if size(G,2)>9
         error=1;
     end

else                            %do not calculate calc if file does not exist
     c=[];
     c2=[];
     c3=[];
     c4=[];
     polyin = polyshape({lum(:,1)'},{lum(:,2)'}, 'simplify', false);
     [x,y] = centroid(polyin);
     A_c=0;
     A_c2=0;
     A_c3=0; 
     A_c4=0;
     A_lum=polyarea(lum(:,1)',lum(:,2)');
end

P=[lum;c;c2];                   %combine all points
k=convhull(P);                  %calculate index of convex hull points

%% vessel contour
d_mean=sqrt((A_lum+A_c+A_c2)*4/pi);  %mean diameter of total area (lumen + calc)
t=-ratio/(2*ratio-1)*d_mean;         %thickness of wall
D_mean=2*t+d_mean;                   %mean diameter of vessel contour
lum_shift=[P(k,1)-x,P(k,2)-y];       %move to centroid

[lum_shift_phi,lum_shift_r]=cart2pol(lum_shift(:,1)',lum_shift(:,2)'); %polar coords
ves_r_o=lum_shift_r+t;                                                 %convex hull + mean thickness t
%ves_r instead of ves_r_o if you don't want to interpolate

%interpolating
lum_shift_phi(end)=[];
ves_r_o(end)=[];
Test=[lum_shift_phi',ves_r_o'];
Test=sortrows(Test);

angle=linspace(min(Test(:,1)),max(Test(:,1)),100);
ves_r=interp1(lum_shift_phi,ves_r_o,angle);
ves_r(end+1)=ves_r(1);
angle(end+1)=angle(1);
                      %lum_shift_phi instead of angle if you don't want to interpolate
[ves_x,ves_y]=pol2cart(angle,ves_r);                           %cartesian coords
ves_z=linspace(image_spacing*i-image_spacing,image_spacing*i-image_spacing,length(ves_x));
if i==1
    ves_x_new=[];
    ves_y_new=[];
    ves_z_new=[];
end
ves_x_new=[ves_x_new;(ves_x+x)'];   %move back to original
ves_y_new=[ves_y_new;(ves_y+y)'];   %move back to original
ves_z_new=[ves_z_new;(ves_z)'];

%% progress
progress=round(100*(filenum(i)-filenum(1))/(2*length(filenum)));
fprintf('Progress: %d percent\n',progress);
close all
end
%%
shp = alphaShape(ves_x_new,ves_z_new,-ves_y_new,threshold,'HoleThreshold',1E10); %change csys
plot(shp)
[F,V] = boundaryFacets(shp);                                             %Faces and Coordinates

%% save wall
i=1;
fprintf('Create images\n');
for i=1:length(filenum)

    YPlane=i*image_spacing-image_spacing;
    [VV,EE]=Polyhedron_Plane(V,F,YPlane);

    vesVV_x=VV(:,1)-5;
    vesVV_y=-VV(:,3)-5;
    [theta,rho] = cart2pol(vesVV_x,vesVV_y);
    O=[theta,rho];
    O=sortrows(O);
    [vesVV_x_sorted,vesVV_y_sorted] = pol2cart(O(:,1),O(:,2));
    vesVV_x_sorted=vesVV_x_sorted+5;
    vesVV_x_sorted(end+1)=vesVV_x_sorted(1);
    vesVV_y_sorted=vesVV_y_sorted+5;
    vesVV_y_sorted(end+1)=vesVV_y_sorted(1);

    %lumen
    b=find(lum_z_new==double(image_spacing*i-image_spacing));
    l_x=lum_x_new(b);
    l_y=lum_y_new(b);
    
    %create image
    I = zeros(1000, 1000);
    BW_wall = roipoly(I,100*vesVV_x_sorted,100*vesVV_y_sorted);
    %BW_wall = imcomplement(BW_wall);
    
windowSize=100;  % Decide smoothing contour
kernel=ones(windowSize)/windowSize^2;
result=conv2(single(BW_wall),kernel,'same');
result=result>0.5;

BW_wall_smooth = imcomplement(result); 
    
    BW_lum = roipoly(I,100*l_x,100*l_y);
    BW = imadd(BW_wall_smooth,BW_lum);
    h=imshow(BW);

    %save image
    cd(path_image_folder)
    filename = ['Wall_' num2str(filenum(i)) '.bmp'];
    imwrite(BW, filename, 'bmp')

%% save calc
    d=find(calc_z_new==double(image_spacing*i-image_spacing));
    c_x=calc_x_new(d);
    c_y=calc_y_new(d);
    
    d2=find(calc_z_new2==double(image_spacing*i-image_spacing));
    c_x2=calc_x_new2(d2);
    c_y2=calc_y_new2(d2);  
    
    d3=find(calc_z_new3==double(image_spacing*i-image_spacing));
    c_x3=calc_x_new3(d3);
    c_y3=calc_y_new3(d3); 
    
    d4=find(calc_z_new4==double(image_spacing*i-image_spacing));
    c_x4=calc_x_new4(d4);
    c_y4=calc_y_new4(d4); 
    
if isempty(d)==1
     BW_blank = imcomplement(I);

     h=imshow(BW_blank);

     cd(path_image_folder)
     filename = ['Calc_' num2str(filenum(i)) '.bmp'];
     imwrite(BW_blank, filename, 'bmp')
else
     BW_calc = roipoly(I,100*c_x,100*c_y);
     BW_calc = imcomplement(BW_calc);
     
     BW_calc2 = roipoly(I,100*c_x2,100*c_y2);
     BW_calc2 = imcomplement(BW_calc2);
     
     BW_calc3 = roipoly(I,100*c_x3,100*c_y3);
     BW_calc3 = imcomplement(BW_calc3);
     
     BW_calc4 = roipoly(I,100*c_x4,100*c_y4);
     BW_calc4 = imcomplement(BW_calc4);
     
     BW_calc10=imfuse(BW_calc,BW_calc2);      %fuse calc 1&2
     BW_calc10=rgb2gray(BW_calc10);
     BW_calc10=im2bw(BW_calc10,0.9);
     
     BW_calc11=imfuse(BW_calc10,BW_calc3);    %fuse calc 3
     BW_calc11=rgb2gray(BW_calc11);
     BW_calc11=im2bw(BW_calc11,0.9);  
     
     BW_calc12=imfuse(BW_calc11,BW_calc4);    %fuse calc 4
     BW_calc12=rgb2gray(BW_calc12);
     BW_calc12=im2bw(BW_calc12,0.9);  
     
     h=imshow(BW_calc12);
     
     %save image     
     cd(path_image_folder)
     filename = ['Calc_' num2str(filenum(i)) '.bmp'];
     imwrite(BW_calc12, filename, 'bmp')
end
end

%% save unsegmented
i=1;
filenum_new=filenum(1)+1:2:filenum(end)-1;
fprintf('Create blank images for interpolation\n');
for i=1:length(filenum_new)
     BW_blank = imcomplement(I);

     h=imshow(BW_blank);

     %save image 
     cd(path_image_folder)
     filename = ['Calc_' num2str(filenum_new(i)) '.bmp'];
     imwrite(BW_blank, filename, 'bmp')
     
     filename = ['Wall_' num2str(filenum_new(i)) '.bmp'];
     imwrite(BW_blank, filename, 'bmp')
end
%%
fprintf('Progress: 100 percent\n');
if error==1
  fprintf('More than 4 calcifications detected\n');
end
