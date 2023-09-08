clear all;
close all;
global filename horz_spacing vert_spacing backgroundWidth num_peaks snr_threshold

%%
%%Image rotation and region captures
filename = '0001_532.tif';
horz_spacing=60;
vert_spacing=200;
[struct1]=roiGeneration(filename,horz_spacing,vert_spacing);

%%
%%background substraction

backgroundWidth=50;
[struct2]=intProf(struct1,backgroundWidth);

%%
%%fitpeaks
num_peaks=1;
snr_threshold=1;
[struct3]=fitPeaks(struct2, num_peaks, snr_threshold);

%% 
%%Imageplotting
[struct4] = goodProfiles(struct3);
AUC=struct4.AUC;
AUC=AUC(:);
dev_to_analyze = struct4.index_dev_to_analyze; % list of  right peaks
center=struct4.center(:); %peak center
width=struct4.width(:);%peak width
sigma=struct4.sigma(:);%peak sigma

AUC_after = AUC(dev_to_analyze);
center_after = center(dev_to_analyze);
width_after = width(dev_to_analyze);
sigma_after = sigma(dev_to_analyze);

figure (5)
histogram(AUC_after,'FaceColor','g');
figure (6)
histogram(center,'FaceColor','r')
figure (7)
scatter(center,AUC)
figure (8)
histogram(width,'FaceColor','b')

iptwindowalign(figure (5),'right',figure (6),'left')


%% Data Saving
%% EpCAM
filename='Data.xlsx';
protein_type={'EpCAM'};
col_header={'center','width','sigma','AUC','center_after','width_after','sigma_after','AUC_after'};
writecell(protein_type,filename,'Sheet','Saving','Range','B1');
writecell(col_header,filename,'Sheet','Saving','Range','B2');
writematrix(center,filename,'Sheet','Saving','Range','B3');
writematrix(width,filename,'Sheet','Saving','Range','C3');
writematrix(sigma,filename,'Sheet','Saving','Range','D3');
writematrix(AUC,filename,'Sheet','Saving','Range','E3');
writematrix(center_after,filename,'Sheet','Saving','Range','F3');
writematrix(width_after,filename,'Sheet','Saving','Range','G3');
writematrix(sigma_after,filename,'Sheet','Saving','Range','H3');
writematrix(AUC_after,filename,'Sheet','Saving','Range','I3');

% %% Vimentin
% filename='Data.xlsx';
% protein_type={'Vimentin'};
% col_header={'center','width','sigma','AUC'};
% writecell(protein_type,filename,'Sheet','MCF7_30s','Range','G1');
% writecell(col_header,filename,'Sheet','MCF7_30s','Range','G2');
% writematrix(center,filename,'Sheet','MCF7_30s','Range','G3');
% writematrix(width,filename,'Sheet','MCF7_30s','Range','H3');
% writematrix(sigma,filename,'Sheet','MCF7_30s','Range','I3');
% writematrix(AUC,filename,'Sheet','MCF7_30s','Range','J3');

% %% GAPDH
% filename='Data.xlsx';
% protein_type={'GAPDH'};
% col_header={'center','width','sigma','AUC'};
% writecell(protein_type,filename,'Sheet','MCF7_30s','Range','L1');
% writecell(col_header,filename,'Sheet','MCF7_30s','Range','L2');
% writematrix(center,filename,'Sheet','MCF7_30s','Range','L3');
% writematrix(width,filename,'Sheet','MCF7_30s','Range','M3');
% writematrix(sigma,filename,'Sheet','MCF7_30s','Range','N3');
% writematrix(AUC,filename,'Sheet','MCF7_30s','Range','O3');

% %% Her2
% filename='Data.xlsx';
% protein_type={'Her2'};
% col_header={'center','width','sigma','AUC'};
% writecell(protein_type,filename,'Sheet','MCF7_30s','Range','Q1');
% writecell(col_header,filename,'Sheet','MCF7_30s','Range','Q2');
% writematrix(center,filename,'Sheet','MCF7_30s','Range','Q3');
% writematrix(width,filename,'Sheet','MCF7_30s','Range','R3');
% writematrix(sigma,filename,'Sheet','MCF7_30s','Range','S3');
% writematrix(AUC,filename,'Sheet','MCF7_30s','Range','T3');

%%
function [struct] = roiGeneration(filename,horzspacing,vertspacing,struct)
global filename horz_spacing vert_spacing
% This function rotates and aligns a raw fluorescent image of a single-cell
% Western blot array and segments the image into regions-of-interest (ROIs)
% for downstream analysis. Each region of interest encompasses single-cell Western blot protein peak(s) in an area
% defined by the horizontal and vertical spacing between microwells in the
% array.

%   Outputs:  
% Struct [structure]: A data structure containing objects:
%   struct.rois: 3D matrix with each ROI contained in a different z.
%   struct.angle: The angle of rotation to straighten the image (number, in degrees).
%   struct.rotate: The angle of rotation required to display the image with
%   separations running vertically instead of horizontally (number, in
%   degrees).
%   struct.array_bounds: User selected boundaries of the array as a 3x2
%   matrix (rows contain upper left, upper right, and lower left
%   coordinates respectively; first column contains x-coordinates; second column contains y-coordinates). 
%   struct.name: The name of the protein target entered by the user
%   (string).
%   struct.wells_per_row: The number of wells per row based on the user
%   selected array bounds and horizontal well spacing.
%   struct.rows: Number of rows in the array 

%   Inputs:
% filename [string]: A string containing the name of the fluorescence image
%                   to be processed.
% horzspacing [num]: Well-to-well spacing (horizontal, in pixels)
% vertspacing [num]: Well-to-well spacing (vertical, in pixels)
% struct [structure] (optional): A structure containing "angle" and
% "array_bounds" if the same image has already been analyzed by
% roiGeneration. The same ROIs will automatically be generated.
%% versions
% 0.1-Created April, 2016
% 0.2 (5.15.16): Updated to apply same transform for ROI generation if user
% inputs a struct with the fields "angle" and "rotate".
%0.3 (5.20.16): Added "rows" and "wells per row" fields to structure.

%% Check input arguments
switch nargin
    % If the user only provides the image, horizontal and vertical spacing
    case 3
        transform = 0;
    case 4
        transform = 1;
        
        tf = isstruct(struct);
            if tf == 0
                 error('Input argument "struct" is not a structure.');
            
            return
            
            end
            
        % retrieve previously determined angle for transformation of image
        angle = struct.angle;
        
        % retrieve previously determined array boundaries
        array_bounds = struct.array_bounds;
        
        % extract the individual x an y coordinates of the array boundaries
        x_upperleftwell = array_bounds(1, 1);
        y_upperleftwell = array_bounds(1, 2);
        
        x_upperrightwell = array_bounds(2, 1);
        y_upperrightwell = array_bounds(2, 2);
        
        x_lowerrightwell = array_bounds(3, 1);
        y_lowerrightwell = array_bounds(3, 2);
        
    
    otherwise
        
        error('Invalid number of input arguments');
            
        return
        
end
%% 
% ask the user the name of their protein target
prompt = 'What is the name of your protein target? ';
str = input(prompt, 's');
struct.name = str;


% Load the image file in MATLAB
img = imread(filename);
    
    if transform == 0
        
        % Display more contrasted image in window
        contrasted_img = histeq(img);
        imshow(contrasted_img);

        % Display a message to the user asking them to look at the array
        title('Take a look at the array and determine if the wells are oriented left of the bands or right of the bands. Then press any key');
        pause()

        % Construct a questdlg to ask the user how the image is currently oriented
        % for coarse rotation
        choice = questdlg('Are the wells currently left of the bands or right of the bands?', ...
        'Current array orientation', ...
        'Wells are left of bands','Wells are right of bands','Wells are right of bands');
       
        % Handle response
        switch choice
            
            case 'Wells are left of bands'; 
                disp([choice 'Okay, the image will be rotated to the right!'])
                rotate = -90;
                
            case 'Wells are right of bands';
                disp([choice 'Okay, the image will be rotated to the left!'])
                rotate = 90;
        end
        
        % Store the course rotation angle to orient the array vertically to
        % the struct
        
        struct.rotate = rotate;
    else
         
        rotate = struct.rotate;
    end
  
  % Display the course-rotated image
  imgrotated = imrotate(img, rotate);
  contrasted_img_r = histeq(imgrotated);
  imshow(contrasted_img_r);
  
  %If struct was not an input argument (and there is no previous
  %angle/array boundary values to draw from), the user will now manually
  %select the array boundaries.
  
while transform == 0
  test = 1;
    
  while test == 1
        % Prompt user to select the upper right well of the array. 
        title('Please zoom in on the the middle of the upper left well and press any key.');
        
        % use mouse button to zoom in or out
        zoom on;   
        pause()
        zoom off;
        
        % preallocate array bounds matrix
        array_bounds = zeros(3, 2);
        
        % prompt user to click on the middle of the upper left well
        title('Please click on the middle of the upper left well.');
        
        [x_click,y_click] = ginput(1);
        
        % store the coordinates the user selected for the upper left well
        x_upperleftwell = x_click;
        y_upperleftwell = y_click;
        zoom out;
        
        array_bounds(1,:) = [x_upperleftwell, y_upperleftwell];
        
        % Change message displayed in figure window to indicate the user should zoom in on the
        % upper right well
        title('Please zoom in on the middle of the upper right well and press any key.')
        
        % use mouse button to zoom in or out
        zoom on;   
        pause()
    
        zoom off;
        
        % prompt user to click on the middle of the upper right well
        title('Please click on the middle of the upper right well.');
        
        [x_click,y_click] = ginput(1);
        
        % store the coordinates of the user-selected upper right well
        x_upperrightwell = x_click;
        y_upperrightwell = y_click;
        zoom out;
        
        array_bounds(2,:) = [x_upperrightwell, y_upperrightwell];
        
        % Change display in imaged window to indicate user should zoom in
        % on the middle of the lower right well
        title('Please zoom in on the middle of the lower right well and press any key.')
        
        % use mouse button to zoom in or out
        zoom on;   
        pause()
        zoom off;
    
        % prompt user to click on the middle of the lower right well
        title('Please click on the middle of the lower right well.');
        
        [x_click,y_click] = ginput(1);
        
        % store the user-selected coordinates of the lower right well
        x_lowerrightwell = x_click;
        y_lowerrightwell = y_click;
        
        array_bounds(3,:) = [x_lowerrightwell, y_lowerrightwell];
        
        % store all of the coordinates of the array bounds to the struct
        struct.array_bounds = array_bounds;
        
        % Construct a questdlg to ask the user if they are happy with their
         % well selection
        choice = questdlg('Are you happy with your well selections?', ...
        'Well selections for array boundaries', ...
        'Yes','No','Yes');
        
        % Handle response
        switch choice
            
            case 'Yes';
                disp([choice 'Great, let''s keep going then!'])
                test = 0;
            
            case 'No';
                disp([choice 'That''s okay, try again!'])
                test = 1;
        end
     
    % check whether the user selected array boundaries are correct    
    if (x_upperrightwell<x_upperleftwell || y_upperrightwell>y_lowerrightwell)        
        test = 1;
        
        title('Oh no! We detected you selected the wells in the wrong order. Please try again. Press any key to continue')
        pause()
    else
        test = 0;
    end
  end
  
    % store the coordinates of the direction vector that extends from the upper left well to the right most point of the array
    dir_vector1 = [x_upperrightwell,y_upperleftwell] - [x_upperleftwell,y_upperleftwell];

    % store the coordinates of the direction vector that extends from the upper left well to the upper right well 
    dir_vector2 = [x_upperrightwell,y_upperrightwell] - [x_upperleftwell,y_upperleftwell];

    % Find angle between the two direction vectors [angle in degrees]
    cosangle = dot(dir_vector1, dir_vector2) / (norm(dir_vector1) * norm(dir_vector2));
    angle = acosd(cosangle);
    
        if (y_upperrightwell<y_upperleftwell)
            angle=-angle;
        end
    
    % store the angle used to straigten the image in the struct
    struct.angle=angle;  
    transform=1;
end


% Display the rotated image so the array is aligned
b = imrotate(imgrotated, angle, 'nearest','crop');
b_contrasted = histeq(b);
imshow(b_contrasted);
hold on
sz = size(b) / 2;

% Generate a rotation matrix to multiply by the array boundary coordinates
% to attain the new array boundaries in the rotated image
rotation_matrix = [cosd(-angle), -sind(-angle);sind(-angle), cosd(-angle)];

% Multiply the rotation matrix by the upper left well coordinates
new_upper_left = rotation_matrix * [(x_upperleftwell - (sz(2)));(y_upperleftwell - sz(1))];

% Multiply the rotation matrix by the upper right well coordinates
new_upper_right = rotation_matrix * [(x_upperrightwell - sz(2));(y_upperrightwell - sz(1))];

%Multiply the rotation matrix by the lower right well coordinates
new_lower_right = rotation_matrix * [(x_lowerrightwell - sz(2));(y_lowerrightwell - sz(1))];

% store the new upper left x and y coordinates
x_new_upper_left = new_upper_left(1) + sz(2);
y_new_upper_left = new_upper_left(2) + sz(1);

% store the new upper right x and y coordinates
x_new_upper_right = new_upper_right(1) + sz(2);
y_new_upper_right = new_upper_right(2) + sz(1);

% store the new lower right x and y coordinates
x_new_lower_right = new_lower_right(1) + sz(2);
y_new_lower_right = new_lower_right(2) + sz(1);


% Determine number of wells per row
wells_per_row = round((x_new_upper_right - x_new_upper_left) / horzspacing);
struct.wells_per_row = wells_per_row;

% Determine number of rows
rows = round((y_new_lower_right - y_new_upper_right) / vertspacing);
struct.rows = rows;

% Determine total number of wells
total_wells = wells_per_row * rows;


% for loop to fill in the 3D matrix with ROIs from the image (proceeds row by row of the microwell array from left to right)
% pre-allocate 3D matrix with zeros
mat = zeros(vertspacing, horzspacing, total_wells);

for i = 1:rows
    for j = 1:wells_per_row
        
        % determine z-coordinate for the current ROI
        z = (wells_per_row) * (i-1)+j;
        
        % set row start and end boundaries 
        row_start = (round(x_new_upper_left) - horzspacing/2) + ((j-1)*horzspacing);
        row_end = row_start + horzspacing;
        
        % set column start and end boundaries
        col_start = (round(y_new_upper_left) + ((i-1)*vertspacing));
        col_end = col_start + vertspacing;
        
        %generate lines that span the x and y coordinates of all the ROIs
        %to overlay over image to show the ROIs
        x = row_start:1:(row_end - 1);
        y = repmat(col_start, 1, length(x));
        y2 = col_start:1:(col_end - 1);
        x2 = repmat((row_end-1), 1, length(y2));
        
        % fill the matrix with the image pixels within the current ROI
        % boundaries
        mat(: ,: ,z) = b(col_start:(col_end - 1), row_start:(row_end - 1));
        
        % plot the ROI grid overlay on the image
        figure(1)
        plot(x', y', 'Color', 'w', 'LineStyle','-');
        plot(x', y', 'Color', 'k', 'LineStyle',':');
        plot(x2', y2', 'Color', 'w', 'LineStyle','-');
        plot(x2', y2', 'Color', 'k', 'LineStyle',':');
    end
end
saveas(gcf,'figure(1)');
% store the 3D matrix of ROIs to the struct
struct.rois = mat;
end

function [struct] = intProf(struct,backgroundwidth)
global backgroundWidth
% Generate intensity profiles from the ROI stacks in the output of roiGeneration
% and perform background subtraction on the profiles
% 
% Outputs
% Struct [structure]: A data structure containing objects (Intensity 
%                     profiles for each ROI, 3D matrix with each ROI 
%                     contained in a different z, and coordinates of 
%                     the ROIs).
%   struct.rois: 3D matrix with each ROI contained in a different z.
%   struct.angle: The angle of rotation to straighten the image (number, in degrees).
%   struct.rotate: The angle of rotation required to display the image with
%   separations running vertically instead of horizontally (number, in
%   degrees).
%   struct.array_bounds: User selected boundaries of the array as a 3x2
%   matrix (rows contain upper left, upper right, and lower left
%   coordinates respectively; first column contains x-coordinates; second column contains y-coordinates). 
%   struct.name: The name of the protein target entered by the user
%   (string).
%   struct.wells_per_row: The number of wells per row based on the user
%   selected array bounds and horizontal well spacing.
%   struct.rows: Number of rows in the array
%   struct.int_prof: a 3D matrix containing an array of the intensity 
%                    profiles [x, intensity value] indexed by the 
%                    third dimension
%                    
% 
% Inputs
%   Struct [structure]: The data structure from roiGeneration containing the
%                       ROIs for each lane
%   backgroundWidth [int]: width of the background region for axial
%                          background subtraction (in pixels)
%                          

% Load the matrix of ROIs

mat=struct.rois;

% Determine the number and size of the intensity profiles
[x_dim,y_dim,z_dim]=size(mat);

% Preallocate the intensity profile matrix with zeros

int_profiles=zeros(x_dim,2,z_dim);

% pix_conversion is the number of microns per pixel
pix_conversion=5;

% for loop to generate the intensity profiles for each ROI in the z-stack 
%of the matrix Mat. The intensity profile is an average of the pixel intensities 
%across the short-axis of the ROI. The background regions are defined by the
%parameter backgroundWidth, and the average pixel intensity in the left and
%right background regions are calculated. The background subtracted intensity
%profile is generated by subtracting the mean background intensity at each 
%point along the long-axis of the ROI from the average pixel intensity.
figure (2)
 for i=1:z_dim
    
    % Get the image of the lane
    lane = mat(:,:,i);
    
    % Get a list of the x coordinates
    dist = (0:pix_conversion:pix_conversion*(x_dim-1));
    
    % Sum the image along the y-axis (transverse to separation axis) to
    % generate an intensity profile
    int = sum(lane,2);
    avg_int = int / (y_dim);
    
    % Get the background regions to the left and right of the lane
    left_backgroundregion = lane(:, (1:backgroundwidth)); 
    right_backgroundregion = lane(:, (((end + 1) - backgroundwidth):end));
    
    % Calculate average background intensity in the left and right regions
    left_background_int = (sum(left_backgroundregion, 2)) / backgroundwidth;
    right_background_int = sum(right_backgroundregion, 2)/ backgroundwidth;
    
    % Create a vector containing the background intensity for each
    % x-coordinate along the lane
    mean_background = (left_background_int + right_background_int) / 2;
    
    % Subtract the background vector from the intensity profile
    bsub_int = avg_int - mean_background;
    
    % Create a matrix with one column containing the x-coordinates and a
    % second column containing the background subtracted intensity profile
    lane_profile=[dist',bsub_int];
    
    % Add the intensity profile to the matrix of intensity profiles
    int_profiles(:,:,i) = lane_profile;
    
    % Add the intensity profile to the plot
    figure (2) 
    plot(dist',bsub_int);
    
    hold on
 end

% Save the intensity profiles to the data structure
struct.int_prof = int_profiles;



end


function data_struct = fitPeaks(data_struct, num_peaks, snr_threshold)
   global num_peaks snr_threshold
%% Check input arguments
switch nargin
    
    % If only the data_structure is provided, set num_peaks = 1
    case 1
        
        num_peaks = 1;
        
    % If provided, ensure the number of peaks is valid
    case 2
        % Exit function if an invalid number of peaks is input
        if ((num_peaks > 3) || (num_peaks < 1))
            
            error('Invalid number of peaks');
            
            return
            
        end
        %If only 2 input arguments provided, user does not want to run the
        %SNR threshold
        apply_snr_threshold=0;
    case 3
        apply_snr_threshold=1;
        %
    otherwise
        
        error('Invalid number of input arguments');
            
        return
    
    
end


%% Get the peaks

% Get the intensity profiles
try 
    intensity_profiles = data_struct.int_prof;

catch 
    
    error('Error accessing data_struct.int_prof');
    
end

% Find all of the good wells
[x_dim,y_dim,z_dim]=size(intensity_profiles);

%store to the structure the starting number of wells analyzed
struct.total_wells=z_dim;


if apply_snr_threshold==1    
    % for loop to filter out SNR<3 lanes with a conservative SNR estimate 
    %calculated from the max intensity of a smooth data set and the standard 
    %deviation of the last 5 pixels of the lane.
    snr3_devices=zeros(z_dim,1);
    snr_est=zeros(z_dim,1);

    %figure
        for i=1:z_dim
            device=intensity_profiles(:,:,i);
            xval=device(:,1);
            yval=device(:,2);
            yvalsmooth=smooth(yval);
	
            noise_est=std(yval(end-5:end));
	
            signal_est=max(yvalsmooth);
	
            snr_est(i)=signal_est/noise_est;
    
                if snr_est(i)<snr_threshold
                    snr3_devices(i)=0;
                else
                    snr3_devices(i)=1;
                    %plot(xval,yval);
                    %hold on
                end
        end
    struct.snr_est=snr_est;
    
    % Get the number of good wells
    num_good_devices = sum(snr3_devices);
    good_indices=find(snr3_devices==1);
    % Exit if there are no good wells
    if (num_good_devices == 0)
    
        error('No good wells in data_struct');
    
    end
else
    num_good_devices=ones(z_dim,1);
    good_indices=find(num_good_devices==1);
end

% Save the good indices
data_struct.good_indices = good_indices;

%% Get the seed parameters



% Let the user select the points for the parameter estimation

bounds_set = false;

while (~bounds_set)
    
    % Plot the good devices
    figure(3);
    hold on

    for i = 1:num_good_devices

        device_index = good_indices(i);

        plot(intensity_profiles(:,1, device_index),...
            intensity_profiles(:,2, device_index), '-k');

    end

%     hold off
    
    uiwait(msgbox('Please select left and right boundaries of each peak'));
    
    % Get the limits of the plot
    y_lim = get(gca, 'YLim');
    
    % Preallocate the nx2 matrix to hold the peak bounds, where n is the
    % number of peaks. Col 1 is the left bound, col 2 is the right bound
    peak_bounds = zeros(num_peaks, 2);
    
    for peak = 1:num_peaks
       
        % Get the left peak boundary
        [x1, y1] = ginput(1);
        
        % Draw the selected peak boundary
        line([x1, x1], y_lim, 'Color', [0, 0, 1]);
        
        
        % Get the right peak boundary
        [x2, y2] = ginput(1);
        
        % Draw the selected peak boundary
        line([x2, x2], y_lim, 'Color', [0, 1, 0]);
        
        % Save the selected bounds
        
        peak_bounds(peak, :) = [x1, x2];
        
          
    end
    
    % Ask if the peaks are correct
    choice = questdlg('Are the peak bounds correct?', ...
	'Done with bound selection?', ...
	'Yes', 'No','No');

    close
   
    % If they are done, exit loop
    if (strcmp(choice, 'Yes'))
       
        bounds_set = true;
        
    end
    
    

end


% Save the bounds
data_struct.fit_bounds = peak_bounds;

%% Create the fit options

% Create the fit options object with the specified number of peaks
switch num_peaks
    
    case 1
        fit_type = 'gauss1';
         
    case 2
        fit_type = 'gauss2';
        
    case 3
        fit_type = 'gauss3';

end

fit_options = fitoptions(fit_type);  

% Assign the locations to the fit options object
for peak = 1:num_peaks
    
    % Get the left and right bound for the peak
    left_bound = peak_bounds(peak, 1);
    right_bound = peak_bounds(peak, 2);
    
    % Set the sigma bounds
    sigma_min = 0;
    sigma_max = right_bound - left_bound;
    
    % Set the peak center bounds
    x_min = left_bound;
    x_max = right_bound;
    
    % Set the ampitude bounds
    a_min = 0;
    a_max = y_lim(2);
    
    % set the upper and lower bounds. correct for difference in c and
    % sigma terms
    
    lb = peak*3 - 2;
    ub = peak*3;

    fit_options.Lower(lb:ub) = [a_min, x_min, (sigma_min * sqrt(2))];
    fit_options.Upper(lb:ub) = [a_max, x_max, (sigma_max * sqrt(2))];

    %% fit_options.Lower = [a_min, x_min, (sigma_min * sqrt(2))];
    %% fit_options.Upper = [a_max, x_max, (sigma_max * sqrt(2))];
    
end


%% Fit each peak

% Preallocate the m x 3 x p matrix for the fit coefficients were m is the 
% number of peaks per roi and p is the number of ROIs and col 1 is the 
% amplitude, col 2 is the peak center, and col 3 is sigma
data_struct.fit_coefficients = zeros(num_peaks, 3, length(good_indices));

% Preallocate the m x 1 matrix for the R^2 values for each fit where
% m is the number of good devices
data_struct.R2 = zeros(num_good_devices, 1);

for i = 1:num_good_devices 
    
    device_index = good_indices(i);
    
         %Display the device number every 50 devices
    %if(mod(i, 10) == 0)
        
        %fprintf('Fitting lane %d/%d\n', i, num_good_devices);
        
    %end
    
    % Get the x and y values
    x = intensity_profiles(:,1, device_index);
    y = intensity_profiles(:,2, device_index);
    
    %Determine index of x_min and x_max for selection of x and y values in
    %the region of the peak
    left_diff=abs(x-x_min);
    left_data=find(left_diff==min(left_diff));
    
    right_diff=abs(x-x_max);
    right_data=find(right_diff==min(right_diff));
    
    %Get the x and y values in the peak region
    x_fit=x(left_data:right_data);
    y_fit=y(left_data:right_data);
    
    % Fit the peaks
    
    [fit_object, gof] = fit(x_fit, y_fit, fit_type, fit_options);
    
    % Get the coefficients
    fit_coeffs = coeffvalues(fit_object);
    
    % Save the coefficients
    for peak = 1:num_peaks
        
        coeff_index_start = (peak - 1)*3 + 1;
        coeff_index_end = peak * 3;
        
        %get peak center and width for AUC calculation
        center=fit_coeffs(coeff_index_end-1);
        sigma=fit_coeffs(coeff_index_end);
        width=sigma/sqrt(2);
        
        %determine location of +/- 2 peak widths from the peak center
        auc_left_bound=center-2*width;
        auc_right_bound=center+2*width;
        
        %Determine index of auc_left_bound and auc_right_bound for selection of x and y values in
        %the region of the peak
        left_diff_auc=abs(x-auc_left_bound);
        left_data_auc=find(left_diff_auc==min(left_diff_auc));
    
        right_diff_auc=abs(x-auc_right_bound);
        right_data_auc=find(right_diff_auc==min(right_diff_auc));
        
        % Make sure the left bound is within the array
        if (left_data_auc < 1)
           
           left_data_auc = 1; 
            
        end
        
        
        % Check to make sure the AUC bounds are within the bounds of the
        % array
        if (right_data_auc > length(y))
            
            right_data_auc = length(y);
            
            
        end
    
        %Sum data within the peak bounds
        peak_region_intensities=y(left_data_auc:right_data_auc);
        AUC(peak,1,i)=sum(peak_region_intensities);
  
        data_struct.fit_coefficients(peak, :, i) =...
            fit_coeffs(coeff_index_start:coeff_index_end);
        data_struct.AUC(peak,1,i)=AUC(peak,1,i);
        data_struct.width(peak,1,i)=width;
        data_struct.sigma(peak,1,i)=sigma;
        data_struct.center(peak,1,i)=center;
    end
    
    % Save the R^2 value
    data_struct.R2(i) = gof.rsquare;
    
end

end


function [data_struct] = goodProfiles(data_struct,r2_threshold)
%perform quality control on intensity profiles, removing lanes with SNR<3
%and allowing the user to select lanes to remove upon visual inspection
%   Outputs
%Struct [structure]: A data structure containing objects (indices of “good” lanes, 
%intensity profiles for each ROI, 3D matrix with each ROI contained in a 
%different z, and coordinates of the ROIs)
%Inputs
%struct [structure]: A structure containing the intensity profiles generated in intProf

% v04

%% Check input arguments
switch nargin
    
    % If only the data_structure is provided, set r2_threshold = 0.7
    case 1
        
        r2_threshold = 0.7;
        
    % If provided, ensure the r2 value is valid
    case 2
        % Exit function if an invalid r2 value is input
        if ((r2_threshold<0) || (r2_threshold > 1))
            
            error('Invalid R^2 value');
            
            return
            
        end
end       

int_prof_all=data_struct.int_prof;
[x_dim,y_dim,z_dim]=size(int_prof_all);
r2=data_struct.R2;
good_r2=find(r2>=r2_threshold);

%determine array position of high r2 value fit lanes
good_indices=data_struct.good_indices;
good_fits=good_indices(good_r2);

good_int_profiles=zeros(x_dim,y_dim,length(good_r2));

for i=1:length(good_r2)
    good_int_profiles(:,:,i)=int_prof_all(:,:,good_fits(i));
end

%set number of rows/columns of subplots to display in each figure window
n=5;
num_subplots=n*n;

plots_display=length(good_r2);
good_devices=ones(length(good_r2),1);
number_subplots=ceil(plots_display/(n*n));
dev_to_analyze=zeros(z_dim,1);

good_subplots = ones(plots_display,1);

disp(number_subplots);

[xfit,yfit,zfit] = size(data_struct.fit_coefficients);

colors = {'-r','-g'};

% for loop to generate subplots for user inspection of the intensity profiles
for i=1:number_subplots
    
    disp(i);
    
    figure
    if i==1
        if (length(good_r2) > (n*n))
          devices_subplot=(1:(n*n));
        else
          devices_subplot=(1:length(good_r2));
        end
    elseif i*n*n<=plots_display
        devices_subplot=((i*n*n)-(n*n)+1):((i*n*n));
    else 
        devices_subplot=((i*n*n)-(n*n)):(plots_display);
    end
    next=0;
    plot_colors = zeros(length(devices_subplot),2);

	for j=1:length(devices_subplot)
             
            dev_number=devices_subplot(j);
            device=good_int_profiles(:,:,dev_number);
            xval=device(:,1);
            yval=device(:,2);
            subplot(n, n, j);
            hold all
            
            % Plotting fit
            fit_extract = data_struct.fit_coefficients(1,:,good_r2(dev_number));
            gauss_y = gaussReturn(fit_extract(1),fit_extract(2),fit_extract(3),xval);


            plot(xval,yval,'Color',[0 0 1],'LineWidth',2,'Tag', sprintf('%d', dev_number), 'buttondownfcn',@clickTest);
            plot(xval,gauss_y,'Color',[0 1 1 0.9],'LineWidth',2);

            hold off
    end
    next=0;
    
    btn = uicontrol('Style', 'pushbutton', 'String', 'Next',...
        'Position', [500 15 50 30],...
        'Callback',@continueButton2);

    while next==0   

        pause(0.01);
    end
    
    %good_devices(devices_subplot(1):devices_subplot(end))=good_subplots;
    
    
close(gcf)    
end 

data_struct.dev_to_analyze = good_subplots;
good_subplot_ind=find(good_subplots==1);
data_struct.index_dev_to_analyze=good_r2(good_subplot_ind);
end
 

function [next]=continueButton2(qstring,title,str1,str2,default)
%UNTITLED5 Summary of this function goes here
qstring='Are you done selecting devices to throw out?';
title='Device Quality Control';
str1='Yes';
str2='No';
default='Yes';
choice = questdlg(qstring,title,str1,str2,default);
                % Handle response
                    switch choice
                        case 'Yes';
                            disp([choice 'Great, let''s keep going then!'])
                            next=1;
                        case 'No';
                            disp([choice 'Okay, please finish selecting devices to throw out'])
                            next=0;
                    end
                    
assignin('caller', 'next', next);

end





function clickTest(line_handle, event)

  good_subplots = evalin('caller', 'good_subplots');
    
  current_tag = get(line_handle, 'Tag');
  
  %split_tag = strsplit(current_tag, ',');
  
  
  subplot_num = str2num(current_tag);
  %subplot_state = str2num(split_tag{1, 2});
  subplot_state = good_subplots(subplot_num);
  
  disp(sprintf('%d, %d', subplot_num, subplot_state));
    
  % Toggle the selection based on the last character in the tag 
  % (0 = off, 1 = on)
  
  p1 = get(line_handle,'Color');


  if (subplot_state)
      
    set(line_handle, 'Color', [1 0 0]);

     good_subplots(subplot_num) = 0;
     
  else
      
    set(line_handle, 'Color',[0 0 1]);

     
     good_subplots(subplot_num) = 1;
     
      
  end
  
%  disp(good_subplots);
  assignin('caller', 'good_subplots', good_subplots);
    
end

function good_wells = getCommonGoodWells(struct_array)

% Get the number of targets
num_structs = length(struct_array);


% Make the first comparison
good_wells = intersect(struct_array{1}.index_dev_to_analyze,...
    struct_array{2}.index_dev_to_analyze);

% If there is more than one good well, make the comparison to each well
if (num_structs > 2)
    
    for i = 3:num_structs
        
        good_wells = intersect(good_wells,...
              struct_array{i}.index_dev_to_analyze);
        
    end
    
    
end


end

function  makeBoxplot(data_structs, normalizer_struct, wells_to_analyze)

%% Get the variables from the data structures
% Get the number of targets
num_targets = length(data_structs);

% Get the number of wells to correlate
num_wells = length(wells_to_analyze);

% Preallocate the AUC matrix
all_auc = zeros(num_wells, num_targets);

% Preallocate the target names cell array
target_names = cell(num_targets, 1);

% Get the normalizer AUC
norm_auc = normalizer_struct.AUC(wells_to_analyze);

% Get the AUC and target names for each dataset
% AUC should be a vertical vector. Slicing with only the wells_to_analyze
for struct_index = 1:num_targets
    
   % Get the AUC vector from the data structure 
   all_auc(:, struct_index) =...
       data_structs{struct_index}.AUC(wells_to_analyze) ./ norm_auc;
   
   % Get the target name from the data structure
   target_names{struct_index, 1} = data_structs{struct_index}.name;
    
   
end


%% Make the boxplot
boxplot(all_auc,'Labels', target_names);


end % function (corrPlot)

function gauss_y = gaussReturn(a,b,c,x)
	gauss_y = a.*exp(-((x-b)./c).^2);
end

