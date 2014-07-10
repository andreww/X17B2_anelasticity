% ReadPositionChange.
%
% Read position change files
%
% syntx: [position_change, time_stamps, box_positions, number_boxes, number_images] = ReadPositionChangeFile('file_name')
%
% Simon Hunt 2014

function varargout = ReadPositionChangeFile(file_name)


%reads in the csv output from the image_analysis.m script
data = importdata(file_name,',');

%extracts the time stamps of the images from the input file.
time_stamps = data.textdata(:,3);
time_stamps = time_stamps(10:length(time_stamps));
time_stamps = char(time_stamps);
time_stamps = str2num(time_stamps); %#ok<ST2NM>

%moves the time stamps so that t=0 for the first image
% time_stamps = time_stamps - min(time_stamps);


%sorts the numeric data into the box positions in the image and their displacements with time.
all_numbers = data.data;
size_all_numbers = size(all_numbers);

foil_change_data = all_numbers(5:size_all_numbers(1), 1:size_all_numbers(2));
box_positions = all_numbers(1:4, 1:size_all_numbers(2));


%define variables describing amount of data
number_images = length(time_stamps);

number_boxes = size(foil_change_data,2);


varargout = {foil_change_data, time_stamps, box_positions, number_boxes, number_images};