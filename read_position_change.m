% READ_POSITION_CHANGE
%
% This function reads one of Simon Hunt's position change files
% and returns the data.
% 
% Usage: [image1, image2, time, box1, box2, box3, metadata] ...
%                                      = read_position_change(filename)
%
%     filename - full path or relative path to position change file
%     image1 - integer number of first (reference) image 
%     image2 - integer number of second image
%     time   - timestamp of second image
%     box1   - movement of first box
%     box2   - movement of second box
%     box3   - movemebt of third box
%     metadata - structure with header info
% 
% Andrew Walker - 6/3/2014

function [image1, image2, time, box1, box2, box3, metadata] ...
                                      = read_position_change(filename)

    fid = fopen(filename);
    
    % Read metadata lines
    in_file_line = fgetl(fid); 
    ab_by_line = fgetl(fid);
    opts_line = fgetl(fid);
    fgetl(fid); % Blank line
    fgetl(fid); % Table header line
    posline1 = fgetl(fid);
    posline2 = fgetl(fid);
    posline3 = fgetl(fid);
    posline4 = fgetl(fid);
    fgetl(fid); % Blank line
    
    % Fill in the metadata - FIXME: needs finishing
    tok = regexp(posline1, 'top,\s+(\d+),\s+(\d+),\s+(\d+)', 'tokens');
    tops = [str2double(tok{1}{1}) str2double(tok{1}{2}) str2double(tok{1}{3})];
    tok = regexp(posline2, 'bottom,\s+(\d+),\s+(\d+),\s+(\d+)', 'tokens');
    bots = [str2double(tok{1}{1}) str2double(tok{1}{2}) str2double(tok{1}{3})];
    tok = regexp(posline3, 'left,\s+(\d+),\s+(\d+),\s+(\d+)', 'tokens');
    lefts = [str2double(tok{1}{1}) str2double(tok{1}{2}) str2double(tok{1}{3})];
    tok = regexp(posline4, 'right,\s+(\d+),\s+(\d+),\s+(\d+)', 'tokens');
    rights = [str2double(tok{1}{1}) str2double(tok{1}{2}) str2double(tok{1}{3})];
    
    % \d*p?_?
    tok = regexp(filename, '(\w+\d+)_(\d+)tons_(\d+)C_(\d+)s_', 'tokens');
    if (~isempty(tok))
        expt_name = tok{1}{1};
        nominal_load = str2double(tok{1}{2});
        nominal_temperature = str2double(tok{1}{3});
        nominal_period = str2double(tok{1}{4});
        nominal_strain = 0.0;
    else
        tok = regexp(filename, '(\w+\d+)_(\d+)tons_(\d+)C_(\d+)p_(\d+)s_', 'tokens');
        if (~isempty(tok))
            expt_name = tok{1}{1};
            nominal_load = str2double(tok{1}{2});
            nominal_temperature = str2double(tok{1}{3});
            nominal_period = str2double(tok{1}{5});
            nominal_strain = str2double(tok{1}{4});
        else
            % Cannot parse file name!
            expt_name = filename;
            nominal_load = NaN;
            nominal_temperature = NaN;
            nominal_period = 30;
            nominal_strain = NaN;
        end
    end
    metadata = struct('InputFile', in_file_line, ...
                      'AnalysisBy', ab_by_line, ...
                      'BoxTops', tops, 'BoxBotoms', bots, ...
                      'BoxLefts', lefts, 'BoxRights', rights, ...
                      'ExperimentName', expt_name, ...
                      'NominalLoad', nominal_load, ...
                      'NominalTemp', nominal_temperature, ...
                      'NominalPeriod', nominal_period, ...
                      'NominalStrain', nominal_strain);
   
    % Read all the data and 
    data = fscanf(fid, '%f, %f, %f, %f, %f, %f', [6 inf]);
    fclose(fid);
    
    % put data into output arrays
    image1 = round(data(1,:));
    image2 = round(data(2,:));
    time = data(3,:);
    box1 = data(4,:);
    box2 = data(5,:);
    box3 = data(6,:);
    
end