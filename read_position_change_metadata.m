% READ_POSITION_CHANGE_METADATA
%
% This function reads one of Simon Hunt's position change files
% and returns the data.
% 
% Usage: [metadata] = read_position_change_metadata(filename)
%
%     filename - full path or relative path to position change file

%     metadata - structure with header info
% 
% Andrew Walker - 6/3/2014

function [metadata] = read_position_change_metadata(filename)
    

    
    tok = regexp(filename, '(\w+_?\d+)_(\d+)tons?_(\d+)C_(\d+)s_', 'tokens');
    if (~isempty(tok))
        expt_name = tok{1}{1};
        nominal_load = str2double(tok{1}{2});
        nominal_temperature = str2double(tok{1}{3});
        nominal_period = str2double(tok{1}{4});
        nominal_strain = 0.0;
    else
        tok = regexp(filename, '(\w+_?\d+)_(\d+)tons?_(\d+)C_(\d+)p_(\d+)s_', 'tokens');
        if (~isempty(tok))
            expt_name = tok{1}{1};
            nominal_load = str2double(tok{1}{2});
            nominal_temperature = str2double(tok{1}{3});
            nominal_period = str2double(tok{1}{5});
            nominal_strain = str2double(tok{1}{4});
        else
            tok = regexp(filename, '(\w+_?\d+)_(\d+)tons?_ramping_(\d+)s_', 'tokens');
            if (~isempty(tok))
                expt_name = tok{1}{1};
                nominal_load = str2double(tok{1}{2});
                nominal_temperature = NaN;
                nominal_period = str2double(tok{1}{3});
                nominal_strain = 0.0;
            else
                % Cannot parse file name!
                expt_name = filename;
                nominal_load = NaN;
                nominal_temperature = NaN;
                nominal_period = NaN;
                nominal_strain = NaN;
            end
        end
    end
    metadata = struct('ExperimentName', expt_name, ...
                      'NominalLoad', nominal_load, ...
                      'NominalTemp', nominal_temperature, ...
                      'NominalPeriod', nominal_period, ...
                      'NominalStrain', nominal_strain);
end