function pos_data = read_pos(pos_file)

% initialize data structure
pos_data = struct('year',[],'month',[],'day',[],'hour',[],'min',[],'sec',[],'x',[],'y',[],'z',[]);

file = fopen(pos_file);
% skip header lines
for i=1:22
    fgetl(file);
end

while ~feof(file)
    % collect time data
    timedata = fgetl(file);
    % check for last line
    if length(timedata) == 3
        break
    end
    
    pos_data.year = [pos_data.year str2num(timedata(4:7))];
    pos_data.month = [pos_data.month str2num(timedata(10))];
    pos_data.day = [pos_data.day str2num(timedata(12:13))];
    pos_data.hour = [pos_data.hour str2num(timedata(15:16))];
    pos_data.min = [pos_data.min str2num(timedata(18:19))]; %#ok<*ST2NM>
    pos_data.sec = [pos_data.sec str2double(timedata(21:30))];
    
    % collect position data
    posdata = fgetl(file);
    pos_data.x = [pos_data.x str2double(posdata(7:18))];
    pos_data.y = [pos_data.y str2double(posdata(21:32))];
    pos_data.z = [pos_data.z str2double(posdata(35:46))];
end

end