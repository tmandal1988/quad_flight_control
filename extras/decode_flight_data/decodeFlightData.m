function dataStruct = decodeFlightData(varargin)
%DECODEFLIGHTDATA decodes flight data

%Assumes data_file.dat is in the path and only one file is with that name

p = inputParser;
addParameter(p,'data_format', 'ff9f3f3f11f16f3f6f7f4f88f', @ischar);
addParameter(p,'top_count', 100000, @isnumeric);

parse(p, varargin{:});

fullSimDir = getProjectDir('fullSim');
fid = fopen(fullfile(fullSimDir, 'data_file.dat'),'r');
idx = 1;
top_count = p.Results.top_count;
counter = 1;

data_format = p.Results.data_format;

% Loop through the data file and convert the data from binary to specified
% data format
while(1)
    counter = counter + 1;
    D = struct_read(fid, data_format);
    if( ~strcmp(D{end}, 'eof') )
        time(idx)               = D{1}*1e-6;
        count(idx)              = D{2};
        imuRaw(idx, :)          = D{3}(1:9);
        posRaw(idx, :)          = D{4}(1:3);
        velRaw(idx, :)          = D{5}(1:3);
        baro(idx, :)            = D{6}(1:11);  
        ekfStates(idx, :)       = D{7}(1:16);
        mahFilter(idx, :)       = D{8}(1:3);
        gpsValid(idx, :)        = D{9}(1:6);
        rcCmds(idx, :)          = D{10}(1:7);       
        actCmds(idx, :)         = D{11}(1:4);        
        fcsDebug(idx, :)        = D{12}(1:88);
        idx = idx + 1;
    else
        break;
    end

    if counter > top_count
        break;
    end
end

dataStruct.time = time;
dataStruct.count = count;
dataStruct.imuRaw = imuRaw;
dataStruct.posRaw = posRaw;
dataStruct.velRaw = velRaw;
dataStruct.baro = baro;
dataStruct.ekfStates = ekfStates;
dataStruct.mahFilter = mahFilter;
dataStruct.gpsValid = gpsValid;
dataStruct.rcCmds = rcCmds;
dataStruct.actCmds = actCmds;
dataStruct.fcsDebug = fcsDebug;
end

