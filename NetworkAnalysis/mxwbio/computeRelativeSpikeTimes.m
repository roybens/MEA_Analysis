function [ values ] = computeRelativeSpikeTimes( fileManagerObj, varargin )
    % COMPUTERELATIVESPIKETIMES computes the spike times on the channels of
    % one recording file in 'fileManagerObj' with respect to the first 
    % frame number of that recording file. In case the first frame number
    % was not stored by the recording software, the spike times are
    % computed with respect to the first spike.
    % 
    % [values] = mxw.util.computeRelativeSpikeTimes(fileManagerObj);
    %
    %   -The input parameters for this function are:
    %    -fileManagerObj: object of the class 'mxw.fileManager'
    %    -varargin: ...
    %    -'file': in case the object from the class 'mxw.fileManager'
    %             contains more than one recording file, the file number to
    %             use has to be declared  
    %
    %   -The output parameter for this method is:
    %    -values: struct containing the 'time' and 'channel' where the 
    %             spike occured.
    %                                                
    %  -Examples
    %     -Considering we want to have the relative spike times of file 
    %     number 12 in the object 'networkAnalysisFiles':
    %
    %     [values] = ...
    %       mxw.util.computeRelativeSpikeTimes(networkAnalysisFiles, ...
    %       'file', 12);
    %
    %     -Considering 'networkAnalysisFiles' has only one file, and we 
    %     want its relative spike times:
    %
    %     [values] = ...
    %       mxw.networkActivity.computeNetworkAct(networkAnalysisFiles);
    %
    %

p = inputParser;

p.addParameter('file', [])

p.parse(varargin{:});
args = p.Results;

if fileManagerObj.nFiles > 1
    if isempty(args.file)
        error('Please specify the file to use, or use a fileManager object containing only one file')
    else
        spikeTimes = double(fileManagerObj.rawMap(args.file).spikes.frameno) / double(fileManagerObj.fileObj(args.file).samplingFreq);
        
        if ischar(fileManagerObj.fileObj(args.file).firstFrameNum)
            firstFrameTime = min(spikeTimes);
        else
            firstFrameTime = double(fileManagerObj.fileObj(args.file).firstFrameNum) / double(fileManagerObj.fileObj(args.file).samplingFreq);
        end
        
        relativeSpikeTimes = spikeTimes - firstFrameTime;
        channel = fileManagerObj.rawMap(args.file).spikes.channel;
    end
else
    spikeTimes = double(fileManagerObj.rawMap.spikes.frameno) / double(fileManagerObj.fileObj.samplingFreq);
    
    if ischar(fileManagerObj.fileObj.firstFrameNum)
        firstFrameTime = min(spikeTimes);
    else
        firstFrameTime = double(fileManagerObj.fileObj.firstFrameNum) / double(fileManagerObj.fileObj.samplingFreq);
    end
    
    relativeSpikeTimes = spikeTimes - firstFrameTime;
    channel = fileManagerObj.rawMap.spikes.channel(relativeSpikeTimes > 0);
    relativeSpikeTimes = relativeSpikeTimes(relativeSpikeTimes > 0);
end

values.time = single(relativeSpikeTimes);
values.channel = uint16(channel);
end