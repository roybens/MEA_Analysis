function plotNetworkActivityModified( networkActivityStruct, varargin )
    % NETWORKACTIVITY plots the network activity previously computed in 
    % 'mxw.networkActivity.computeNetworkAct'.
    % 
    % mxw.plot.networkActivity(networkActivityStruct);
    %
    %   -The input parameters for this function are:
    %    -networkActivityStruct: structure with the 'time' and 'firingRate'
    %                            values to be plotted. These values are
    %                            previously computed in
    %                            'mxw.networkActivity.computeNetworkAct'
    %    -varargin: ...
    %    -'Figure': boolean that defines whether or not a new figure is 
    %               created for this plot. This is useful when the user 
    %               wants to show several plots in one figure by using the
    %               subplot option
    %    -'Threshold': threshold to be plotted. This threshold should be 
    %                  the same as the threshold output given by 
    %                  'mxw.networkActivity.computeNetworkStats'
    %    -'Axis': boolean that defines if the axis are plotted
    %    -'Xlimits': defines the minimum and maximum values on the x-axis
    %    -'Ylimits': defines the minimum and maximum values on the y-axis
    %    -'Title': title of the plot
    %    -'Ylabel': y-axis label
    %    -'Xlabel': x-axis label
    %                                   
    %  -Examples
    %     -Considering we want to plot the network activity contained in 
    %     the structure 'networkAct' and the threshold given by 'thr':
    %
    %     mxw.plot.networkActivity(networkAct, 'Threshold', thr);
    %
    %
    
p = inputParser;

p.addRequired('struct', @(x) isstruct(x));
p.addParameter('Figure', true);
p.addParameter('ThresholdFunction', 'rms');
p.addParameter('Threshold',[])
p.addParameter('Axis', true);
p.addParameter('Xlimits', []);
p.addParameter('Ylimits', []);
p.addParameter('Title', 'Network Activity');
p.addParameter('Ylabel', 'Firing Rate [Hz]');
p.addParameter('Xlabel', 'Time [s]');
p.addParameter('Color', '#135ba3');

p.parse(networkActivityStruct, varargin{:});
args = p.Results;

if args.Figure
    figure;
end

plot(networkActivityStruct.time, networkActivityStruct.firingRate, 'Color', args.Color)
hold



switch args.ThresholdFunction

    case 'Fixed'
    plot(args.Threshold*ones(ceil(networkActivityStruct.time(end)),1))
    case 'RMS'
    plot(args.Threshold*ones(ceil(networkActivityStruct.time(end)),1))

    case 'Adaptive'
    plot(args.Threshold)
end




if ~(args.Axis)
    axis off;
end

if ~isempty(args.Xlimits)
    xlim(args.Xlimits)
end

if ~isempty(args.Ylimits)
    ylim(args.Ylimits)
end

title(args.Title);
xlabel(args.Xlabel);
ylabel(args.Ylabel);
end