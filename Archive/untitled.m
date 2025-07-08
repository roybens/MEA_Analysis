% Define a sample signal with two Gaussian peaks and some noise
x = linspace(0, 20, 1000);
y = exp(-0.2*(x-5).^2) + 0.8*exp(-0.1*(x-12).^2) + 0.05*randn(size(x));

% Find peaks and compute their properties
[peaks, locs, widths, prominences, left_bases, right_bases] = findpeaks(y, x);

% Calculate the width at full prominence
full_prominence_widths = right_bases - left_bases;

% Plot the original signal
figure;
plot(x, y, '-b'); hold on;
plot(locs, peaks, 'vr', 'MarkerFaceColor', 'r'); % Mark the peaks

% Mark the left and right bases of each peak
for i = 1:length(peaks)
    % Left base
    plot(left_bases(i), y(x == left_bases(i)), 'go', 'MarkerFaceColor', 'g');
    
    % Right base
    plot(right_bases(i), y(x == right_bases(i)), 'mo', 'MarkerFaceColor', 'm');
    
    % Draw a line connecting the left and right bases to show full prominence width
    line([left_bases(i), right_bases(i)], [y(x == left_bases(i)), y(x == right_bases(i))], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    
    % Annotate the width at full prominence
    mid_point = (left_bases(i) + right_bases(i)) / 2;
    text(mid_point, min(y), sprintf('%.2f', full_prominence_widths(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Color', 'black', 'FontSize', 8);
end

% Add annotations
legend('Original Signal', 'Peaks', 'Left Base', 'Right Base', 'Full Prominence Width');
title('Peak Analysis: Width at Full Prominence for Each Peak');
xlabel('Time');
ylabel('Amplitude');
grid on;
hold off;

% Display the width at full prominence for each peak
disp('Full Prominence Width for Each Peak:');
disp(full_prominence_widths);