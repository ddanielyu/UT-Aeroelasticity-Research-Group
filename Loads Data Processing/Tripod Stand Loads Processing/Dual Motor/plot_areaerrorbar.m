function plot_areaerrorbar(x_axis,y_mean,y_err,color)

    % Default options        
    if color == [0,0.447,0.741]
        options.color_area = [128 193 219]./255;    % Blue theme
    elseif color == [0.850,0.325,0.098]
        options.color_area = [0.92,0.53,0.37];    % Orange theme
    elseif color == [0.929000000000000,0.694000000000000,0.125000000000000]; %yellow theme
        options.color_area = [0.96,0.76,0.30];
    end
    options.alpha      = 0.5;
    options.line_width = 2;
    
    x_vector = [x_axis,fliplr(x_axis)];
    patch = fill(x_vector, [y_mean+y_err,fliplr(y_mean-y_err)], options.color_area);
%     hold on
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', options.alpha);
    mean_plot = plot(x_axis, y_mean, '-','color', color, ...
        'LineWidth', options.line_width);
    set(mean_plot,'HandleVisibility','off');
%     hold off
    
end