function y = medianfilter( wind,x)
    N = length(x) +1 - wind;
    y = zeros(1,numel(x) - wind+1);
    for j = 1:N
        y(j) = median(x(j:(j+wind-1)));
    end
end

