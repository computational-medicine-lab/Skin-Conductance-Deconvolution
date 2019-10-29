function [out] = pickhighest(in, num)
% pickhighest -- Picks the highest num elements of the vector in (abs)
%
% [out] = pickhighest(in, num)
%
% JFM   10/6/2000
% Rev:  10/6/2000

cols = length(in);
out = zeros(1, cols);

if(num >= cols)
    out = in;
    return;
end

inabs = abs(in);
minin = min(inabs);

for i = 1:num
    [x, index] = max(inabs);
    inabs(index) = minin;
    out(index) = in(index);
end

return;