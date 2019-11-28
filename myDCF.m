function [dcf] = myDCF(x, y)
%dcf = zeros(1, x*y);

v = ones(1, x);
for i = 1:16
    v(1, i) =  i/16;
    v(1, 512-i) = 1 - (16-i)/16;

end
dcf = v
end
