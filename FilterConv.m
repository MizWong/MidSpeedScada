function [ filtData ] = FilterConv( data, coef )

dataLong  = [ data(end-(length(coef)-1)/2 + 1 : end), data, data(1:(length(coef)-1)/2) ];

filtData = conv(dataLong, coef, 'valid');

end
