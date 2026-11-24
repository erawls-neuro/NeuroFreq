function y = stat_restore(a, shape, order)
shape(order(1)) = size(a,1);
y = ipermute(reshape(a,shape(order)),order);
end
