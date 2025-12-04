function y = FIR(filter, X)
    y = sum(filter.*X(1:length(filter)));
end