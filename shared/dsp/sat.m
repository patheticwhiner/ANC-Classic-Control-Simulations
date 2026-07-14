function x = sat(x, a)
    a = abs(a);
    if(x>a)
        x = a;
    elseif(x<-a)
        x = -a;
    end
end