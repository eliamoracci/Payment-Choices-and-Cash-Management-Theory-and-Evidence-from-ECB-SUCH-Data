function c=cost(m,bar_m,r,gamma)
if m<bar_m
    c=r*m;
else
    c=r*m+(m-bar_m)^(gamma);
end
end
    