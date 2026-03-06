function j = no_terms(r)

j = 0;
for i = 1:length(r)
    if r(i) > 5*10^(-16)
        j = j+1;
    end
end