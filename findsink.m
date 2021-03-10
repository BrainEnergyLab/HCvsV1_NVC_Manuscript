function rateconsump = findsink(region,state,Vmax,Km)

%% Function to make a conc dependent consumption (i.e. sink)

nr = length(region.y); % Number of columns
rateconsump = zeros(1,nr);

try
    ncProp = state.u;
    %formula for rate of consumption 
    rateconsump = -Vmax.*ncProp./(Km+ncProp);
    
    if isnan(rateconsump)
        rateconsump(1,nr) = 0;
    end
    
catch
    rateconsump(1,nr) = 0;
end

end %end of func