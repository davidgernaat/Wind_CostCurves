%% Prep
% there are a few cells in L-Amerika with extreme high wind speeds. The neighbouring cells have a wind speed of 8
% the high wind speeds are reduced to 8 also.

if vcorrection==1
    for m=1:12
        for r=1:nr
            for c=1:nc
                if V{m}(r,c) > 10
                    V{m}(r,c) = 8;
                end
            end
        end
    end
end

if areafile==1
    Area = AreaHoogwijk;
end

% figure(1);clf;imagesc(V{13});axis image