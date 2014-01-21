
    
    [nv vx vy] = size(v);
    [ns sx sy] = size(s);
    figure(1)
    for c = 1:nv
        hold on
        plot(v{c,:,:},v{c,:,:},'.');
        hold off
    end

figure(2)
    for c = 1:nv
        hold on
        plot(s{c,:,:},s{c,:,:});
        hold off
    end
