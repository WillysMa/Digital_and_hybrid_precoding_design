function [Nth, Ntv] = array_dimension(Nt)
    % Find the configuration of UPA that minimizes beam squint effect
    n = ceil(sqrt(Nt));
    for i = n+1:-1:2
        if mod(Nt, i) == 0
            Nth = i;
            Ntv = Nt / i;
            break;
        end
    end
end
