% Using Kepler's equation to solve eccentric anomaly given mean anomaly
% Pass:
%     e = eccentricity
%     M = mean anomaly
%   tol = tolerance of stop condition
% Return:
%     E = eccenctric anomaly
function E = kepler_E(e, M)
    
    % Determine an initial guess for E
    if (M < -pi && 0 ) || M > pi
        E1 = M - e;
    else
        E1 = M + e;
    end

    % Begin interation
    while 1
        % Set value of this step to previous step and iterate
        E = E1;
        E1 = E - ((E-(e*sin(E))-M)/(1-(e*cos(E))));

        % Iterate until the error between the iterations is less than tol
        if abs(E - E1) < (10^(-12))
            break
        end
    end
end