% Convert Keplerian orbital elements to Cartesian states
% Pass:
%      x = Keplerian state vector
%     mu = gravitational parameter of the orbited body
% Return
%      r = xyz position
%      v = xyz velocity
function xcart = kep2cart(xkep, mu)
    % Define orbital elements from state vector input
    a=xkep(1); e=xkep(2); i=xkep(3); O=xkep(4); w=xkep(5); M=xkep(6);

    % Calculate Eccentric Anomaly and True Anomaly from Mean Anomaly
    E = kepler_E(e, M);
    nu = wrapTo2Pi(atan2(sqrt(1-(e^2))*sin(E),cos(E)-e));

    % Calculate orbital radius
    rmag = (a*(1-(e^2)))/(1+(e*cos(nu)));

    % Calculate velocity
    vmag = sqrt(((2*mu)/rmag)-(mu/a));

    % Calculate flightpath angle
    h = sqrt(mu*a*(1-(e^2)));       % Angular momentum
    fpaFind = h/(rmag*vmag);        % Accounting for numerical drift
    if fpaFind > 1
        fpaFind = 1;
    elseif fpaFind < -1
        fpaFind = -1;
    end
    gamma = acos(fpaFind);         % FPA
    if nu > pi                      % Resolve quadrant ambiguity
        gamma = -gamma;
    end

    % Define the position and velocity vectors in the rtn frame
    r = [rmag; 0; 0];
    v = [vmag*sin(gamma); vmag*cos(gamma); 0];

    % Define the Argument of Latitude using the Argument of Periapsis and
    % True Anomaly
    u = (w + nu); 

    % Convert position and velocity from RTN frame to XYZ frame
    r = rtn2xyz(r,O,i,u);
    v = rtn2xyz(v,O,i,u);

    % Define Cartesian state vector
    xcart = [r;v];
end