% Equations of motion for the lunar third-body problem without extending
% the phase space (time appears explicitly through the Moon's position)
% Pass:
%      t = time
%      x = state vector [q; p]
%     mu = Earth gravitational parameter
%    muM = Moon gravitational parameter
%     dM = Moon orbital radius
%     nM = Moon mean motion
%   phiM = Moon initial phase
% Return
%   xdot = time derivative of the state vector
function xdot = threeBodyNonExtended(t, x, mu, muM, dM, nM, phiM)
    % Unpack state x = [q; p]
    q = x(1:3);
    p = x(4:6);

    % Moon position at time t
    R = dM*[cos(nM*t + phiM); sin(nM*t + phiM); 0];

    % Equations of motion
    qdot = p;
    pdot = -mu*q/norm(q)^3 - muM*((q - R)/norm(q - R)^3 + R/norm(R)^3);

    % Pack derivative
    xdot = [qdot; pdot];
end