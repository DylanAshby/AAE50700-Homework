function xdot = J2drag(x, mu, J2, ro, CD, A, rho0, H)

    % Pull terms from the state
    r = x(1); theta = x(3); rd = x(4); phid = x(5); thetad = x(6);

    % Compute the velocity magnitude
    v = sqrt((rd^2)+((r^2)*(cos(theta)^2)*(phid^2))+((r^2)*(thetad^2)));
    
    % Compute atmospheric density
    rho = rho0*exp(-(r - 6978)/H);

    % Compute acceleration terms
    rdd = -((CD*A*rho*v*rd)/(2)) + r*(phid^2)*(cos(theta)^2) +...
          r*(thetad^2) - mu/(r^2) +...
          (3*mu*J2*(ro^2))/(2*(r^4))*(3*(sin(theta)^2)-1);
    phidd = -((CD*A*rho*v*phid)/(2)) -((2*rd*phid)/r) +...\
            2*phid*thetad*tan(theta);
    thetadd = -((CD*A*rho*v*thetad)/(2)) - ((2*rd*thetad)/r) - ...
              cos(theta)*sin(theta)*((phid^2)+((3*mu*J2*(ro^2))/(r^5)));

    % Assemble state derivative as a column vector
    xdot = [x(4:6); rdd; phidd; thetadd];
end