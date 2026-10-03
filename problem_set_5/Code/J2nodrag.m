function xdot = J2nodrag(x, mu, J2, ro)

    % Pull terms from the state
    r = x(1); theta = x(3); rd = x(4); phid = x(5); thetad = x(6);

    % Compute acceleration terms
    rdd = r*(phid^2)*(cos(theta)^2) + r*(thetad^2) - mu/(r^2) + ...
          (3*mu*J2*(ro^2))/(2*(r^4))*(3*(sin(theta)^2)-1);
    phidd = -((2*rd*phid)/r) + 2*phid*thetad*tan(theta);
    thetadd = -((2*rd*thetad)/r) - ...
              cos(theta)*sin(theta)*((phid^2)+((3*mu*J2*(ro^2))/(r^5)));

    % Assemble state derivative as a column vector
    xdot = [x(4:6); rdd; phidd; thetadd];
end