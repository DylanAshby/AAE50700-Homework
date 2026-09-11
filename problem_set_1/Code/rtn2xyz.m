% Convert cartesian rtn frame coordinates to xyz frame
% Pass:
%      xrtn = rtn frame Cartesian coordinates
%      O = Right Ascension of the Ascending Node
%      i = Inclination
%      u = Argument of Latitude
% Return
%      xxyz = xyz frame Cartesian coordinates
function xxyz = rtn2xyz(xrtn, O, i, u)
    % Define half rotation coefficients for the rotation quaternion
    O = 0.5*O;          % Half RAAN
    i = 0.5*i;          % Half inclination
    u = 0.5*u;          % Half Argument of Latitude

    % Define the RTN to XYZ rotation quaternion
    qo = (cos(O)*cos(i)*cos(u))-(sin(O)*cos(i)*sin(u));
    q1 = (cos(O)*sin(i)*cos(u))+(sin(O)*sin(i)*sin(u));
    q2 = (-cos(O)*sin(i)*sin(u))+(sin(O)*sin(i)*cos(u));
    q3 = (cos(O)*cos(i)*sin(u))+(sin(O)*cos(i)*cos(u));

    % Define the quaternion rotation operator
    Q = [((2*(qo^2))-1+(2*(q1^2))),...
        ((2*q1*q2)+(2*qo*q3)),...
        ((2*q1*q3)-(2*qo*q2));...
        ((2*q1*q2)-(2*qo*q3)),...
        ((2*(qo^2))-1+(2*(q2^2))),... 
        ((2*q2*q3)+(2*qo*q1));...
        ((2*q1*q3)+(2*qo*q2)),...
        ((2*q2*q3)-(2*qo*q1)),...
        ((2*(qo^2))-1+(2*(q3^2)))]';

    % Convert xrtn from RTN frame to XYZ frame using the quarternion
    % rotation operator
    xxyz = Q*xrtn;
end