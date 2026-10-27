function xd = rwEoM(x, t, u, bt, bdt, bddt, IB, IA0, Jw, mA, d)
    
    % Break the state vector into its components
    th = x(1:3); thd = x(7:9); ad = x(10:12);

    % Find the value of the beta angles for this time step
    b = bt(t); bd = bdt(t); bdd = bddt(t);

    % Define 2nd axis body frame vector and its skew-symmetric matrix
    bh2 = [0; 1; 0];
    bh2x = [0 0 1; 0 0 0; -1 0 0];

    % Define the 2nd axis rotation about b2 body frame vector
    R2 = [cos(b)  0  sin(b);
          0          1  0;
         -sin(b)  0  cos(b)];
    
    % Find IA matrix, the IAC matrix, and its derivative
    IA = R2*IA0*R2';
    IAC = 2*IA + 2*mA*(d^2)*(eye(3) - bh2*bh2');
    IACd = 2*bd * (bh2x*IA - IA*bh2x);

    % Find the Euler sequence matrix and its derivative
    E = [sin(th(3))*sin(th(2)) cos(th(3))  0;
         sin(th(2))*cos(th(3)) -sin(th(3)) 0;
         cos(th(2))            0           1];
    Ed = [cos(th(3))*sin(th(2))*thd(3) + sin(th(3))*cos(th(2))*thd(2) ...
          -sin(th(3))*thd(3) 0;
          cos(th(2))*cos(th(3))*thd(2) - sin(th(2))*sin(th(3))*thd(3) ...
          -cos(th(3))*thd(3) 0;
          -sin(th(2))*thd(2) 0 0];

    % Find the partials of E with respect to theta
    dE1 = zeros(3);
    dE2 = [sin(th(3))*cos(th(2)) 0 0;
           cos(th(2))*cos(th(3)) 0 0;
          -sin(th(2))            0 0];
    dE3 = [cos(th(3))*sin(th(2)) -sin(th(3)) 0;
          -sin(th(2))*sin(th(3)) -cos(th(3)) 0;
           0                     0           0];
    
    % Assemble partials times thetadot column by column
    dEthd = [dE1*thd dE2*thd dE3*thd];

    % Compute IS and HS matricies
    IS = IB + IAC + Jw;
    HS = IB*E*thd + IAC*(E*thd + bd*bh2) + Jw*(E*thd + ad);

    % Compute the f vector
    f = [-E' * (IS*Ed*thd + IACd*(E*thd + bd*bh2) + IAC*bdd*bh2) + ...
         (dEthd - Ed)'*HS; u - Jw*Ed*thd];

    % Compute the M matrix
    M = [E'*IS*E E'*Jw; Jw*E Jw];

    % Find the angle accelerations and pass through the angular rates
    xd = [thd; ad; M\f];

end