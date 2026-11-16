function xdot = threeBodyPhaseSpace(t, x, mu, muM, dM, nM, phiM)

    % Unpack extended state x = [q; p; q0; p0]
    q  = x(1:3);
    p  = x(4:6);
    q0 = x(7);

    % Moon position at time coordinate q0
    R = dM*[cos(nM*q0 + phiM); sin(nM*q0 + phiM); 0];

    % Equations of motion
    qdot  = p;
    pdot  = -mu*q/norm(q)^3 - muM*((q - R)/norm(q - R)^3 + R/norm(R)^3);
    q0dot = 1;
    p0dot = muM*nM*(1/norm(R - q)^3 - 1/dM^3)*dot(q, cross([0;0;1], R));

    % Pack derivative
    xdot = [qdot; pdot; q0dot; p0dot];
end