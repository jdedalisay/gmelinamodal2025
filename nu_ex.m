function [c, ceq] = nu_ex(x)
    % Define the nonlinear inequality constraints
    ey = 10^(x(1)+11);
    ex = 10^(x(2)+11);
    ez = 343713300;
    nxy = x(3);
    nyz = 0.246;
    nxz = 0.34;
    c(1) = -1*(1 -nxy^2*ey/ex -nyz^2*ez/ex -nxz^2*ez/ex -2*nxy*nyz*nxz*ez/ex);
    % No equality constraints
    ceq = [];
end