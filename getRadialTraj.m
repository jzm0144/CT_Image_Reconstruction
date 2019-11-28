function ktraj = getRadialTraj(readout_pts, no_of_projections)

step_size = 1/readout_pts;
Kx = (-.5:step_size:.5-step_size)';
Ky = zeros(size(Kx));
Ko = Kx + 1i*Ky;

ktraj = zeros(readout_pts,no_of_projections);
th = linspace(0, pi, no_of_projections+1);

for nn = 1:no_of_projections
    ktraj(:,nn) = Ko*exp(1i*th(nn));
end
end