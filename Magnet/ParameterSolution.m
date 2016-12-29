nxp1=nx-1;
nyp1=ny-1;

% TEz components
eps_r_x     = ones (nx  , nyp1, nz);
eps_r_y     = ones (nxp1, ny,   nz);
mu_r_z      = ones (nx  , ny,   nz);
k_r_z       = ones (nx  , ny,   nz);
sigma_e_x   = zeros(nx  , nyp1, nz);
sigma_e_y   = zeros(nxp1, ny,   nz);
sigma_m_z   = zeros(nx  , ny,   nz);

% TMz components
eps_r_z     = ones (nxp1, nyp1, nz);
mu_r_x      = ones (nxp1, ny,   nz);
mu_r_y      = ones (nx  , nyp1, nz);
k_r_x       = ones (nxp1, ny,   nz);
k_r_y       = ones (nx  , nxp1, nz);
sigma_e_z   = zeros(nxp1, nyp1, nz);
sigma_m_x   = zeros(nxp1, ny,   nz);
sigma_m_y   = zeros(nx  , nyp1, nz);

eps_r_x(1:nx,2:ny, :) = 0.5 * (epsilon(1:nx,2:ny, :) + epsilon(1:nx,1:ny-1, :));
eps_r_y(2:nx,1:ny, :) = 0.5 * (epsilon(2:nx,1:ny, :) + epsilon(1:nx-1,1:ny, :));
eps_r_z(2:nx,2:ny, :) = 0.25 * (epsilon(2:nx,2:ny, :) + epsilon(1:nx-1,2:ny, :) + epsilon(2:nx,1:ny-1, :) + epsilon(1:nx-1,1:ny-1, :));

sigma_e_x(1:nx,2:ny,:) = 0.5 * (sigmax(1:nx,2:ny) + sigmax(1:nx,1:ny-1, :));
sigma_e_y(2:nx,1:ny,:) = 0.5 * (sigmay(2:nx,1:ny) + sigmay(1:nx-1,1:ny, :));
sigma_e_z(2:nx,2:ny,:) = 0.25 * (sigmaz(2:nx,2:ny, :) + sigmaz(1:nx-1,2:ny, :) + sigmaz(2:nx,1:ny-1, :) + sigmaz(1:nx-1,1:ny-1, :));

mu_r_x(2:nx,1:ny,:) = 2 * (mu(2:nx,1:ny, :) .* mu(1:nx-1,1:ny, :)) ./(mu(2:nx,1:ny, :) + mu(1:nx-1,1:ny, :));
mu_r_y(1:nx,2:ny,:) = 2 * (mu(1:nx,2:ny, :) .* mu(1:nx,1:ny-1, :)) ./(mu(1:nx,2:ny, :) + mu(1:nx,1:ny-1, :));
mu_r_z(1:nx,1:ny,:) = mu(1:nx,1:ny, :);

k_r_x(2:nx,1:ny,:) = 2 * (mu(2:nx,1:ny, :) .* mu(1:nx-1,1:ny, :)) ./(mu(2:nx,1:ny, :) + mu(1:nx-1,1:ny, :));
k_r_y(1:nx,2:ny,:) = 2 * (mu(1:nx,2:ny, :) .* mu(1:nx,1:ny-1, :)) ./(mu(1:nx,2:ny, :) + mu(1:nx,1:ny-1, :));
k_r_z(1:nx,1:ny,:) = mu(1:nx,1:ny, :);

sigma_m_x(2:nx,1:ny,:) = 2 * (sigma_starx(2:nx,1:ny, :) .* sigma_starx(1:nx-1,1:ny, :)) ./(sigma_starx(2:nx,1:ny, :) + sigma_starx(1:nx-1,1:ny, :));
sigma_m_y(1:nx,2:ny,:) = 2 * (sigma_stary(1:nx,2:ny, :) .* sigma_stary(1:nx,1:ny-1, :)) ./(sigma_stary(1:nx,2:ny, :) + sigma_stary(1:nx,1:ny-1, :));
sigma_m_z(1:nx,1:ny,:) = sigma_starz(1:nx,1:ny, :);