%% Function: Cole-Cole dielectric model

function epsd = cole_model(omega, params)

    epsd = params.eps_inf ...
        + params.eps_delta ./ (1 + (1j * params.tau .* omega) .^ (1 - params.alpha)) ...
        + params.sigma ./ (1j * params.eps0 .* omega);

end
