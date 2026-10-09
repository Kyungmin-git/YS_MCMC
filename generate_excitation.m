function [candidatet0, candidateqn_n, gastype, logJacobian, logproposal] = ...
    generate_excitation(currentt0, currentqn_n)
%GENERATELPEXC13_10_N4 Perturb excitation times or normalized masses.
%   Inputs are nonempty row vectors of equal length:
%     currentt0   Excitation times within [0, tau] (s).
%     currentqn_n Nonnegative excitation mass fractions summing to one.
%   Set global tau to the observation-window duration before calling.
%
%   gastype = 1: perturb times, then sort time/mass pairs by time.
%   gastype = 2: perturb masses, reflect at zero, then renormalize.
%   The excitation count remains fixed; birth/death moves are not used.
%   Outputs and active proposal calculations retain the original interface.

global tau gastype

% Select either active proposal with equal probability.
gastype = randsample(2, 1);
step_fraction = 0.01*abs(normrnd(0, 1));
logJacobian = 0;
logproposal = 0;

n_exc = size(currentt0, 2);
n_perturb = max(1, round(step_fraction*n_exc));

if gastype == 1
    %% Excitation-time perturbation
    indices = randperm(n_exc, n_perturb);
    candidatet0 = currentt0;
    candidatet0(indices) = candidatet0(indices) + ...
        normrnd(0, 0.013, [1, n_perturb]);

    % Preserve the original reflection at the window boundaries.
    candidatet0 = abs(candidatet0);
    outside_window = candidatet0 > tau;
    candidatet0(outside_window) = 2*tau-candidatet0(outside_window);

    % Sort paired times and masses, preserving their association.
    candidateqn_n = currentqn_n;
    gas_pairs = sortrows([candidatet0; candidateqn_n]');
    candidatet0 = gas_pairs(:, 1)';
    candidateqn_n = gas_pairs(:, 2)';

else
    %% Normalized excitation-mass perturbation
    indices = randperm(n_exc, n_perturb);
    candidateqn_n = currentqn_n;
    for j = 1:n_perturb
        candidateqn_n(indices(j)) = currentqn_n(indices(j)) + ...
            normrnd(0, 0.00002, 1);
    end

    candidateqn_n = abs(candidateqn_n);
    candidateqn_n = candidateqn_n/sum(candidateqn_n);
    candidatet0 = currentt0;

    % Retain the supplied Jacobian calculation without changing the kernel.
    mass_difference = sum(candidateqn_n(indices)-currentqn_n(indices));
    n_unperturbed = length(currentqn_n)-n_perturb;
    Jacobian = (1-mass_difference)^n_unperturbed;
    logJacobian = log(Jacobian);
end
end
