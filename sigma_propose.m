
function sigma_k_log = sigma_propose(current_sigma_k_log)
%SIGMA_PROPOSE Generate a bounded MCMC proposal for log10(sigma).
%
%   sigma_k_log = SIGMA_PROPOSE(current_sigma_k_log)
%
%   Generates a Gaussian random-walk proposal for the base-10
%   logarithm of the noise standard deviation, sigma.
%
%   INPUT
%       current_sigma_k_log  Current value of log10(sigma).
%
%   OUTPUT
%       sigma_k_log          Proposed value of log10(sigma).
%
%   GLOBAL VARIABLE
%       stepsize2            Standard deviation of the Gaussian
%                            proposal in log10 space.
%
%   PARAMETER BOUNDS
%       -2 <= log10(sigma) <= 1
%       Equivalent to 0.01 <= sigma <= 10.
%
%   PROPOSAL SCHEME
%       1. Generate a Gaussian random-walk perturbation.
%       2. Apply reflective boundary conditions.
%       3. Return the bounded candidate value.

    global stepsize2

    % Bounds in log10 space
    lower = -2;
    upper =  1;

    % Validate inputs
    validateattributes(current_sigma_k_log, {'numeric'}, ...
        {'real', 'finite', 'scalar', '>=', lower, '<=', upper});

    validateattributes(stepsize2, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'nonnegative'});

    % Gaussian random-walk proposal
    proposal = current_sigma_k_log + stepsize2 * randn;

    % Reflect proposal at parameter boundaries
    width = upper - lower;
    reflected = mod(proposal - lower, 2 * width);

    if reflected <= width
        sigma_k_log = lower + reflected;
    else
        sigma_k_log = upper - (reflected - width);
    end

end

