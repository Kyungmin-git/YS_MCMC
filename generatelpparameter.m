function y = generatelpparameter(x)
%GENERATELPPARAMETER Generate a bounded MCMC parameter proposal.
%
%   y = GENERATELPPARAMETER(x) randomly selects one model parameter
%   and perturbs it using a Gaussian random-walk proposal.
%
%   Reflective boundary conditions ensure that the proposed value
%   remains within its prescribed parameter bounds.
%
%   INPUT
%       x        Current parameter vector (8 elements).
%
%   OUTPUT
%       y        Proposed parameter vector (same shape as x).
%
%   GLOBAL VARIABLE
%       stepsize Proposal standard deviation for each parameter.
%
%   PARAMETER DEFINITIONS
%       1. Temperature (K)
%       2. Gas mass flow rate (kg/s)
%       3. Gas pocket radius (m)
%       4. Cap thickness (m)
%       5. Negative log10 permeability (m^2)
%       6. Negative log10 gas pocket thickness (m)
%       7. Log10 porosity (%)
%       8. Source-receiver distance (m)
%
%   PROPOSAL SCHEME
%       1. Select one parameter uniformly at random.
%       2. Add a Gaussian perturbation.
%       3. Reflect the proposal at the prescribed boundaries.
%       4. Return the candidate parameter vector.
%
%   NOTE
%       Acceptance or rejection is handled by the main MCMC routine.

    global stepsize

    % Parameter bounds (columns: lower, upper)
    bounds = [
        458.15,          463.15;        % Temperature (K)
        0,               80;            % Gas flow rate (kg/s)
        35,              100;           % Gas pocket radius (m)
        10,              20;            % Cap thickness (m)
        log10(2) + 6,    17;            % -log10 permeability
        -log10(5),       2;             % -log10 pocket thickness
        -1,              1 + log10(2);  % log10 porosity (%)
        0,               50             % Source-receiver distance (m)
    ];

    % Validate inputs
    validateattributes(x, {'numeric'}, ...
        {'real', 'finite', 'vector', 'numel', 8});

    validateattributes(stepsize, {'numeric'}, ...
        {'real', 'finite', 'vector', 'numel', 8, 'nonnegative'});

    if any(x(:) < bounds(:,1) | x(:) > bounds(:,2))
        error('Current parameters must lie within their bounds.');
    end

    % Initialize candidate with current parameters
    y = x;

    % Randomly select one parameter
    idx = randi(numel(x));

    % Gaussian random-walk perturbation
    proposal = x(idx) + stepsize(idx) * randn;

    % Apply reflective boundary conditions
    lower = bounds(idx, 1);
    upper = bounds(idx, 2);
    width = upper - lower;

    % Periodic reflection handles arbitrarily large overshoots
    reflected = mod(proposal - lower, 2 * width);

    if reflected <= width
        y(idx) = lower + reflected;
    else
        y(idx) = upper - (reflected - width);
    end

end