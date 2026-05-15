function tau_bar = calibrate_tau_bar(target_share)
% Calibrate tau_bar so that (tau * co2_agg / Y) = target_share in steady state.

  % Use Dynare globals that are created when you run the .mod file:
  global M_ oo_

  % Initial guess for per-ton price (goods per ton)
  set_param_value('tau_bar', 0.04);

  for it = 1:50
      steady;   % solve SS with current tau_bar

      % Pull SS values:
      Y_ss   = oo_.steady_state(strmatch('Y',       M_.endo_names,'exact'));
      CO2_ss = oo_.steady_state(strmatch('co2_agg', M_.endo_names,'exact'));
      tau_ss = oo_.steady_state(strmatch('tau',     M_.endo_names,'exact'));  % = tau_bar at SS

      share  = (tau_ss * CO2_ss) / Y_ss;

      if ~isfinite(share) || share<=0
          error('Calib error: nonpositive/ill-defined share in SS. Check SS first.');
      end

      if abs(share - target_share) < 1e-12
          break;
      end

      % Multiplicative update (stable, scale-free)
      set_param_value('tau_bar', get_param_by_name('tau_bar') * (target_share / share));
  end

  tau_bar = get_param_by_name('tau_bar');
  fprintf('Calibrated tau_bar (good/ton) = %.8f; SS share = %.8f\n', tau_bar, share);
end
