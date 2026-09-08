# L1 closure card

## Review grades

- **H-FN1_soft_modes**: {'grade': 'computed_galerkin', 'accepted': False, 'evidence': [0.051857765435603144, 0.7202379026698027, 1.3240904405927962, 2.051758165237811, 2.3674631990072488, 2.685780384381832, 3.0113008064065614, 3.2459814701560292], 'limit': 'Galerkin 2N basis, not continuum Goldstone count'}
- **G_ij_from_hessian**: {'grade': 'computed', 'g0_eff_meV': 0.01853792309120878, 'rel_frobenius': 15.729136554439238, 'limit': 'Amplitude-block Hessian off-diagonals; not torsion'}
- **H-CAT1**: {'grade': 'computed_ode', 'accepted': True, 'ratio': 65.36567842728033, 'limit': 'γ∝ω² ansatz on Galerkin eigenvalues'}
- **H-TOR1**: {'grade': 'ansatz_future', 'accepted': False, 'parked': True, 'rel_frobenius': 12.138311467779474, 'limit': '2D proxy parked; future 3D Weitzenböck, not used for G_ij'}
- **CS_EOM**: {'grade': 'covariant_stepper', 'a_in_stepper': True, 'flux_rel_err': 0.00037990626235976917, 'j_cancellation_ratio': 1.001584725766571, 'limit': 'background a in kinetic EOM; Hall remains algebra', 'analytic_attachment_rms': 3.81174592563034e-16}
- **C1_quantum_metric**: {'grade': 'toy_band_physical_branch', 'chern': 0.9925353199430855, 'kappa_geom_meV': 0.04297708289033414, 'KT_hit_Tstar': False, 'K_fold_Tstar': 5.437513778749539, 'limit': 'QWZ C=1, not MATBG continuum; n0 ansatz'}
- **Track_B**: {'grade': 'matched_not_hubbard', 't_over_U': 1.0, 't_lambda_over_U': 0.37869822485207105, 'U_over_W': 0.11142857142857142}
- **Villain_C7**: {'grade': 'checkerboard_xy_scan', 'eta_uv_fit': 0.05001242264844266, 'eta_KT_fit': 2.3950552656161124, 'limit': 'Finite-L XY, not thermodynamic Villain'}
- **Hopf_integral**: {'grade': 'S3_whitehead_derived', 'Q_S3': 0.9997991943200188, 'Q_stereo': -0.5942549114759039, 'Q_spinor': -0.38292713514304394, 'Q1': 0.20120864812853712, 'QN': 0.004284310029988481, 'lock': 0.2, 'limit': 'string-free FFT; lock remains CS/linking 1/N'}

## Locks
{
  "g0_star_meV": 0.09816962498385381,
  "Delta0_meV": 0.6,
  "E_vac_meV": 0.12870089825559464,
  "lambda_star": 0.7627734251620608,
  "Q_H": 0.2,
  "N": 5
}
