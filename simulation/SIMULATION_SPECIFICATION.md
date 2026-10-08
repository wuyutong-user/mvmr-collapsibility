# Simulation Specification Table — Frozen v1.1

This file is the implementation contract for the revised simulation study. The production code must not silently alter these settings after the strength pilot has been frozen.

| Domain | Frozen specification |
|---|---|
| Main aims | (1) Exact structural theorem calibration; (2) adjustment-related bias and precision; (3) empirical type I error, power and coverage; (4) conditional estimability; (5) robustness to sample overlap, instrument geometry and practical instrument selection |
| Roles | Six Fig.1 roles only: confounder; collider; independent cause; upstream surrogate of exposure; downstream surrogate of outcome; downstream surrogate of exposure |
| Scientific target | Total causal X→Y target; beta_X = 0 or 0.10 in the primary grid; beta_X = 0.05 is a prespecified sensitivity |
| Structural path coefficients | 0.30 whenever a Fig.1 non-target directed edge is present |
| Common U | U enters X, Z and Y with coefficient 0.30. In binary outcomes U explicitly enters the logistic linear predictor |
| Residuals | Independent N(0,1) residuals for continuous X/Z/Y equations as applicable |
| SNP generation | Independent biallelic SNPs, primary MAF=0.30 |
| SNP total | K=100 in every primary instrument configuration |
| Direct genetic source classes | `direct_X`, `direct_Z`, `direct_shared`; used only for the individual-level six-DAG DGP |
| Theorem summary classes | `X_only`, `Z_only`, `shared`; defined by population SNP-X/SNP-Z summary associations, not by biological/direct source |
| Primary direct-source configs | 30/30/40; 0/60/40; 60/0/40; 50/50/0 for direct-X/direct-Z/direct-shared |
| Exact theorem summary configs | 30/30/40; 0/60/40; 60/0/40; 50/50/0 for X-only/Z-only/shared |
| Primary direct shared alignment | +0.80 in the direct genetic source construction; the achieved summary-space rho is stored rather than assumed equal to +0.80 |
| Alignment sensitivities | direct-source correlation 0, -0.80, +0.95; exact theorem module separately constructs weighted summary-space rho 0, -0.80, +0.80, +0.95 |
| Theorem same set | Exact identical SNP rows, outcome vector and W in unadjusted and adjusted regressions |
| Strict finite-sample same set | All K=100 prespecified SNPs must have finite GWAS estimates. A failed binary SNP-outcome fit makes the replicate fail rather than silently reducing K |
| Practical set sensitivity | Within each replicate, UVMR uses SNPs significant for X and MVMR uses the union significant for X or Z; this is explicitly secondary to the strict same-set analysis |
| Primary sample design | Non-overlapping two-sample: one exposure sample estimates SNP-X and SNP-Z; an independent outcome sample estimates SNP-Y |
| Exposure covariance | Because X and Z are estimated in the same exposure GWAS sample, SNP-X/SNP-Z sampling covariance is estimated and used in conditional-F calculations |
| Sample-overlap sensitivity | 100% same-sample; common random numbers pair this sensitivity with the primary exposure-sample realization |
| Primary N | N_E=N_Y=10,000 |
| N sensitivity | 5,000 and 50,000 |
| Primary binary prevalence | approximately 10%; scenario-level logistic intercept is calibrated once using a reference population and then frozen |
| Binary prevalence sensitivity | 50% |
| Primary estimator | No-intercept inverse-outcome-variance weighted summary regression; point estimates exactly match the theorem regression |
| Primary finite-sample SE | Weighted-lm SE, matching `summary(lm(..., weights=1/seY^2))` and the MVMR `TwoSampleMR::mv_multiple` convention; fixed-effect SE is also retained |
| Primary p-value | Two-sided normal-reference p-value based on the primary SE; t-reference p-value is retained as an audit field |
| Primary Theorem-1 check | response = alpha_Z × delta_XZ; identity error < 1e-10 |
| Six-DAG continuous population check | gamma_Y = beta_X gamma_X + alpha_Z gamma_Z with alpha_Z=0.30 for confounder/independent cause and 0 for the other four roles |
| Binary realistic simulation | Individual-level logistic outcome; marginal SNP-Y logORs estimated by SNP-wise logistic regression; no rare-outcome approximation |
| Binary interpretation | Finite-sample binary error is target-relative total-logOR error and may include scale-dependent variation; structural and scale components are not claimed to be separately identifiable from the logistic simulation |
| Binary exact theorem check | Separate summary-level module constructs ell = S^-1(beta* x + alpha* z)+r_scale and verifies delta_total=delta_struct+delta_scale |
| Scale factors | S=0.75, 1.00, 1.50 |
| Binary residual modes | none; scale-only; reinforce; partial-offset; exact-offset |
| Precision outputs | empirical variance; empirical SD; replicate model SE; mean/median model SE; variance ratio Var(adj)/Var(unadj); relative efficiency; RMSE |
| Performance outputs | continuous bias or binary target-relative total-logOR error; RMSE; 95% coverage |
| Inference outputs | empirical type I error under beta_X=0; empirical power under beta_X≠0; Wilson 95% Monte Carlo intervals |
| Reliability | primary Eq.7 F_X|Z, F_Z|X, minimum F; package-compatible MVMR crosswalk; weighted rho; delta_XZ; condition number |
| Conditional F primary | Sanderson-Spiller-Bowden Eq.7: F_TS,k = Q_xk/[L-(K-1)], giving Q/(L-1) for two exposures; nuisance delta obtained by no-intercept IVW |
| Conditional F package crosswalk | Also stores the public `MVMR::strength_mvmr` convention (no-intercept OLS nuisance delta and Q/L) for reproducibility/crosswalk |
| Conditional F role | Reliability only; never assigns a causal role or adjustment class |
| Strength pilot | Before production, one prespecified reference scenario evaluates an RMS grid using common random numbers; values closest to median min-F targets 25, 10 and 5 are frozen to an RDS file |
| Primary strength | Production loads the frozen RMS closest to target median min-F≈25. No retuning after final results are inspected |
| Strength sensitivity | Frozen RMS values closest to conditional-F targets 25, 10 and 5 |
| Primary replicates | 2,000 continuous; 1,000 binary per scenario |
| Sensitivity replicates | 1,000 per scenario |
| Primary scenario count | 6 roles × 2 outcomes × 2 beta values × 4 direct-source configurations = 96 scenarios |
| RNG | Deterministic scenario/replicate seeds. Common random numbers are used across controlled beta/strength/alignment/N/prevalence/sample-overlap contrasts. Results do not depend on parallel scheduling order |
| HPC | One persistent PSOCK cluster is reused across scenarios. BLAS/OpenMP threads are forced to 1 inside workers to avoid 100×nested oversubscription |
| Checkpoints | Replicates are checkpointed in chunks (default 100); production supports deterministic resume after interruption |
| Failure handling | No silent deletion. Each failed replicate records a reason; >1% failure is flagged/warned and must be investigated |
| Raw output | Replicate-level RDS files sufficient to reconstruct all summaries |
| Summary output | Separate primary and sensitivity CSV/RDS summaries; empirical SD and model SE are never conflated |
| Egger | Optional sensitivity only. It is not a direct validation of Theorem 1 because it adds an intercept |
| Production gate | Full production will not run without a successful server-side parse/validation/smoke preflight performed after the frozen strength calibration |
