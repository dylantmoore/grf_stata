*! grf_generate_causal_survival_data.ado -- Synthetic right-censored survival data
*! Version 0.1.0
*! Generates survival data with known treatment effects for benchmarking

program define grf_generate_causal_survival_data, rclass
    version 14.0

    syntax , N(integer) P(integer)  ///
        [                           ///
            DGP(string)             ///
            SEED(integer -1)        ///
        ]

    /* ---- Defaults ---- */
    if "`dgp'" == "" {
        local dgp "simple1"
    }
    local dgp = lower("`dgp'")

    /* ---- Validate inputs ---- */
    if `n' < 10 {
        display as error "n() must be at least 10"
        exit 198
    }
    if `p' < 1 {
        display as error "p() must be at least 1"
        exit 198
    }

    /* Validate DGP name */
    local valid_dgps "simple1 type1 type2 type3 type4 type5"
    local dgp_ok 0
    foreach d of local valid_dgps {
        if "`dgp'" == "`d'" {
            local dgp_ok 1
        }
    }
    if !`dgp_ok' {
        display as error "dgp(`dgp') not recognized"
        display as error "valid options: `valid_dgps'"
        exit 198
    }

    /* Minimum p requirements */
    if "`dgp'" == "type3" | "`dgp'" == "type4" | "`dgp'" == "type5" {
        if `p' < 2 {
            display as error "dgp(`dgp') requires p >= 2"
            exit 198
        }
    }

    /* ---- Set seed ---- */
    if `seed' >= 0 {
        set seed `seed'
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "Generating Causal Survival Data"
    display as text "{hline 55}"
    display as text "DGP:                   " as result "`dgp'"
    display as text "Observations:          " as result `n'
    display as text "Covariates:            " as result `p'
    if `seed' >= 0 {
        display as text "Seed:                  " as result `seed'
    }
    display as text "{hline 55}"
    display as text ""

    /* ---- Clear and set observation count ---- */
    drop _all
    quietly set obs `n'

    /* ---- Generate covariates X1..Xp ---- */
    forvalues j = 1/`p' {
        quietly gen double X`j' = runiform()
    }

    /* ---- Generate data by DGP ---- */

    if "`dgp'" == "simple1" {
        /* Basic proportional hazards with heterogeneous treatment effect
         * Baseline hazard: lambda_0 = 0.1
         * Covariate effect: exp(beta * X1)
         * Treatment effect: exp(gamma(X) * W) where gamma(X) = -0.5 * X1
         * T = -log(U) / lambda, C ~ Exp(0.1)
         */
        quietly {
            gen byte W = rbinomial(1, 0.5)
            gen double gamma_x = -0.5 * X1
            gen double tau_true = gamma_x
            gen double lambda = 0.1 * exp(0.5 * X1 + gamma_x * W)
            gen double U_surv = runiform()
            gen double T_latent = -log(U_surv) / lambda
            gen double C = -log(runiform()) / 0.1
            gen double T = min(T_latent, C)
            gen byte D = (T_latent <= C)
            drop gamma_x lambda U_surv T_latent C
        }
        local dgp_desc "Simple PH: lambda=0.1*exp(0.5*X1+gamma*W), gamma=-0.5*X1"
    }

    else if "`dgp'" == "type1" {
        /* Weibull baseline hazard with linear treatment heterogeneity
         * Baseline: Weibull with shape=2, scale depends on X1
         * Treatment: additive effect on log-hazard, heterogeneous in X1
         */
        quietly {
            gen byte W = rbinomial(1, 0.5)
            gen double tau_true = 0.5 * (X1 - 0.5)
            gen double log_lambda = -1 + 0.5 * X1 + tau_true * W
            gen double lambda = exp(log_lambda)
            gen double U_surv = runiform()
            /* Weibull: T = (-log(U)/lambda)^(1/shape) with shape=2 */
            gen double T_latent = (-log(U_surv) / lambda)^(0.5)
            gen double C = -log(runiform()) / 0.15
            gen double T = min(T_latent, C)
            gen byte D = (T_latent <= C)
            drop log_lambda lambda U_surv T_latent C
        }
        local dgp_desc "Type 1: Weibull baseline, tau = 0.5*(X1-0.5)"
    }

    else if "`dgp'" == "type2" {
        /* Log-normal survival with non-linear treatment effect
         * T = exp(mu + sigma*Z) where mu depends on X and W
         */
        quietly {
            gen byte W = rbinomial(1, 0.5)
            gen double tau_true = 1 / (1 + exp(-5 * (X1 - 0.5)))
            gen double mu = 1 + 0.5 * X1 + tau_true * W
            gen double T_latent = exp(mu + 0.5 * rnormal())
            gen double C = -log(runiform()) / 0.08
            gen double T = min(T_latent, C)
            gen byte D = (T_latent <= C)
            drop mu T_latent C
        }
        local dgp_desc "Type 2: Log-normal, tau = sigmoid(X1)"
    }

    else if "`dgp'" == "type3" {
        /* Proportional hazards with interaction treatment effect
         * Treatment heterogeneity depends on X1 and X2
         */
        quietly {
            gen byte W = rbinomial(1, 0.5)
            gen double tau_true = X1 * X2
            gen double lambda = 0.2 * exp(0.3 * X1 + 0.3 * X2 ///
                + tau_true * W)
            gen double U_surv = runiform()
            gen double T_latent = -log(U_surv) / lambda
            gen double C = -log(runiform()) / 0.15
            gen double T = min(T_latent, C)
            gen byte D = (T_latent <= C)
            drop lambda U_surv T_latent C
        }
        local dgp_desc "Type 3: PH with interaction tau = X1*X2"
    }

    else if "`dgp'" == "type4" {
        /* Heterogeneous propensity and treatment effect
         * Propensity depends on X1, treatment effect on X2
         * Creates confounding
         */
        quietly {
            gen double e_x = 0.25 + 0.5 * X1
            gen byte W = rbinomial(1, e_x)
            gen double tau_true = 0.5 + 0.5 * X2
            gen double lambda = 0.1 * exp(0.5 * X1 + 0.5 * X2 ///
                + tau_true * W)
            gen double U_surv = runiform()
            gen double T_latent = -log(U_surv) / lambda
            gen double C = -log(runiform()) / 0.12
            gen double T = min(T_latent, C)
            gen byte D = (T_latent <= C)
            drop e_x lambda U_surv T_latent C
        }
        local dgp_desc "Type 4: Confounded propensity, tau = 0.5+0.5*X2"
    }

    else if "`dgp'" == "type5" {
        /* Step function treatment effect with informative censoring
         * Treatment effect is a step function of X1
         * Censoring rate depends on X2 (informative censoring)
         */
        quietly {
            gen byte W = rbinomial(1, 0.5)
            gen double tau_true = cond(X1 > 0.5, 1, 0) ///
                * cond(X2 > 0.5, 0.5, -0.5)
            gen double lambda = 0.15 * exp(0.4 * X1 + 0.4 * X2 ///
                + tau_true * W)
            gen double U_surv = runiform()
            gen double T_latent = -log(U_surv) / lambda
            /* Informative censoring: censoring rate depends on X2 */
            gen double cens_rate = 0.1 + 0.2 * X2
            gen double C = -log(runiform()) / cens_rate
            gen double T = min(T_latent, C)
            gen byte D = (T_latent <= C)
            drop lambda U_surv T_latent C cens_rate
        }
        local dgp_desc "Type 5: Step tau, informative censoring"
    }

    /* ---- Label variables ---- */
    label variable T "Observed survival time"
    label variable D "Event indicator (1=event, 0=censored)"
    label variable W "Treatment indicator"
    label variable tau_true "True treatment effect"
    forvalues j = 1/`p' {
        label variable X`j' "Covariate `j'"
    }

    /* ---- Summary statistics ---- */
    quietly summarize T
    local t_mean = r(mean)
    local t_sd = r(sd)
    quietly summarize tau_true
    local tau_mean = r(mean)
    local tau_sd = r(sd)
    quietly summarize W
    local w_mean = r(mean)
    quietly count if D == 1
    local n_events = r(N)
    local n_censored = `n' - `n_events'
    local cens_pct = 100 * `n_censored' / `n'

    /* ---- Display results ---- */
    display as text "Data generated successfully."
    display as text ""
    display as text "DGP: `dgp_desc'"
    display as text ""
    display as text "Variables created:"
    display as text "  X1 ... X`p'   covariates"
    display as text "  W             treatment  (mean = " as result %5.3f `w_mean' as text ")"
    display as text "  T             obs. time  (mean = " as result %7.4f `t_mean' ///
        as text ", sd = " as result %7.4f `t_sd' as text ")"
    display as text "  D             event ind. (events = " as result `n_events' ///
        as text ", censored = " as result `n_censored' ///
        as text " [" as result %4.1f `cens_pct' as text "%])"
    display as text "  tau_true      true tau   (mean = " as result %7.4f `tau_mean' ///
        as text ", sd = " as result %7.4f `tau_sd' as text ")"
    display as text ""

    /* ---- Store r() results ---- */
    return scalar N          = `n'
    return scalar p          = `p'
    return scalar tau_mean   = `tau_mean'
    return scalar tau_sd     = `tau_sd'
    return scalar t_mean     = `t_mean'
    return scalar t_sd       = `t_sd'
    return scalar w_mean     = `w_mean'
    return scalar n_events   = `n_events'
    return scalar n_censored = `n_censored'
    return scalar cens_pct   = `cens_pct'
    return local  dgp          "`dgp'"
    return local  dgp_desc     "`dgp_desc'"
end
