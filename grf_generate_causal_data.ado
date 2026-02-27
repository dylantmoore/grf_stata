*! grf_generate_causal_data.ado -- Generate synthetic causal data for benchmarking
*! Version 0.1.0
*! Implements DGPs from the R grf package for Monte Carlo experiments

program define grf_generate_causal_data, rclass
    version 14.0

    syntax , N(integer) P(integer)  ///
        [                           ///
            DGP(string)             ///
            SEED(integer -1)        ///
        ]

    /* ---- Defaults ---- */
    if "`dgp'" == "" {
        local dgp "simple"
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
    local valid_dgps "simple aw1 aw2 aw3 ai1 ai2 kunzel nw1 nw2 nw3 nw4"
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

    /* Some DGPs need minimum p */
    if "`dgp'" == "aw2" | "`dgp'" == "aw3" | "`dgp'" == "nw2" ///
        | "`dgp'" == "nw3" | "`dgp'" == "nw4" {
        if `p' < 2 {
            display as error "dgp(`dgp') requires p >= 2"
            exit 198
        }
    }
    if "`dgp'" == "simple" | "`dgp'" == "kunzel" {
        if `p' < 3 {
            display as error "dgp(`dgp') requires p >= 3"
            exit 198
        }
    }
    if "`dgp'" == "nw1" {
        if `p' < 5 {
            display as error "dgp(nw1) requires p >= 5"
            exit 198
        }
    }

    /* ---- Set seed ---- */
    if `seed' >= 0 {
        set seed `seed'
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "Generating Causal Data"
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
    /* Most DGPs use uniform covariates; some use normal */
    local x_dist "uniform"
    if "`dgp'" == "nw1" | "`dgp'" == "nw2" | "`dgp'" == "nw3" | "`dgp'" == "nw4" {
        local x_dist "normal"
    }

    forvalues j = 1/`p' {
        if "`x_dist'" == "uniform" {
            quietly gen double X`j' = runiform()
        }
        else {
            quietly gen double X`j' = rnormal()
        }
    }

    /* ---- Generate treatment, outcome, and true CATE by DGP ---- */

    if "`dgp'" == "simple" {
        /* Y = max(X1+X2+X3, 0) + eps, tau(x) = X1+X2, W ~ Bernoulli(0.5) */
        quietly {
            gen double tau_true = X1 + X2
            gen byte W = rbinomial(1, 0.5)
            gen double eps = rnormal()
            gen double mu = max(X1 + X2 + X3, 0)
            gen double Y = mu + tau_true * W + eps
            drop eps mu
        }
        local dgp_desc "Simple: Y = max(X1+X2+X3,0) + tau*W + eps, tau = X1+X2"
    }

    else if "`dgp'" == "aw1" {
        /* Athey-Wager 1: tau = 1 + 1/(1+exp(-20*(X1-1/3))) */
        quietly {
            gen double e_x = 0.5
            gen double mu = 2 * X1 - 1
            gen double tau_true = 1 + 1 / (1 + exp(-20 * (X1 - 1/3)))
            gen byte W = rbinomial(1, e_x)
            gen double eps = rnormal()
            gen double Y = mu + tau_true * W + eps
            drop eps mu e_x
        }
        local dgp_desc "Athey-Wager 1: tau = 1 + 1/(1+exp(-20*(X1-1/3)))"
    }

    else if "`dgp'" == "aw2" {
        /* Athey-Wager 2: tau = 1 + 1/(1+exp(-20*(X1-1/3))) * 1/(1+exp(-20*(X2-1/3))) */
        quietly {
            gen double e_x = 0.5
            gen double mu = 2 * X1 - 1
            gen double tau_true = 1 ///
                + 1 / (1 + exp(-20 * (X1 - 1/3))) ///
                * 1 / (1 + exp(-20 * (X2 - 1/3)))
            gen byte W = rbinomial(1, e_x)
            gen double eps = rnormal()
            gen double Y = mu + tau_true * W + eps
            drop eps mu e_x
        }
        local dgp_desc "Athey-Wager 2: tau with X1 and X2 interaction"
    }

    else if "`dgp'" == "aw3" {
        /* Athey-Wager 3: tau = (X1>0.5)*(X2>0.5) */
        quietly {
            gen double e_x = 0.5
            gen double mu = 2 * X1 - 1
            gen double tau_true = (X1 > 0.5) * (X2 > 0.5)
            gen byte W = rbinomial(1, e_x)
            gen double eps = rnormal()
            gen double Y = mu + tau_true * W + eps
            drop eps mu e_x
        }
        local dgp_desc "Athey-Wager 3: tau = I(X1>0.5)*I(X2>0.5)"
    }

    else if "`dgp'" == "ai1" {
        /* Athey-Imbens 1: heterogeneous propensity */
        quietly {
            gen double e_x = 0.25 * (1 + dbeta(X1, 2, 4))
            gen double mu = 2 * X1 - 1
            gen double tau_true = 1 + 1 / (1 + exp(-20 * (X1 - 1/3)))
            gen byte W = rbinomial(1, e_x)
            gen double eps = rnormal()
            gen double Y = mu + tau_true * W + eps
            drop eps mu e_x
        }
        local dgp_desc "Athey-Imbens 1: heterogeneous propensity"
    }

    else if "`dgp'" == "ai2" {
        /* Athey-Imbens 2: non-constant variance */
        quietly {
            gen double e_x = 0.5
            gen double mu = 2 * X1 - 1
            gen double tau_true = 1 + 1 / (1 + exp(-20 * (X1 - 1/3)))
            gen byte W = rbinomial(1, e_x)
            gen double sigma = 1 + 2 * X1
            gen double eps = rnormal() * sigma
            gen double Y = mu + tau_true * W + eps
            drop eps mu e_x sigma
        }
        local dgp_desc "Athey-Imbens 2: non-constant variance"
    }

    else if "`dgp'" == "kunzel" {
        /* Kunzel et al. (2019): complex response surface */
        quietly {
            gen double e_x = 0.5
            gen double mu_0 = sin(_pi * X1 * X2) ///
                + 2 * (X3 - 0.5)^2
            gen double tau_true = (X1 + X2) / 2
            gen byte W = rbinomial(1, e_x)
            gen double eps = rnormal()
            gen double Y = mu_0 + tau_true * W + eps
            drop eps mu_0 e_x
        }
        local dgp_desc "Kunzel: mu0 = sin(pi*X1*X2) + 2*(X3-0.5)^2"
    }

    else if "`dgp'" == "nw1" {
        /* Nie-Wager 1: complex mu, moderate tau */
        quietly {
            gen double e_x = 0.5
            gen double mu = 2 * X1 - 1 + X2 + X3 + X4 + X5
            gen double tau_true = 0.5 * X1
            gen byte W = rbinomial(1, e_x)
            gen double eps = rnormal()
            gen double Y = mu + tau_true * W + eps
            drop eps mu e_x
        }
        local dgp_desc "Nie-Wager 1: mu = 2*X1-1+X2+X3+X4+X5, tau = 0.5*X1"
    }

    else if "`dgp'" == "nw2" {
        /* Nie-Wager 2: non-linear response, interaction */
        quietly {
            gen double e_x = 1 / (1 + exp(-X1))
            gen double mu = max(X1 + X2, 0) + max(X1 + X2, 0)
            gen double tau_true = X1 + log(1 + exp(X2))
            gen byte W = rbinomial(1, e_x)
            gen double eps = rnormal()
            gen double Y = mu + tau_true * W + eps
            drop eps mu e_x
        }
        local dgp_desc "Nie-Wager 2: confounded propensity, non-linear tau"
    }

    else if "`dgp'" == "nw3" {
        /* Nie-Wager 3: cosine response */
        quietly {
            gen double e_x = 1 / (1 + exp(-X1) + exp(-X2))
            gen double mu = cos(X1 + X2)
            gen double tau_true = 1
            gen byte W = rbinomial(1, e_x)
            gen double eps = rnormal()
            gen double Y = mu + tau_true * W + eps
            drop eps mu e_x
        }
        local dgp_desc "Nie-Wager 3: constant tau=1, cos response"
    }

    else if "`dgp'" == "nw4" {
        /* Nie-Wager 4: polynomial response */
        quietly {
            gen double e_x = 1 / (1 + exp(-X1) + exp(-X2))
            gen double mu = (1 + 1 / (1 + exp(-20 * (X1 - 1/3)))) ///
                * (1 + 1 / (1 + exp(-20 * (X2 - 1/3))))
            gen double tau_true = ///
                1 / (1 + exp(-20 * (X1 - 1/3))) ///
                * 1 / (1 + exp(-20 * (X2 - 1/3)))
            gen byte W = rbinomial(1, e_x)
            gen double eps = rnormal()
            gen double Y = mu + tau_true * W + eps
            drop eps mu e_x
        }
        local dgp_desc "Nie-Wager 4: sigmoid tau with confounded propensity"
    }

    /* ---- Label variables ---- */
    label variable Y "Outcome"
    label variable W "Treatment indicator"
    label variable tau_true "True CATE"
    forvalues j = 1/`p' {
        label variable X`j' "Covariate `j'"
    }

    /* ---- Summary statistics ---- */
    quietly summarize Y
    local y_mean = r(mean)
    local y_sd = r(sd)
    quietly summarize tau_true
    local tau_mean = r(mean)
    local tau_sd = r(sd)
    quietly summarize W
    local w_mean = r(mean)

    /* ---- Display results ---- */
    display as text "Data generated successfully."
    display as text ""
    display as text "DGP: `dgp_desc'"
    display as text ""
    display as text "Variables created:"
    display as text "  X1 ... X`p'   covariates"
    display as text "  W             treatment (mean = " as result %5.3f `w_mean' as text ")"
    display as text "  Y             outcome   (mean = " as result %7.4f `y_mean' ///
        as text ", sd = " as result %7.4f `y_sd' as text ")"
    display as text "  tau_true      true CATE (mean = " as result %7.4f `tau_mean' ///
        as text ", sd = " as result %7.4f `tau_sd' as text ")"
    display as text ""

    /* ---- Store r() results ---- */
    return scalar N          = `n'
    return scalar p          = `p'
    return scalar tau_mean   = `tau_mean'
    return scalar tau_sd     = `tau_sd'
    return scalar y_mean     = `y_mean'
    return scalar y_sd       = `y_sd'
    return scalar w_mean     = `w_mean'
    return local  dgp          "`dgp'"
    return local  dgp_desc     "`dgp_desc'"
end
