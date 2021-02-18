

#===============================================================================
solve_pH_from_AT <- function (p_alktot, p_dicvar, p_bortot, p_po4tot, p_siltot, 
                             p_nh4tot, p_h2stot, p_so4tot, p_flutot, p_pHscale, p_dicsel, 
                             p_askVal=FALSE, p_dissoc, p_temp=18, p_sal=35, p_pres=0, p_hini)
#===============================================================================

# Determines [H+] from Total alkalinity and dissolved total elements in sea water. 
# Universal pH solver that converges from any given initial value,
# 
{
    #--------------------#
    # Argument variables #
    #--------------------#

    # p_dicvar           :   stands for one of the four carbonate system related variables,
    #                        depending on the value of p_dicsel:
    # p_dicsel = "DIC":   p_dicvar = DIC
    # p_dicsel = "CO2":   p_dicvar = [CO2*]
    # p_dicsel = "HCO3":  p_dicvar = [HCO3-]
    # p_dicsel = "CO3":   p_dicvar = [CO3--]
    # p_hini             :  (optionnal) initial value(s) for iterative algorithm 
    #                       if pdicsel = "CO3", it is an array of 2 values
    #                       else            it is a scalar

    # -----------------------------------------------------------------------
    
    
    i_dicsel = factor(toupper(p_dicsel), c("DIC", "CO2", "HCO3", "CO3"))
    if (is.na(i_dicsel))
        stop("Invalid parameter p_dicsel= ", p_dicsel)
    p_dicsel = toupper(p_dicsel)
    
    
    #-----------#
    # Constants #
    #-----------#


    # Threshold value for switching from
    # pH space to [H^+] space iterations.
    pz_exp_threshold = 1.0

    # Greatest acceptable ratio of variation for a Newton iterate zh relative to its
    # predecessor zh_prev:
    # exp(-pz_exp_upperlim)  <  zh/zh_prev  <  exp(pz_exp_upperlim)
    # i. e.,
    # abs(log(zh/zh_prev)) < pz_exp_upperlim
    pz_exp_upperlim  = 4.6   # exp(4.6) = 100.

    # Lowest limit for a Regula Falsi iterate
    # zh relative to zh_min, as a fraction of the length of the [zh_min, zh_max] interval:
    # zh_min + pz_rf_thetahmin * (zh_max - zh_min)
    #    <= zh <= zh_max - pz_rf_thetahmin * (zh_max - zh_min)
    pz_rf_thetamin   = 0.10
    pz_rf_thetamax   = 1.00 - pz_rf_thetamin

    # Maximum number of successive H_min/H_max changes
    jp_nmaxsucc      = 3

    # NaN for [H^+] results
    pp_hnan = -1.

    # Maximum number of iterations for each method
    jp_maxniter    = 100
    

    # -----------------------------------------------------------------------

    
    # Bookkeeping variable
    niter    = c(jp_maxniter+ 2, jp_maxniter+ 2)
    # pH value(s) that is(are)  solution
    p_solution = vector("numeric", 2)
    # value(s) of the rational function form of AT-pH equation at these pH
    p_value    = vector("numeric", 2)
    
    if (missing (p_dissoc)) 
    {
        p_dissoc <- compute_dissoc_constants (p_temp, p_sal, p_pres, p_pHscale)
    }

    api = list()

    # Compute conversion factor from [H+] free scale to chosen scale
    
    if (p_pHscale == "T") {
        # Conversion factor from [H+]free to [H+] total
        api$aphscale <- 1.0 + p_so4tot/p_dissoc$K_HSO4
    } else if (p_pHscale == "SWS") {
        # Conversion factor from [H+]free to [H+] sea-water
        api$aphscale <- 1.0 + p_so4tot/p_dissoc$K_HSO4 + + p_flutot/p_dissoc$K_HF
    } else {
        api$aphscale <- 1.0 
    }

    # Compute PI factors from dissociation constants K1, K2, Kb, Kw, ...
    
    # For each acid system A, 
    # - api1_aaa <-- K_A1
    # - api2_aaa <-- K_A1*K_A2
    # - api3_aaa <-- K_A1*K_A2*K_A3
    # - ...

    api$api1_dic <- p_dissoc$K1_DIC 
    api$api2_dic <- p_dissoc$K1_DIC * p_dissoc$K2_DIC
    api$api1_bor <- p_dissoc$K_BT
    api$api1_po4 <- p_dissoc$K1_PO4
    api$api2_po4 <- p_dissoc$K1_PO4 * p_dissoc$K2_PO4
    api$api3_po4 <- p_dissoc$K1_PO4 * p_dissoc$K2_PO4 * p_dissoc$K3_PO4
    api$api1_sil <- p_dissoc$K_Sil
    api$api1_nh4 <- p_dissoc$K_NH4
    api$api1_h2s <- p_dissoc$K_H2S
    api$api1_wat <- p_dissoc$K_H2O
    # K_HSO4 and K_HF given on free pH scale: convert to chosen scale
    api$api1_so4 <- p_dissoc$K_HSO4 * api$aphscale
    api$api1_flu <- p_dissoc$K_HF   * api$aphscale
        

    if (missing (p_hini)) {

        zh_ini = c(pp_hnan, pp_hnan)
    } else {

        zh_ini = p_hini

    }

    result = HINFSUPINI (p_alktot,  p_dicvar, p_bortot,
                        p_po4tot, p_siltot,
                        p_nh4tot, p_h2stot,
                        p_so4tot, p_flutot,
                        p_dicsel, api,  zh_ini)
    zh_inf   = result$p_hinf
    zh_sup   = result$p_hsup
    zh_ini   = result$p_hini
    k_nroots = result$k_nroots


    if (exists("DEBUG_PHSOLVERS") && DEBUG_PHSOLVERS)
    {
        print (c('[solve_pH_from_AT] n_roots  :', k_nroots))
        print (c('[solve_pH_from_AT] h_inf(:) :', zh_inf))
        print (c('[solve_pH_from_AT] h_sup(:) :', zh_sup))
        print (c('[solve_pH_from_AT] h_ini(:) :', zh_ini))
    }


    for (i_root in seq_len(k_nroots))
    {
        niter[i_root] = 0           # Reset counters of iterations

        zh_min = zh_inf[i_root]
        zh_max = zh_sup[i_root]

        zeqn_hmin = equation_at(p_alktot, zh_min,  
                                p_dicvar, p_bortot,
                                p_po4tot, p_siltot,
                                p_nh4tot, p_h2stot,
                                p_so4tot, p_flutot,
                                p_dicsel, api  )

        zeqn_hmax = equation_at(p_alktot, zh_max,  
                                p_dicvar, p_bortot,
                                p_po4tot, p_siltot,
                                p_nh4tot, p_h2stot,
                                p_so4tot, p_flutot,
                                p_dicsel, api  )

        zh     = zh_ini[i_root]

        if (abs(zeqn_hmin) < abs(zeqn_hmax)) {
            zh_absmin   = zh_min
            zeqn_absmin = zeqn_hmin
        }
        else {
            zh_absmin   = zh_max
            zeqn_absmin = zeqn_hmax
        }

        
        zh_prev = zh

        result = equation_at (p_alktot, zh,      
                            p_dicvar, p_bortot,
                            p_po4tot, p_siltot,
                            p_nh4tot, p_h2stot,
                            p_so4tot, p_flutot,
                            p_dicsel, api, TRUE )
        zeqn    = result[1]
        zdeqndh = result[2]

        nsucc_max = 0
        nsucc_min = 0

                                    # The second root, if any, is on an
        if (i_root == 2) {          # increasing branch of the EQUATION_AT
                                    # function. solve_pH_from_AT requires
            zeqn    = -zeqn         # that EQUATION_AT(..., zh_min, ...) > 0
            zdeqndh = -zdeqndh      # and  EQUATION_AT(..., zh_max, ...) < 0.
                                    # We therefore change the sign of the
        }                           # function for the second root.


        repeat {

            if (niter[i_root] > jp_maxniter) {
                zh = pp_hnan
                break
            }


            if (zeqn > 0.) {
                zh_min    = zh
                zeqn_hmin = zeqn
                nsucc_min = nsucc_min + 1
                nsucc_max = 0
            } else if (zeqn < 0.) {
                zh_max = zh
                zeqn_hmax = zeqn
                nsucc_max = nsucc_max + 1
                nsucc_min = 0
            } else {
                # zh is the root; unlikely but, one never knows
                break
            }


            # Now determine the next iterate zh
            niter[i_root] = niter[i_root] + 1



            if (abs(zeqn) >= 0.5*abs(zeqn_absmin)) {

                # If the function evaluation at the current point is not
                # decreasing faster than expected with a bisection step
                # (at least linearly) in absolute value take one regula falsi
                # step, except if either the minimum or the maximum value has
                # been modified more than three times (default - can be
                # overridden by modifying the paraemeter value jp_nmaxsucc)
                # in a row. This increases chances that the maximum, resp.
                # minimum, is also revised from time to time.
                
                if ((nsucc_min > jp_nmaxsucc) || (nsucc_max > jp_nmaxsucc)) {

                    # Bisection step in pH-space:
                    # ph_new = (ph_min + ph_max)/2d0
                    # In terms of [H]_new:
                    # [H]_new = 10**(-ph_new)
                    #         = 10**(-(ph_min + ph_max)/2d0)
                    #         = sqrt(10**(-(ph_min + phmax)))
                    #         = sqrt(zh_max * zh_min)
                    zh = sqrt(zh_max * zh_min)
                    nsucc_min = 0
                    nsucc_max = 0

                } else {

                    # Regula falsi step  on [H_min, H_max] (too expensive on [pH_min, pH_max])
                    zrf_hmin = -zeqn_hmax/(zeqn_hmin - zeqn_hmax)
                    zrf_hmax =  zeqn_hmin/(zeqn_hmin - zeqn_hmax)

                    if (zrf_hmin < pz_rf_thetamin) {
                        zh = pz_rf_thetamin * zh_min + pz_rf_thetamax * zh_max
                    } else if (zrf_hmin >  pz_rf_thetamax) {
                        zh = pz_rf_thetamax * zh_min + pz_rf_thetamin * zh_max
                    } else {
                        zh = zrf_hmin*zh_min + zrf_hmax*zh_max
                    }

                }


            } else {

                # dzeqn/dpH = dzeqn/d[H] * d[H]/dpH
                #           = -zdeqndh * log(10) * [H]
                # \Delta pH = -zeqn/(zdeqndh*d[H]/dpH) = zeqn/(zdeqndh*[H]*log(10))

                # pH_new = pH_old + \deltapH

                # [H]_new = 10**(-pH_new)
                #         = 10**(-pH_old - \Delta pH)
                #         = [H]_old * 10**(-zeqn/(zdeqndh*[H]_old*log(10)))
                #         = [H]_old * exp(-log(10)*zeqn/(zdeqndh*[H]_old*log(10)))
                #         = [H]_old * exp(-zeqn/(zdeqndh*[H]_old))

                zh_lnfactor = -zeqn/(zdeqndh*zh_prev)

                if (abs(zh_lnfactor) < pz_exp_threshold) {
                    zh_delta    = zh_lnfactor*zh_prev
                    zh          = zh_prev + zh_delta
                } else if (abs(zh_lnfactor) < pz_exp_upperlim) {
                    zh          = zh_prev*exp(zh_lnfactor)
                } else {
                    zh_lnfactor = pz_exp_upperlim * sign(zh_lnfactor)
                    zh          = zh_prev*exp(zh_lnfactor)
                }


                if ( ( zh < zh_min ) || ( zh > zh_max ) ) {
                    # if [H]_new < [H]_min or [H]_new > [H]_max, then take
                    # one regula falsi step on [H_min, H_max]
                    zrf_hmin = -zeqn_hmax/(zeqn_hmin - zeqn_hmax)
                    zrf_hmax =  zeqn_hmin/(zeqn_hmin - zeqn_hmax)

                    if (zrf_hmin < pz_rf_thetamin) {
                        zh = pz_rf_thetamin * zh_min + pz_rf_thetamax * zh_max
                    } else if (zrf_hmin >  pz_rf_thetamax) {
                        zh = pz_rf_thetamax * zh_min + pz_rf_thetamin * zh_max
                    } else {
                        zh = zrf_hmin*zh_min + zrf_hmax*zh_max
                    }

                }

            }


            if (abs(zeqn_absmin) > abs(zeqn)) 
            {
                if (exists("DEBUG_PHSOLVERS") && DEBUG_PHSOLVERS)
                {
                    print (c('[solve_pH_from_AT] adjusting absmin     :', zh_prev, zeqn))
                }
                zh_absmin   = zh_prev
                if (i_root == 2) {
                    zeqn_absmin = -zeqn
                } else {
                    zeqn_absmin =  zeqn
                }
            }


            result = equation_at(p_alktot, zh,      
                            p_dicvar, p_bortot,
                            p_po4tot, p_siltot,
                            p_nh4tot, p_h2stot,
                            p_so4tot, p_flutot,
                            p_dicsel, api, TRUE )
            zeqn    = result[1]
            zdeqndh = result[2]

                                            # Exit if the length of [H_min, H_max] is of
                                            # the same order as the required precision
            if ((zh_max - zh_min) < (0.5*(zh_max + zh_min) * pp_rdel_ah_target)) {

                                            # Check if the previously registered
                                            # absolute minimum eqn value was lower than
                                            # the current one - if so, return that one.
                if (abs(zeqn_absmin) < abs(zeqn)) {
                    zh   = zh_absmin
                    zeqn = zeqn_absmin
                }

                break                          # Done!
            }

                                            # Prepare the next iteration
            if (i_root == 2) { 
                zeqn    = -zeqn             #  - revert equation signs if we
                zdeqndh = -zdeqndh          #    are searching for the second root
            }

            zh_prev = zh                    #  - save the previous iterate

        }


        p_solution[i_root] = zh

        if (p_askVal) {

            if (zh != pp_hnan) { 
                p_value[i_root] = zeqn
            } else {
                p_value[i_root] = .Machine$double.xmax   #Biggest floating point number
            }

        }

    }   # end for (i_root in seq_len(k_nroots))


                                        # Finally, initialize the variables
                                        # related to the non-used roots.
    while (i_root < 2) {

        i_root = i_root + 1
        niter[i_root]               = 0

        p_solution[i_root] = pp_hnan
        if (p_askVal) p_value[i_root] = .Machine$double.xmax   #Biggest floating point number

    }

    if (p_dicsel != "CO3")
    {
        # There is at most one root : no need to return two values
        p_solution = p_solution[1]
        p_value    = p_value[1]
    }
    
    if (p_askVal) {
        result = data.frame(zh = p_solution, val = p_value)
        return (result)
    } else {
        return(p_solution)
    }
}
