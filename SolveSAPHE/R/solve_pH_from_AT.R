#===============================================================================
solve_pH_from_AT <- function (p_alktot, p_dictot, p_bortot, p_po4tot, p_siltot, 
                              p_nh4tot, p_h2stot, p_so4tot, p_flutot, p_pHscale,
                              p_askVal=FALSE, p_dissoc, p_hini)
#===============================================================================
{
    # Determines [H+] from Total alkalinity and dissolved total elements in sea water. 
    # Universal and robust algorithm from Munhoven (2013) with Newton-Raphson iterations
    # that converges from any given [H+] initial value.

    # constant threshold
    pz_exp_threshold = 1.0
    # Maximum number of iterations
    jp_maxniter_atgen    = 50
    # stopping criterion = tolerance on relative change of [H+]
    pp_rdel_ah_target = 1.E-8

    #==============================================================================


    if (missing (p_dissoc)) {
    
    } else {
    
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
        
    }
    
    if (missing (p_hini)) {

        if (exists("DEBUG_PHSOLVERS") && DEBUG_PHSOLVERS)
            print ('[solve_pH_from_AT] Calling HINI_FOR_AT for h_ini')

        zh_ini = hini_for_at(p_alktot, p_dictot, p_bortot, api)

        if (exists("DEBUG_PHSOLVERS") && DEBUG_PHSOLVERS)
            print (c('[solve_pH_from_AT] h_ini :', zh_ini))

    } else {

        zh_ini = p_hini

    }

    result = nw_infsup (p_dictot, p_bortot,
                        p_po4tot, p_siltot,  p_nh4tot, p_h2stot,
                        p_so4tot, p_flutot)
    zalknw_inf = result[1]
    zalknw_sup = result[2]

    zdelta = (p_alktot-zalknw_inf)^2 + 4.*api$api1_wat/api$aphscale

    if (p_alktot >= zalknw_inf) {
        zh_min = 2.*api$api1_wat /( p_alktot-zalknw_inf + sqrt(zdelta) )
    } else {
        zh_min = api$aphscale*(-(p_alktot-zalknw_inf) + sqrt(zdelta) ) / 2.
    }


    zdelta = (p_alktot-zalknw_sup)^2 + 4.*api$api1_wat/api$aphscale

    if (p_alktot <= zalknw_sup) {
        zh_max = api$aphscale*(-(p_alktot-zalknw_sup) + sqrt(zdelta) ) / 2.
    } else {
        zh_max = 2.*api$api1_wat /( p_alktot-zalknw_sup + sqrt(zdelta) )
    }

    if (exists("DEBUG_PHSOLVERS") && DEBUG_PHSOLVERS) {
        print (c('[solve_pH_from_AT] h_min :', zh_min))
        print (c('[solve_pH_from_AT] h_max :', zh_max))
    }
    
    zh = max(min(zh_max, zh_ini), zh_min)
    #zh = sqrt(zh_max*zh_min)              # Uncomment this line for the
                                          # "safe" initialisation test

    niter_atgen        = 0                 # Reset counters of iterations

    zeqn_absmin        = .Machine$double.xmax    #Biggest floating point number


    repeat
    {

        if (niter_atgen >= jp_maxniter_atgen) {
            zh = -1.
            break
        }

        zh_prev = zh

        result = equation_at (p_alktot, zh,   p_dictot, p_bortot,
                            p_po4tot, p_siltot, p_nh4tot, p_h2stot,
                            p_so4tot, p_flutot, api, TRUE)
        zeqn = result[1]
        zdeqndh = result[2]
        
        # Adapt bracketing interval
        if (zeqn > 0.) {
            zh_min = zh_prev
        } else if (zeqn < 0.) {
            zh_max = zh_prev
        } else {
            # zh is the root; unlikely but, one never knows
            break
        }


        # Now determine the next iterate zh
        niter_atgen = niter_atgen + 1



        if (abs(zeqn) >= 0.5*zeqn_absmin) {

            # if the function evaluation at the current point is
            # not decreasing faster than with a bisection step (at least linearly)
            # in absolute value take one bisection step on [ph_min, ph_max]
            # ph_new = (ph_min + ph_max)/2d0
            # In terms of [H]_new:
            # [H]_new = 10^(-ph_new)
            #         = 10^(-(ph_min + ph_max)/2d0)
            #         = sqrt(10^(-(ph_min + phmax)))
            #         = sqrt(zh_max * zh_min)
          
            zh = sqrt(zh_max * zh_min)

            zh_lnfactor = (zh - zh_prev)/zh_prev # Required to test convergence below

        } else {

            # dzeqn/dpH = dzeqn/d[H] * d[H]/dpH
            #           = -zdeqndh * LOG(10) * [H]
            # \Delta pH = -zeqn/(zdeqndh*d[H]/dpH) = zeqn/(zdeqndh*[H]*LOG(10))

            # pH_new = pH_old + \deltapH

            # [H]_new = 10^(-pH_new)
            #         = 10^(-pH_old - \Delta pH)
            #         = [H]_old * 10^(-zeqn/(zdeqndh*[H]_old*LOG(10)))
            #         = [H]_old * EXP(-LOG(10)*zeqn/(zdeqndh*[H]_old*LOG(10)))
            #         = [H]_old * EXP(-zeqn/(zdeqndh*[H]_old))

            zh_lnfactor = -zeqn/(zdeqndh*zh_prev)

            if (abs(zh_lnfactor) > pz_exp_threshold) {
                zh          = zh_prev*exp(zh_lnfactor)
            } else {
                zh_delta    = zh_lnfactor*zh_prev
                zh          = zh_prev + zh_delta
            }

            if (exists("DEBUG_PHSOLVERS") && DEBUG_PHSOLVERS)
                print (c('[solve_pH_from_AT] testing zh :', zh, zeqn, zh_lnfactor))


            if ( zh < zh_min ) {
                # if [H]_new < [H]_min
                # i.e., if ph_new > ph_max then
                # take one bisection step on [ph_prev, ph_max]
                # ph_new = (ph_prev + ph_max)/2d0
                # In terms of [H]_new:
                # [H]_new = 10^(-ph_new)
                #         = 10^(-(ph_prev + ph_max)/2d0)
                #         = sqrt(10^(-(ph_prev + phmax)))
                #         = sqrt([H]_old*10^(-ph_max))
                #         = sqrt([H]_old * zh_min)

                zh                = sqrt(zh_prev * zh_min)

                zh_lnfactor       = (zh - zh_prev)/zh_prev # Required to test convergence below
            }

            if ( zh > zh_max ) {
                # if [H]_new > [H]_max
                # i.e., if ph_new < ph_min, then
                # take one bisection step on [ph_min, ph_prev]
                # ph_new = (ph_prev + ph_min)/2d0
                # In terms of [H]_new:
                # [H]_new = 10^(-ph_new)
                #         = 10^(-(ph_prev + ph_min)/2d0)
                #         = sqrt(10^(-(ph_prev + ph_min)))
                #         = sqrt([H]_old*10^(-ph_min))
                #         = sqrt([H]_old * zhmax)

                zh                = sqrt(zh_prev * zh_max)

                zh_lnfactor       = (zh - zh_prev)/zh_prev # Required to test convergence below
            }

        }

        zeqn_absmin = min( abs(zeqn), zeqn_absmin)


        # Stop iterations once |\delta{[H]}/[H]| < rdel
        # <=> |(zh - zh_prev)/zh_prev| = |EXP(-zeqn/(zdeqndh*zh_prev)) -1| < rdel
        # |EXP(-zeqn/(zdeqndh*zh_prev)) -1| ~ |zeqn/(zdeqndh*zh_prev)|

        # Alternatively:
        # |\Delta pH| = |zeqn/(zdeqndh*zh_prev*LOG(10))|
        #             ~ 1/LOG(10) * |\Delta [H]|/[H]
        #             < 1/LOG(10) * rdel

        # Hence |zeqn/(zdeqndh*zh)| < rdel

        # rdel <-- pp_rdel_ah_target

        l_exitnow = (abs(zh_lnfactor) < pp_rdel_ah_target)

        if (l_exitnow) break

    }   # end repeat



    if (p_askVal) {

        if (zh > 0.) {

            val = equation_at (p_alktot, zh,       p_dictot, p_bortot,
                                p_po4tot, p_siltot, p_nh4tot, p_h2stot,
                                p_so4tot, p_flutot, api, FALSE)

        } else {
            val = .Machine$double.xmax   #Biggest floating point number
        }

        return (data.frame(zh, val))
    }
    else
    {
        return (zh)
    }

}

