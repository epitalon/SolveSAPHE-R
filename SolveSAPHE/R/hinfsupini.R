
# General parameters

# Threshold relative difference between successive iterates
# (convergence criterion)
pp_rdel_ah_target = 1.E-08

# NaN for [H^+] results
pp_hnan = -1.

#===============================================================================
HINFSUPINI <- function (p_alktot,   p_dicvar, p_bortot,
                                    p_po4tot, p_siltot,
                                    p_nh4tot, p_h2stot,
                                    p_so4tot, p_flutot,
                                    p_dicsel, api,  p_hini)
#===============================================================================

{

    #--------------------#
    # Argument variables #
    #--------------------#

    #  p_alktot       :  Alkalinity total
    #  p_dicvar       :  Carbonate value
    #  p_bortot       :  Dissolved bore total
    #  p_po4tot       :  Dissolved Phosphate total
    #  p_siltot       :  Dissolved Silicium total
    #  p_nh4tot       :  Dissolved Ammonium total
    #  p_h2stot       :  Dissolved Hydrogene sulfide total
    #  p_so4tot       :  Dissolved Sulfate total
    #  p_flutot       :  Dissolved Fluor total
    #  p_dicsel       :  Carbonate variable selector (default = DIC)
    #                      p_dicsel = "DIC":   p_dicvar = DIC
    #                      p_dicsel = "CO2":   p_dicvar = [CO2*]
    #                      p_dicsel = "HCO3":  p_dicvar = [HCO3-]
    #                      p_dicsel = "CO3":   p_dicvar = [CO3--]
    #  api            :  Pi factors of dissociation constants
    #  p_hini         :  Initial value of pH (optionnal), 1 or 2 element long array

    #-------------------------------------------------------------------------------


    result = NULL

    if (p_dicsel == "DIC")
    {
        result = HINFSUPINI_DIC(p_alktot, p_dicvar, p_bortot,
                                    p_po4tot, p_siltot,
                                    p_nh4tot, p_h2stot,
                                    p_so4tot, p_flutot,
                                    api, p_hini)

    }
    else if (p_dicsel == "CO2")
    {
        result = HINFSUPINI_CO2(p_alktot, p_dicvar, p_bortot,
                                    p_po4tot, p_siltot,
                                    p_nh4tot, p_h2stot,
                                    p_so4tot, p_flutot,
                                    api, p_hini)

    }
    else if (p_dicsel == "HCO3")
    {
        result = HINFSUPINI_HCO3(p_alktot, p_dicvar, p_bortot,
                                    p_po4tot, p_siltot,
                                    p_nh4tot, p_h2stot,
                                    p_so4tot, p_flutot,
                                    api, p_hini)

    }
    else if (p_dicsel == "CO3")
    {
        result = HINFSUPINI_CO3(p_alktot, p_dicvar, p_bortot,
                                    p_po4tot, p_siltot,
                                    p_nh4tot, p_h2stot,
                                    p_so4tot, p_flutot,
                                    api, p_hini)

    }
    
    if (is.null(result))
    {
        k_nroots = 0

        p_hinf = pp_hnan
        p_hsup = pp_hnan
        p_hini = pp_hnan
        
        result = list(p_hinf=p_hinf, p_hsup=p_hsup, p_hini=p_hini, k_nroots=k_nroots)
    }

    return (result)
}

#===============================================================================
HINFSUPINI_DIC <- function (p_alktot, p_dictot, p_bortot,
                            p_po4tot, p_siltot,
                            p_nh4tot, p_h2stot,
                            p_so4tot, p_flutot,
                            api,      p_hini)
#===============================================================================

# Subroutine provides
#  - the infimum and the supremum of Alk_{nW}, obtained from ANW_INFSUP
#  - a valid initial value for any solver, already brought within brackets.
# If the given p_hini is not equal to pp_hnan, the calculation of h_ini from
# the ACB approximation is skipped and the provided value is only brought into
# [p_hinf, p_hsup] if necessary.
{

    #--------------------#
    # Argument variables #
    #--------------------#

    #  p_alktot       :  Alkalinity total
    #  p_dictot       :  DIC total
    #  p_bortot       :  Dissolved bore total
    #  p_po4tot       :  Dissolved Phosphate total
    #  p_siltot       :  Dissolved Silicium total
    #  p_nh4tot       :  Dissolved Ammonium total
    #  p_h2stot       :  Dissolved Hydrogene sulfide total
    #  p_so4tot       :  Dissolved Sulfate total
    #  p_flutot       :  Dissolved Fluor total
    #  api            :  Pi factors of dissociation constants
    #  p_hini         :  Initial value of pH (optionnal)

    #-------------------------------------------------------------------------------



    k_nroots  = 1


    result = anw_infsup(p_dictot, p_bortot,
                        p_po4tot, p_siltot,
                        p_nh4tot, p_h2stot,
                        p_so4tot, p_flutot)
    zalknw_inf = result[1]
    zalknw_sup = result[2]
    
    za1 = p_alktot - zalknw_inf
    zd  = za1^2 + 4.*api$api1_wat/api$aphscale

    if (za1 > 0) {
        z_hinf = 2.*api$api1_wat / ( za1 + sqrt(zd) )
    } else {
        z_hinf =         api$aphscale*(-za1 + sqrt(zd) ) / 2.
    }


    za1 = p_alktot - zalknw_sup
    zd  = za1^2 + 4.*api$api1_wat/api$aphscale

    if (za1 > 0) {
        z_hsup = 2.*api$api1_wat / ( za1 + sqrt(zd) )
    } else {
        z_hsup =       api$aphscale * (-za1 + sqrt(zd) ) / 2.
    }


    if (exists("SAFEGEOMEAN_INIT") && SAFEGEOMEAN_INIT)
    {
        if (   (p_hini[1] == pp_hnan)          # request to calculate H_ini
            || (p_hini[1] < z_hinf)          # the given H_ini is too low
            || (p_hini[1] > z_hsup) )        # the given H_ini too high
        {
            p_hini[1] = sqrt(z_hsup*z_hinf)
        }

    }
    else
    {
        if (p_hini[1] == pp_hnan) {      # request to calculate H_ini

            if (exists("DEBUG_PHSOLVERS") && DEBUG_PHSOLVERS)
            {
                print ('[HINFSUPINI] using hini_ACB_DIC to set H_ini')
            }

            p_hini[1] = hini_ACB_DIC(p_alktot, p_dictot, p_bortot, z_hinf, z_hsup, api)

        } 
        else 
        {    
            # H_ini given: only bring it into [z_hinf, z_hsup]
            p_hini[1] = max(min(z_hsup, p_hini[1]), z_hinf)
        }

    }


    p_hinf = c(z_hinf, pp_hnan)
    p_hsup = c(z_hsup, pp_hnan)
    p_hini[2] = pp_hnan

    return ( list(p_hinf=p_hinf, p_hsup=p_hsup, p_hini=p_hini, k_nroots=k_nroots))
}


#===============================================================================
HINFSUPINI_CO2 <- function (p_alktot, p_co2,   p_bortot,
                            p_po4tot, p_siltot,
                            p_nh4tot, p_h2stot,
                            p_so4tot, p_flutot,
                            api,      p_hini)
#===============================================================================

# Subroutine provides
#  - the infimum and the supremum of Alk_{nW}, obtained from ANW_INFSUP
#  - a valid initial value for any solver, already brought within brackets.
# If the given p_hini is not equal to pp_hnan, the calculation of h_ini from
# the ACB approximation is skipped and the provided value is only brought into
# [p_hinf, p_hsup] if necessary.
{

    #--------------------#
    # Argument variables #
    #--------------------#

    #  p_alktot       :  Alkalinity total
    #  p_co2          :  [CO2*] concentration
    #  p_bortot       :  Dissolved bore total
    #  p_po4tot       :  Dissolved Phosphate total
    #  p_siltot       :  Dissolved Silicium total
    #  p_nh4tot       :  Dissolved Ammonium total
    #  p_h2stot       :  Dissolved Hydrogene sulfide total
    #  p_so4tot       :  Dissolved Sulfate total
    #  p_flutot       :  Dissolved Fluor total
    #  api            :  Pi factors of dissociation constants
    #  p_hini         :  Initial value of pH (optionnal)

    #-------------------------------------------------------------------------------



    k_nroots  = 1


    # Get Alk_nWCinf and Alk_nWCsup:
    z_dictot = 0.                    # set C_T to 0 an call ANW_INFSUP

    result = anw_infsup(z_dictot, p_bortot,
                        p_po4tot, p_siltot,
                        p_nh4tot, p_h2stot,
                        p_so4tot, p_flutot)
    zalknwc_inf = result[1]
    zalknwc_sup = result[2]

    # Lower root bracket
    # would normally require the solution of a cubic equation.  It is 
    # nevertheless sufficient to chose the location of the minimum of the
    # cubic (written such that a_3 > 0).

                                        # Coefficients of cubic polynomial
    za3 = 1./api$aphscale
    za2 = (p_alktot - zalknwc_inf)
    za1 = -(api$api1_dic * p_co2 + api$api1_wat)
    # za0 = -2. * api$api2_dic * p_dicvar

    # The derivative has as constant term za1 < 0 
    # => one positive and one negative root.
    # Determine the positive one,
    # i.e., the greater one of the two
    zd   = za2^2 - 3.*za1*za3

    if (za2 > 0) {
        z_hinf = -za1 / ( za2 + sqrt(zd) )
    } else {
        z_hinf =        (-za2 + sqrt(zd) ) / (3.*za3)
    }


    # Upper root bracket:
    # would normally require the solution of a cubic equation.
    # It is nevertheless sufficient to chose the greater of the roots of
    # the quadratic expansion around the minimum of the cubic (written
    # such that a_3 > 0)

    # Coefficients of cubic polynomial
    za3 = 1./api$aphscale
    za2 = (p_alktot - zalknwc_sup)
    za1 = -(api$api1_dic*p_co2 + api$api1_wat)
    za0 = -2. * api$api2_dic * p_co2

    # The derivative has as constant term za1 < 0 
    # => one positive and one negative root.
    # Determine the positive one,
    # i.e., the greater one of the two
    zd     = za2^2 - 3.*za1*za3
    zsqrtd = sqrt(zd)

    if (za2 > 0) {
        z_hcmin =     -za1 / ( za2 + zsqrtd )
    } else {
        z_hcmin =            (-za2 + zsqrtd ) / (3.*za3)
    }

    z_hsup = z_hcmin + sqrt(-(za0 + z_hcmin*(za1 + z_hcmin*(za2 + z_hcmin*za3)))/zsqrtd)


    if (exists("SAFEGEOMEAN_INIT") && SAFEGEOMEAN_INIT)
    {

        if (   (p_hini[1] == pp_hnan)           # request to calculate H_ini
            || (p_hini[1] < z_hinf)             # the given H_ini is too low
            || (p_hini[1] > z_hsup) )           # the given H_ini too high
        {
            p_hini[1] = sqrt(z_hsup*z_hinf)
        }
    }
    else
    {
        if (p_hini[1] == pp_hnan) {      # request to calculate H_ini

            if (exists("DEBUG_PHSOLVERS") && DEBUG_PHSOLVERS)
            {
                print ('[HINFSUPINI] using hini_ACB_CO2 to set H_ini')
            }

            p_hini[1] = hini_ACB_CO2(p_alktot, p_co2, p_bortot, z_hinf, z_hsup, api)

        } 
        else 
        {
            # H_ini given: only bring it into [z_hinf, z_hsup]
            p_hini[1] = max(min(z_hsup, p_hini[1]), z_hinf)

        }
    }


    p_hinf = c(z_hinf, pp_hnan)
    p_hsup = c(z_hsup, pp_hnan)
    p_hini[2] = pp_hnan

    return ( list(p_hinf=p_hinf, p_hsup=p_hsup, p_hini=p_hini, k_nroots=k_nroots))

}



#===============================================================================
HINFSUPINI_HCO3 <- function (p_alktot, p_hco3, p_bortot,
                            p_po4tot, p_siltot,
                            p_nh4tot, p_h2stot,
                            p_so4tot, p_flutot,
                            api,      p_hini)
#===============================================================================

# Subroutine provides
#  - the infimum and the supremum of Alk_{nW}, obtained from ANW_INFSUP
#  - a valid initial value for any solver, already brought within brackets.
# If the given p_hini is not equal to pp_hnan, the calculation of h_ini from
# the ACB approximation is skipped and the provided value is only brought into
# [p_hinf, p_hsup] if necessary.
{

    #--------------------#
    # Argument variables #
    #--------------------#

    #  p_alktot       :  Alkalinity total
    #  p_hco3         :  [HCO3-] concentration
    #  p_bortot       :  Dissolved bore total
    #  p_po4tot       :  Dissolved Phosphate total
    #  p_siltot       :  Dissolved Silicium total
    #  p_nh4tot       :  Dissolved Ammonium total
    #  p_h2stot       :  Dissolved Hydrogene sulfide total
    #  p_so4tot       :  Dissolved Sulfate total
    #  p_flutot       :  Dissolved Fluor total
    #  api            :  Pi factors of dissociation constants
    #  p_hini         :  Initial value of pH (optionnal)

    #-------------------------------------------------------------------------------

    k_nroots  = 1


    # Get Alk_nWCinf and Alk_nWCsup:
    z_dictot = 0.                    # set C_T to 0 an call ANW_INFSUP

    result = anw_infsup(z_dictot, p_bortot,
                        p_po4tot, p_siltot,
                        p_nh4tot, p_h2stot,
                        p_so4tot, p_flutot)
    zalknwc_inf = result[1]
    zalknwc_sup = result[2]

    # za2 =  1./api$aphscale
    za1 =  p_alktot - zalknwc_inf - p_hco3
    za0 = -(2.*(api$api2_dic/api$api1_dic)*p_hco3 + api$api1_wat)

    zd  = za1^2 - 4. * za0 / api$aphscale

    if (za1 > 0) {
        z_hinf = -2. * za0 / ( za1 + sqrt(zd) )
    } else {
        z_hinf =       api$aphscale*(-za1 + sqrt(zd) ) / 2.
    }


    # za2 =  1./api$aphscale
    za1 =  p_alktot - zalknwc_sup - p_hco3
    # za0 = -(2.*(api$api2_dic/api$api1_dic)*p_hco3 + api$api1_wat)

    zd    = za1^2 - 4. * za0 / api$aphscale

    if (za1 > 0) {
        z_hsup = -2. * za0 / ( za1 + sqrt(zd) )
    } else {
        z_hsup =     api$aphscale * (-za1 + sqrt(zd) ) / 2.
    }

    if (exists("SAFEGEOMEAN_INIT") && SAFEGEOMEAN_INIT)
    {

        if (   (p_hini[1] == pp_hnan)           # request to calculate H_ini
            || (p_hini[1] < z_hinf)             # the given H_ini is too low
            || (p_hini[1] > z_hsup) )           # the given H_ini too high
        {
            p_hini[1] = sqrt(z_hsup*z_hinf)
        }
    }
    else
    {
        if (p_hini[1] == pp_hnan) {      # request to calculate H_ini

            if (exists("DEBUG_PHSOLVERS") && DEBUG_PHSOLVERS)
            {
                print ('[HINFSUPINI] using hini_ACBW_HCO3 to set H_ini')
            }

            p_hini[1] = hini_ACBW_HCO3(p_alktot, p_hco3, p_bortot, z_hinf, z_hsup, api)

        } 
        else 
        {
            # H_ini given: only bring it into [z_hinf, z_hsup]
            p_hini[1] = max(min(z_hsup, p_hini[1]), z_hinf)

        }
    }


    p_hinf = c(z_hinf, pp_hnan)
    p_hsup = c(z_hsup, pp_hnan)
    p_hini[2] = pp_hnan

    return ( list(p_hinf=p_hinf, p_hsup=p_hsup, p_hini=p_hini, k_nroots=k_nroots))

}


#===============================================================================
HINFSUPINI_CO3 <- function (p_alktot, p_co3,   p_bortot,
                            p_po4tot, p_siltot,
                            p_nh4tot, p_h2stot,
                            p_so4tot, p_flutot,
                            api,      p_hini)
#===============================================================================

# Subroutine provides
#  - the infimum and the supremum of Alk_{nW}, obtained from ANW_INFSUP
#  - a valid initial value for any solver, already brought within brackets.
# If the given p_hini is not equal to pp_hnan, the calculation of h_ini from
# the ACB approximation is skipped and the provided value is only brought into
# [infimum, supremum] if necessary.


{
    #--------------------#
    # Argument variables #
    #--------------------#

    #  p_alktot       :  Alkalinity total
    #  p_co3          :  [CO3--] concentration
    #  p_bortot       :  Dissolved bore total
    #  p_po4tot       :  Dissolved Phosphate total
    #  p_siltot       :  Dissolved Silicium total
    #  p_nh4tot       :  Dissolved Ammonium total
    #  p_h2stot       :  Dissolved Hydrogene sulfide total
    #  p_so4tot       :  Dissolved Sulfate total
    #  p_flutot       :  Dissolved Fluor total
    #  api            :  Pi factors of dissociation constants
    #  p_hini         :  Initial value of pH (optionnal), 2 element long array

    #-------------------------------------------------------------------------------



    z_dictot = 0.

    result = anw_infsup(z_dictot, p_bortot,
                        p_po4tot, p_siltot,
                        p_nh4tot, p_h2stot,
                        p_so4tot, p_flutot, api)
    zalknwc_inf = result[1]
    zalknwc_sup = result[2]
    zalknwc_asympt_coeff = result[3]
    

    z_gamma = p_co3 * (api$api1_dic/api$api2_dic)- 1./api$aphscale


    if (z_gamma < 0.) { #-----------------------------------------------------
                        #-----------------------------------------------------
                        #-----------------------------------------------------
                        
        k_nroots  = 1

                                            # H_inf[1] = abscissa of P_UL,
                                            # which is the positive root of
                                            # the following quadratic
        za2 = z_gamma
        za1 = - (p_alktot - (p_co3 + p_co3) - zalknwc_inf)
        za0 = api$api1_wat


        zd  = za1^2 - 4. * api$api1_wat * z_gamma

        if (za1 > 0) {
            z_hinf1 =               -( za1 + sqrt(zd) ) / (za2 + za2)
        } else {
            z_hinf1 = -(za0 + za0) / ( za1 - sqrt(zd) )
        }


                                            # H_sup[1] = abscissa of P_LL,
                                            # which is the positive root of
                                            # the following quadratic
        za2 = z_gamma
        za1 = - (p_alktot - (p_co3 + p_co3) - zalknwc_sup)
        za0 = api$api1_wat

        zd  = za1^2 - 4. * za0 * za2

        if (za1 > 0.) {
            z_hsup1 = -(za0 + za0) / ( za1 - sqrt(zd) )
        } else {
            z_hsup1 =               -( za1 + sqrt(zd) ) / (za2 + za2)
        }

        z_hinf2 = pp_hnan
        z_hsup2 = pp_hnan

    }
    else if (z_gamma > 0.) { #-------------------------------------------------
                            #-------------------------------------------------
                            #-------------------------------------------------

        # 0, 1, or 2 roots ?

        z_hmin = sqrt(api$api1_wat/z_gamma)
        z_lmin = 2. * sqrt(z_gamma * api$api1_wat) + (p_co3 + p_co3)


        zalknwc_hmin = ANW (z_hmin,              
                            p_co3,    p_bortot, 
                            p_po4tot, p_siltot, 
                            p_nh4tot, p_h2stot, 
                            p_so4tot, p_flutot, "CO3", api)


        if (p_alktot >= (z_lmin + zalknwc_hmin)) 
        {

            k_nroots = 2

                                            # H_inf[1] = abscissa of P_UL, 
                                            # which is the lower (positive)
                                            # root of the following quadratic;
                                            # H_sup[2] = abscissa of P_UR, 
                                            # which is the greater (positive)
                                            # root of the following quadratic;
            za2 =  z_gamma
            za1 = -(p_alktot - (p_co3 + p_co3) - zalknwc_inf)
            za0 =  api$api1_wat

            zd  = za1^2 - 4. * za0 * za2


            if (za1 > 0) {
                z_hinf1 =               -( za1 + sqrt(zd) ) / (za2 + za2)
            } else {
                z_hinf1 = -(za0 + za0) / ( za1 - sqrt(zd) )
            }


            if (za1 > 0) {
                z_hsup2 = -(za0 + za0) / ( za1 + sqrt(zd) )
            } else {
                z_hsup2 =                (-za1 + sqrt(zd) ) / (za2 + za2)
            }



            if (p_alktot > (z_lmin + zalknwc_sup)) {

                                            # H_sup[1] = abscissa of P_LL, 
                                            # which is the lower (positive)
                                            # root of the following quadratic;
                                            # H_inf[2] = abscissa of P_LR, 
                                            # which is the greater (positive)
                                            # root of the following quadratic;
                za2 =  z_gamma
                za1 = -(p_alktot - (p_co3 + p_co3) - zalknwc_sup)
                za0 =  api$api1_wat

                zd  = za1^2 - 4. * za0 * za2

                if (za1 > 0) {
                        z_hsup1 =               -( za1 + sqrt(zd) ) / (za2 + za2)
                } else {
                        z_hsup1 = -(za0 + za0) / ( za1 - sqrt(zd) )
                }


                if (za1 > 0) {
                    z_hinf2 = -(za0 + za0) / ( za1 + sqrt(zd) )
                } else {
                    z_hinf2 =                (-za1 + sqrt(zd) ) / (za2 + za2)
                }


            } else {

                z_hsup1 = z_hmin
                z_hinf2 = z_hmin

            }

        }
        else if (p_alktot <= (z_lmin + zalknwc_inf)) 
        {

                k_nroots = 0

                z_hinf1  = pp_hnan
                z_hsup1  = pp_hnan

                z_hinf2  = pp_hnan
                z_hsup2  = pp_hnan


        }
        else 
        {
                                        # Alk_T is in the intermediate region
                                        # where we need to determine H_tan and Alk_tan

                                        # Brackets of the H_tan:
                                        # (H_min, A_min) and (H_UR, A_UR)

                                        # H_UR = abscissa of P_UR, 
                                        # which is the greater (positive)
                                        # root of the following quadratic:
                                        # If Alk_T > Alk_tan, H_UR will also be
                                        # H_sup[2]
            za2 =  z_gamma
            za1 = -(p_alktot - (p_co3 + p_co3) - zalknwc_inf)
            za0 =  api$api1_wat

            zd  = za1^2 - 4. * za0 * za2

            if (za1 > 0) {
                z_hsup2 = -(za0 + za0) / ( za1 + sqrt(zd) )
            } else {
                z_hsup2 =                (-za1 + sqrt(zd) ) / (za2 + za2)
            }



            z_tol = z_hmin * pp_rdel_ah_target      # We assume that H_min is
                                                    # of the same order of magnitude
                                                    # as H_tan for fixing the absolute
                                                    # tolerance on H_tan

            # z_atan = ALK_TAN(z_hmin, z_hsup2, z_tol, z_htan)
            f <- function (z_h) equation_at (0., z_h, p_co3, p_bortot, 
                      p_po4tot, p_siltot, p_nh4tot, p_h2stot, p_so4tot, p_flutot, 4, api)

            # Use Brent algorithm to find the minimum of f(z_h)
            xmin <- optimize(f, c(z_hmin, z_hsup2), tol = z_tol)
            z_htan = xmin$minimum
            z_atan = xmin$objective
     


            if (p_alktot < z_atan) {

                k_nroots = 0

                z_hinf1  = pp_hnan
                z_hsup1  = pp_hnan

                z_hinf2  = pp_hnan
                z_hsup2  = pp_hnan
        
            }
            else if (p_alktot > z_atan) 
            {

                k_nroots = 2

                                                # H_inf[1] = abscissa of P_UL, 
                                                # which is the lower (positive)
                                                # root of the following quadratic;
                                                # H_sup[2] = abscissa of P_UR, which
                                                # has already served as a bracket for H_tan
                za2 =  z_gamma
                za1 = -(p_alktot - (p_co3 + p_co3) - zalknwc_inf)
                za0 =  api$api1_wat

                zd = za1*za1 - 4. * za2 * za0

                if (za1 > 0) {
                    z_hinf1 = -(za0 + za0) / ( za1 - sqrt(zd) )
                } else {
                    z_hinf1 =               -( za1 + sqrt(zd) ) / (za2 + za2)
                }


                z_hsup1  = z_htan
                z_hinf2  = z_htan

            }
            else
            {

                k_nroots = 1

                z_hinf1  = z_htan
                z_hsup1  = z_htan

                z_hinf2  = pp_hnan
                z_hsup2  = pp_hnan

            }

        }


    } 
    else { # z_gamma = 0 -----------------------------------------------------
        #             -----------------------------------------------------
        #             -----------------------------------------------------

        zdiff = p_alktot - (p_co3 + p_co3) - zalknwc_inf

        if (zdiff > 0) {

            k_nroots = 1

            z_hinf1 = api$api1_wat / zdiff
            z_hsup1 = (api$api1_wat + zalknwc_asympt_coeff)  / zdiff

            z_hinf2  = pp_hnan
            z_hsup2  = pp_hnan

        } else {

            k_nroots = 0
            
            z_hinf1  = pp_hnan
            z_hsup1  = pp_hnan

            z_hinf2  = pp_hnan
            z_hsup2  = pp_hnan

        }

    }   #-----------------------------------------------------------------
        #-----------------------------------------------------------------
        #-----------------------------------------------------------------

    z_hinf = c(z_hinf1, z_hinf2)
    z_hsup = c(z_hsup1, z_hsup2)


    if (exists("SAFEGEOMEAN_INIT") && SAFEGEOMEAN_INIT)
    {
        if (   (p_hini[1] == pp_hnan)           # request to calculate H_ini
            || (p_hini[1] < z_hinf1)             # the given H_ini is too low
            || (p_hini[1] > z_hsup1) )           # the given H_ini too high
        {
            p_hini[1] = sqrt(z_hsup1*z_hinf1)
        }

        if (   (p_hini[2] == pp_hnan)           # request to calculate H_ini
            || (p_hini[2] < z_hinf1)             # the given H_ini is too low
            || (p_hini[2] > z_hsup1) )           # the given H_ini too high
        {
            p_hini[2] = sqrt(z_hsup2*z_hinf2)

        }

    }
    else
    {
        z_hini = hini_ACBW_CO3(p_alktot, p_co3, p_bortot, z_hinf, z_hsup, k_nroots, api)

        if (p_hini[1] == pp_hnan) # request to calculate H_ini[1]
        {
            p_hini[1] = z_hini[1]
        }
        else 
        {
            # H_ini[1] given: only bring it into [z_hinf1, z_hsup1]
            p_hini[1] = max(min(z_hsup1, p_hini[1]), z_hinf1)
        }

        if (p_hini[2] == pp_hnan)      # request to calculate H_ini[1]
        {
            p_hini[2] = z_hini[2]
        }
        else 
        {
            # H_ini[2] given: only bring it into [z_hinf2, z_hsup2]
            p_hini[2] = max(min(z_hsup2, p_hini[2]), z_hinf2)
        }

    }


    p_hinf = z_hinf
    p_hsup = z_hsup
    return ( list(p_hinf=p_hinf, p_hsup=p_hsup, p_hini=p_hini, k_nroots=k_nroots))
}

