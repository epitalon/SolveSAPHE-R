#===============================================================================
nw_infsup <- function (p_dictot, p_bortot,
                       p_po4tot, p_siltot,  p_nh4tot, p_h2stot,
                       p_so4tot, p_flutot)
#===============================================================================
{
    # Subroutine returns the lower and upper bounds of "non-water-selfionization"
    # contributions to total alkalinity (the infimum and the supremum), i.e
    # inf(TA - [OH-] + [H+]) and sup(TA - [OH-] + [H+])


    # alknw_inf = -\Sum_i m_i Xtot_i

    # alknw_inf =-p_dictot*0._wp &          # n = 2, m = 0
    #              -p_bortot*0._wp &          # n = 1, m = 0
    #              -p_po4tot*1._wp &          # n = 3, m = 1
    #              -p_siltot*0._wp &          # n = 1, m = 0
    #              -p_nh4tot*0._wp &          # n = 1, m = 0
    #              -p_h2stot*0._wp &          # n = 1, m = 0
    #              -p_so4tot*1._wp &          # n = 1, m = 1
    #              -p_flutot*1._wp            # n = 1, m = 1

    alknw_inf =    -p_po4tot - p_so4tot - p_flutot


    # alknw_sup = \Sum_i (n_i - m_i) Xtot_i

    # alknw_sup = p_dictot*(2._wp-0._wp) &  # n = 2, m = 0
    #               p_bortot*(1._wp-0._wp) &  # n = 1, m = 0
    #               p_po4tot*(3._wp-1._wp) &  # n = 3, m = 1
    #               p_siltot*(1._wp-0._wp) &  # n = 1, m = 0
    #               p_nh4tot*(1._wp-0._wp) &  # n = 1, m = 0
    #               p_h2stot*(1._wp-0._wp) &  # n = 1, m = 0
    #               p_so4tot*(1._wp-1._wp) &  # n = 1, m = 1
    #               p_flutot*(1._wp-1._wp)    # n = 1, m = 1

    alknw_sup =   p_dictot + p_dictot + p_bortot +
                  p_po4tot + p_po4tot + p_siltot +
                  p_nh4tot + p_h2stot

    return (c(alknw_inf, alknw_sup))
}


