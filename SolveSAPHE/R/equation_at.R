#===============================================================================
 equation_at <- function(p_alktot, p_h,       p_dictot, p_bortot,
                         p_po4tot, p_siltot,  p_nh4tot, p_h2stot,
                         p_so4tot, p_flutot,  api, p_deriv)
#===============================================================================
{

    # H2CO3 - HCO3 - CO3 : n=2, m=0
    znumer_dic <- 2.*api$api2_dic + p_h*       api$api1_dic
    zdenom_dic <-       api$api2_dic + p_h*(      api$api1_dic + p_h)
    zalk_dic   <- p_dictot * (znumer_dic/zdenom_dic)

    # B(OH)3 - B(OH)4 : n=1, m=0
    znumer_bor <-       api$api1_bor
    zdenom_bor <-       api$api1_bor + p_h
    zalk_bor   <- p_bortot * (znumer_bor/zdenom_bor)

    # H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
    znumer_po4 <- 3.*api$api3_po4 + p_h*(2.*api$api2_po4 + p_h* api$api1_po4)
    zdenom_po4 <-       api$api3_po4 + p_h*(      api$api2_po4 + p_h*(api$api1_po4 + p_h))
    zalk_po4   <- p_po4tot * (znumer_po4/zdenom_po4 - 1.) # Zero level of H3PO4 = 1

    # H4SiO4 - H3SiO4 : n=1, m=0
    znumer_sil <-       api$api1_sil
    zdenom_sil <-       api$api1_sil + p_h
    zalk_sil   <- p_siltot * (znumer_sil/zdenom_sil)

    # NH4 - NH3 : n=1, m=0
    znumer_nh4 <-       api$api1_nh4
    zdenom_nh4 <-       api$api1_nh4 + p_h
    zalk_nh4   <- p_nh4tot * (znumer_nh4/zdenom_nh4)

    # H2S - HS : n=1, m=0
    znumer_h2s <-       api$api1_h2s
    zdenom_h2s <-       api$api1_h2s + p_h
    zalk_h2s   <- p_h2stot * (znumer_h2s/zdenom_h2s)

    # HSO4 - SO4 : n=1, m=1
    znumer_so4 <-       api$api1_so4
    zdenom_so4 <-       api$api1_so4 + p_h
    zalk_so4   <- p_so4tot * (znumer_so4/zdenom_so4 - 1.)

    # HF - F : n=1, m=1
    znumer_flu <-       api$api1_flu
    zdenom_flu <-       api$api1_flu + p_h
    zalk_flu   <- p_flutot * (znumer_flu/zdenom_flu - 1.)

    # H2O - OH
    zalk_wat   <- api$api1_wat/p_h - p_h/api$aphscale


    zdiffalk <-   (zalk_dic + zalk_bor + zalk_po4 + zalk_sil 
                  + zalk_nh4 + zalk_h2s + zalk_so4 + zalk_flu 
                  + zalk_wat) - p_alktot

    # Do we need to compute derivative dAlk/dH ?
    if (p_deriv) 
    {

      # H2CO3 - HCO3 - CO3 : n=2
      zdnumer_dic <- api$api1_dic*api$api2_dic + p_h*(4.*api$api2_dic 
                                      + p_h*       api$api1_dic)
      zdalk_dic   <- -p_dictot*(zdnumer_dic/zdenom_dic^2)

      # B(OH)3 - B(OH)4 : n=1
      zdnumer_bor <- api$api1_bor
      zdalk_bor   <- -p_bortot*(zdnumer_bor/zdenom_bor^2)

      # H3PO4 - H2PO4 - HPO4 - PO4 : n=3
      zdnumer_po4 <- api$api2_po4*api$api3_po4 + p_h*(4.*api$api1_po4*api$api3_po4
                                      + p_h*(9.*api$api3_po4 + api$api1_po4*api$api2_po4
                                      + p_h*(4.*api$api2_po4
                                      + p_h*       api$api1_po4)))
      zdalk_po4   <- -p_po4tot * (zdnumer_po4/zdenom_po4^2)

      # H4SiO4 - H3SiO4 : n=1
      zdnumer_sil <- api$api1_sil
      zdalk_sil   <- -p_siltot * (zdnumer_sil/zdenom_sil^2)

      # NH4 - NH3 : n=1
      zdnumer_nh4 <- api$api1_nh4
      zdalk_nh4   <- -p_nh4tot * (zdnumer_nh4/zdenom_nh4^2)

      # H2S - HS : n=1
      zdnumer_h2s <- api$api1_h2s
      zdalk_h2s   <- -p_h2stot * (zdnumer_h2s/zdenom_h2s^2)

      # HSO4 - SO4 : n=1
      zdnumer_so4 <- api$api1_so4
      zdalk_so4   <- -p_so4tot * (zdnumer_so4/zdenom_so4^2)

      # HF - F : n=1
      zdnumer_flu <- api$api1_flu
      zdalk_flu   <- -p_flutot * (zdnumer_flu/zdenom_flu^2)

      zddiffalk <-   (zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil
                    + zdalk_nh4 + zdalk_h2s + zdalk_so4 + zdalk_flu
                    ) - api$api1_wat/p_h^2 - 1./api$aphscale
                    
      return (c(zdiffalk, zddiffalk))
    }
    else 
    {
      return (zdiffalk)
    }
}


