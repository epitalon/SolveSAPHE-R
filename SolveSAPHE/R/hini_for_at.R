#===============================================================================
 hini_for_at <- function(p_alkcb, p_dictot, p_bortot, api)
#===============================================================================
{
    # Subroutine returns the root for the 2nd order approximation of the
    # DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic polynomial)
    # around the local minimum, if it exists.

    # Returns * 1E-03 if p_alkcb <= 0
    #         * 1E-10 if p_alkcb >= 2*p_dictot + p_bortot
    #         * 1E-07 if 0 < p_alkcb < 2*p_dictot + p_bortot
    #                    and the 2nd order approximation does not have a solution



    if (p_alkcb <= 0.) {
        hini = 1.e-3
    } else if (p_alkcb >= (2.*p_dictot + p_bortot)) {
        hini = 1.e-10
    } else {
        zca = p_dictot/p_alkcb
        zba = p_bortot/p_alkcb

        # Coefficients of the cubic polynomial
        za2 = api$api1_bor*(1. - zba) + api$api1_dic*(1.-zca)
        za1 = api$api1_dic*api$api1_bor*(1. - zba - zca) + api$api2_dic*(1. - (zca+zca))
        za0 = api$api2_dic*api$api1_bor*(1. - zba - (zca+zca))


        # Taylor expansion around the minimum

        zd = za2*za2 - 3.*za1              # Discriminant of the quadratic equation
                                              # for the minimum close to the root

        if (zd > 0.) {                   # If the discriminant is positive

            zsqrtd = sqrt(zd)

            if (za2 < 0) {
                zhmin = (-za2 + zsqrtd)/3.
            } else {
                zhmin = -za1/(za2 + zsqrtd)
            }

            hini = zhmin + sqrt(-(za0 + zhmin*(za1 + zhmin*(za2 + zhmin)))/zsqrtd)

        } else {

            hini = 1.e-7
        }
    }

    return (hini)
}


