#include "ksgeneral/ksgeneral.h"

#include <math.h>
#include <stdio.h>

int main(void) {
    if (ksg_abi_version() != 1u) return 1;

    int32_t m[] = {1, 1};
    double w[] = {1.0};
    double out = NAN;
    ksg_status_t status = ksg_ks2_probability_summary(
        1, 1, KSG_ALT_TWO_SIDED, m, 2, 0.5, w, 1, 1e-8,
        KSG_PROB_LT, &out);
    if (status != KSG_OK || !isfinite(out)) {
        char err[256];
        ksg_last_error(err, sizeof(err));
        fprintf(stderr, "%s\n", err);
        return 2;
    }
    return 0;
}
