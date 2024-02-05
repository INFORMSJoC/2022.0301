#include <cassert>
#include "util.h"

/* struct to store information on GUBs */
struct GUBinfo {
    double gubweight; /**< total weight of the GUB */
    int    gubindex;  /**< index of GUB in GUB list */
};

/* compares two GUBs based on their weight
 *
 * Return -1, 0, 1 if first GUB has larger, same, smaller weight than the second
 * GUB.
 */
int compare_GUBpair(const void* a, /**< first GUB */
                    const void* b  /**< second GUB */
) {
    const struct GUBinfo* da = (const struct GUBinfo*)a;
    const struct GUBinfo* db = (const struct GUBinfo*)b;

    return da->gubweight > db->gubweight ? -1 : da->gubweight < db->gubweight;
}

/* Separates GUB facet inequalities for a knapsack polytope with knapsack
 * inequality 2 * (x_1 + ... + x_r) + (y_1 + ... + y_s) - (z_1 + ... + z_t) <=
 * rhs, where some of the y- and z-variables are contained in GUB constraints.
 *
 * To avoid handling y- and z-variables that are not contained in a GUB, we
 * assume that also single variables can form a GUB.
 *
 * We encode GUBs in the following way:
 *
 * - posvals contains the values of variables with a +1 coefficient in the
 * knapsack inequality
 * - posgubs is a list containing encoding the GUBs
 * - posgubsbegins encodes where the GUBs in posgubs begin and end
 *
 * for example,
 *
 * if posgubs = [0,3,2,4,1] and posgubsbegins = [0,2,5], then posgubs contains
 * two GUBs, the first one is formed by the entries posgubsbegins[0], ...,
 * posgubsbegins[1] - 1 in posgubs, i.e., {0,3} and the second one corresponds
 * to the entries between posgubsbegins[1], ..., posgubsbegins[2] - 1, i.e.,
 * {2,4,1}.
 *
 * The weights of the GUBs are then posvals[posgubs[posgubsbegins[0]]] + ... +
 * posvals[posgubs[posgubsbegins[1] - 1]] and posvals[posgubs[posgubsbegins[1]]]
 * + ... + posvals[posgubs[posgubsbegins[2] - 1]].
 */
void separateGUBknapsackFacets(
    double* twovals,    /**< values of solution to separate for variables with
                           2-coefficient */
    int     ntwovals,   /**< number of variables with 2-coefficient */
    double* posvals,    /**< values of solution to separate for variables with
                           1-coefficient */
    int  nposvals,      /**< number of variables with 1-coefficient */
    int* posgubs,       /**< array containing indices of the GUBs in posvals */
    int* posgubsbegins, /**< array containing the start positions of GUBs in
                           posgubs */
    int     nposgubs,   /**< number of GUBs encoded in posgubs */
    double* negvals,    /**< values of solution to separate for variables with
                           (-1)-coefficient */
    int  nnegvals,      /**< number of variables with (-1)-coefficient */
    int* neggubs,       /**< array containing indices of the GUBs in negvals */
    int* neggubsbegins, /**< array containing the start positions of GUBs in
                           neggubs */
    int   nneggubs,     /**< number of GUBs encoded in neggubs */
    int   rhs,          /**< right-hand side of knapsack inequality */
    int** twocoeffs,    /**< pointer to coefficients of 2-variables in separated
                           inequality */
    int** poscoeffs,    /**< pointer to coefficients of 1-variables in separated
                           inequality */
    int** negcoeffs, /**< pointer to coefficients of (-1)-variables in separated
                        inequality */
    int*  separhs,   /**< pointer to right-hand side of separated inequality */
    bool* success /**< pointer to store whether an inequality could be separated
                   */
) {
    int             i;
    int             j;
    int             r;
    int             index;
    double          sumtwovals = 0.0;
    double*         posgubweights;
    double*         neggubweights;
    int*            posgublabels;
    int*            neggublabels;
    int             rhsmod;
    int             nallgubs;
    double          facetlhs;
    double*         valsallgubs;
    int*            alllabels;
    struct GUBinfo* posgubweightinfo;
    struct GUBinfo* neggubweightinfo;
    struct GUBinfo* allgubweightinfo;
    int             cnt;
    int             cnt1;
    int             cnt2;

    assert(twovals != NULL);
    assert(ntwovals >= 0);
    assert(posvals != NULL);
    assert(nposvals >= 0);
    assert(posgubs != NULL);
    assert(posgubsbegins != NULL);
    assert(0 <= nposgubs && nposgubs <= nposvals);
    assert(negvals != NULL);
    assert(nnegvals >= 0);
    assert(neggubs != NULL);
    assert(neggubsbegins != NULL);
    assert(0 <= nneggubs && nneggubs <= nnegvals);
    assert(rhs >= 0);
    assert(twocoeffs != NULL);
    assert(poscoeffs != NULL);
    assert(negcoeffs != NULL);
    assert(separhs != NULL);
    assert(success != NULL);

    *success = false;

    /* store sum of all values corresponding to 2-coefficients
     * (all these variables will have a 1-coefficients in an inequality and thus
     * contribute)
     */
    for (i = 0; i < ntwovals; ++i)
        sumtwovals += twovals[i];

    /* for each GUB in posgubs, store its total weight and index (i for the i-th
     * GUB) */
    posgubweightinfo =
        (struct GUBinfo*)malloc(nposgubs * sizeof(struct GUBinfo));
    for (i = 0; i < nposgubs; ++i) {
        posgubweightinfo[i].gubweight = 0.0;
        posgubweightinfo[i].gubindex = i;

        for (j = posgubsbegins[i]; j < posgubsbegins[i + 1]; ++j)
            posgubweightinfo[i].gubweight += posvals[posgubs[j]];
    }

    /* for each GUB in neggubs, store its total negated weight and index ( - (i
     * + 1) for the i-th GUB) */
    neggubweightinfo =
        (struct GUBinfo*)malloc(nneggubs * sizeof(struct GUBinfo));
    for (i = 0; i < nneggubs; ++i) {
        neggubweightinfo[i].gubweight = 1.0;
        neggubweightinfo[i].gubindex = -(i + 1);

        for (j = neggubsbegins[i]; j < neggubsbegins[i + 1]; ++j)
            neggubweightinfo[i].gubweight -= negvals[neggubs[j]];
    }

    /* adapt the right-hand side of the knapsack inequality by negating
     * (-1)-coefficient GUBs */
    rhsmod = rhs + nneggubs;
    nallgubs = nposgubs + nneggubs;

    /* sort GUBs in non-increasing order w.r.t. their weights in a common array
     */
    qsort(posgubweightinfo, nposgubs, sizeof(struct GUBinfo), compare_GUBpair);
    qsort(neggubweightinfo, nneggubs, sizeof(struct GUBinfo), compare_GUBpair);

    /* store sorted GUBs in common array in non-increasing order */
    allgubweightinfo =
        (struct GUBinfo*)malloc(nallgubs * sizeof(struct GUBinfo));
    cnt = 0;
    cnt1 = 0;
    cnt2 = 0;
    for (cnt = 0; cnt < nallgubs; ++cnt) {
        if (cnt1 == nposgubs || posgubweightinfo[cnt1].gubweight <
                                    neggubweightinfo[cnt2].gubweight) {
            allgubweightinfo[cnt].gubweight = neggubweightinfo[cnt2].gubweight;
            allgubweightinfo[cnt].gubindex = neggubweightinfo[cnt2++].gubindex;
        } else {
            allgubweightinfo[cnt].gubweight = posgubweightinfo[cnt1].gubweight;
            allgubweightinfo[cnt].gubindex = posgubweightinfo[cnt1++].gubindex;
        }
    }

    /* iterate over all possible number of (transformed) GUBs with 1-coefficient
     */
    facetlhs = sumtwovals;
    cnt = 0;
    for (r = 0; r < nallgubs && r < rhs; ++r) {
        int facetrhs;

        /* skip dominated inequalities */
        if (r % 2 == rhsmod % 2)
            continue;

        facetrhs = r + (rhsmod - r - 1) / 2;

        if (r > 0) {
            facetlhs += allgubweightinfo[cnt].gubweight;
            ++cnt;
        }

        /* we have found a violated inequality */
        if (facetlhs > facetrhs) {
            *success = true;

            *twocoeffs = (int*)malloc(ntwovals * sizeof(int));
            *poscoeffs = (int*)malloc(nposvals * sizeof(int));
            *negcoeffs = (int*)malloc(nnegvals * sizeof(int));

            /* set coefficients of 2-variables */
            for (i = 0; i < ntwovals; ++i)
                (*twocoeffs)[i] = 1;

            /* init coefficient vectors */
            for (i = 0; i < nposvals; ++i)
                (*poscoeffs)[i] = 0;
            for (i = 0; i < nnegvals; ++i)
                (*negcoeffs)[i] = 0;

            /* set coefficients of (+/- 1)-variables */
            cnt2 = 0;
            for (i = 0; i <= r; ++i) {
                index = allgubweightinfo[i].gubindex;
                if (index >= 0) {
                    for (j = posgubsbegins[index]; j < posgubsbegins[index + 1];
                         ++j)
                        (*poscoeffs)[posgubs[j]] = 1;
                } else {
                    ++cnt2;
                    index = allgubweightinfo[i].gubindex;
                    index = -(index + 1);

                    for (j = neggubsbegins[index]; j < neggubsbegins[index + 1];
                         ++j)
                        (*negcoeffs)[neggubs[j]] = -1;
                }
            }

            *separhs = facetrhs - cnt2;
            break;
        }

        facetlhs += allgubweightinfo[cnt].gubweight;
        ++cnt;
    }

    if (!*success)
        printf("No violated inequality found.\n");

    free(allgubweightinfo);
    free(neggubweightinfo);
    free(posgubweightinfo);
}
