#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// Linear-time quadratic forms and log-determinant for a phylogenetic covariance
// matrix plus a diagonal (Felsenstein 1973; Ho & Ane 2014, "three-point"
// algorithm).
//
// The covariance is  V = C(len) + diag(tip_extra),  where C(len) is the
// Brownian-motion covariance implied by the edge lengths `len` (root edge
// ignored, as in ape::vcv) and tip_extra[i] is added to the variance of tip i.
// For y and the vector of ones, it returns
//   c(log|V|, 1'V^-1 1, 1'V^-1 y, y'V^-1 y).
//
// Edges must be in postorder (every edge below a node comes before the edge
// above it), e.g. from ape::reorder.phylo(phy, "postorder"). parent and child
// use ape's numbering: tips 1..ntip, root ntip + 1, internal nodes after that.
//
// Each subtree carries the four summaries for its tips with V restricted to the
// subtree. Adding the edge of length t above a node adds t * 11' to that block;
// Sherman-Morrison and the matrix determinant lemma update the summaries in O(1)
// per edge. Siblings are independent given their parent, so their summaries add.
//
// Returns NA for all four if a tip's total variance is not positive (V is then
// singular or not a covariance matrix).

// [[Rcpp::export(.three_point)]]
NumericVector three_point(const IntegerVector& parent, const IntegerVector& child,
                          const NumericVector& len, const NumericVector& tip_extra,
                          const NumericVector& y, const int ntip, const int nnode) {
    const int nedge = parent.size();
    if (child.size() != nedge || len.size() != nedge) {
        stop("parent, child and len must have the same length.");
    }
    if (tip_extra.size() != ntip || y.size() != ntip) {
        stop("tip_extra and y must have one value per tip.");
    }
    const int nall = ntip + nnode;
    std::vector<double> logd(nall, 0.0), qq(nall, 0.0), qy(nall, 0.0), yy(nall, 0.0);

    for (int e = 0; e < nedge; e++) {
        const int p = parent[e] - 1;
        const int c = child[e] - 1;
        const double t = len[e];
        double ld, a, b, d;
        if (c < ntip) {
            const double v = t + tip_extra[c];
            if (!(v > 0.0)) {
                return NumericVector::create(NA_REAL, NA_REAL, NA_REAL, NA_REAL);
            }
            ld = std::log(v);
            a = 1.0 / v;
            b = y[c] / v;
            d = y[c] * y[c] / v;
        } else {
            const double den = 1.0 + t * qq[c];
            if (!(den > 0.0)) {
                return NumericVector::create(NA_REAL, NA_REAL, NA_REAL, NA_REAL);
            }
            ld = logd[c] + std::log(den);
            a = qq[c] / den;
            b = qy[c] / den;
            d = yy[c] - t * qy[c] * qy[c] / den;
        }
        logd[p] += ld;
        qq[p] += a;
        qy[p] += b;
        yy[p] += d;
    }

    const int root = ntip;  // node ntip + 1 in ape numbering
    return NumericVector::create(logd[root], qq[root], qy[root], yy[root]);
}
