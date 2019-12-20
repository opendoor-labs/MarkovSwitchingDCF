// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Rginv
arma::mat Rginv(arma::mat m);
RcppExport SEXP _MarkovSwitchingDCF_Rginv(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(Rginv(m));
    return rcpp_result_gen;
END_RCPP
}
// gen_inv
arma::mat gen_inv(arma::mat m);
RcppExport SEXP _MarkovSwitchingDCF_gen_inv(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_inv(m));
    return rcpp_result_gen;
END_RCPP
}
// ss_prob
arma::mat ss_prob(arma::mat mat);
RcppExport SEXP _MarkovSwitchingDCF_ss_prob(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(ss_prob(mat));
    return rcpp_result_gen;
END_RCPP
}
// kalman_filter
Rcpp::List kalman_filter(arma::mat B0, arma::mat P0, arma::mat Dt, arma::mat At, arma::mat Ft, arma::mat Ht, arma::mat Qt, arma::mat Rt, arma::mat yt);
RcppExport SEXP _MarkovSwitchingDCF_kalman_filter(SEXP B0SEXP, SEXP P0SEXP, SEXP DtSEXP, SEXP AtSEXP, SEXP FtSEXP, SEXP HtSEXP, SEXP QtSEXP, SEXP RtSEXP, SEXP ytSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P0(P0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dt(DtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type At(AtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ft(FtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ht(HtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Qt(QtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Rt(RtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type yt(ytSEXP);
    rcpp_result_gen = Rcpp::wrap(kalman_filter(B0, P0, Dt, At, Ft, Ht, Qt, Rt, yt));
    return rcpp_result_gen;
END_RCPP
}
// kalman_smoother
Rcpp::List kalman_smoother(arma::mat B_tl, arma::mat B_tt, arma::cube P_tl, arma::cube P_tt, arma::mat Ft);
RcppExport SEXP _MarkovSwitchingDCF_kalman_smoother(SEXP B_tlSEXP, SEXP B_ttSEXP, SEXP P_tlSEXP, SEXP P_ttSEXP, SEXP FtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type B_tl(B_tlSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B_tt(B_ttSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type P_tl(P_tlSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type P_tt(P_ttSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ft(FtSEXP);
    rcpp_result_gen = Rcpp::wrap(kalman_smoother(B_tl, B_tt, P_tl, P_tt, Ft));
    return rcpp_result_gen;
END_RCPP
}
// kim_filter
Rcpp::List kim_filter(arma::cube B0, arma::cube P0, arma::cube At, arma::cube Dt, arma::cube Ft, arma::cube Ht, arma::cube Qt, arma::cube Rt, arma::mat Tr_mat, arma::mat yt, bool weighted);
RcppExport SEXP _MarkovSwitchingDCF_kim_filter(SEXP B0SEXP, SEXP P0SEXP, SEXP AtSEXP, SEXP DtSEXP, SEXP FtSEXP, SEXP HtSEXP, SEXP QtSEXP, SEXP RtSEXP, SEXP Tr_matSEXP, SEXP ytSEXP, SEXP weightedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type P0(P0SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type At(AtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Dt(DtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Ft(FtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Ht(HtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Qt(QtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Rt(RtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Tr_mat(Tr_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type yt(ytSEXP);
    Rcpp::traits::input_parameter< bool >::type weighted(weightedSEXP);
    rcpp_result_gen = Rcpp::wrap(kim_filter(B0, P0, At, Dt, Ft, Ht, Qt, Rt, Tr_mat, yt, weighted));
    return rcpp_result_gen;
END_RCPP
}
// kim_smoother
Rcpp::List kim_smoother(arma::cube B_tlss, arma::cube B_tts, arma::mat B_tt, arma::field<arma::cube> P_tlss, arma::field<arma::cube> P_tts, arma::mat Pr_tls, arma::mat Pr_tts, arma::cube At, arma::cube Dt, arma::cube Ft, arma::cube Ht, arma::cube Qt, arma::cube Rt, arma::mat Tr_mat);
RcppExport SEXP _MarkovSwitchingDCF_kim_smoother(SEXP B_tlssSEXP, SEXP B_ttsSEXP, SEXP B_ttSEXP, SEXP P_tlssSEXP, SEXP P_ttsSEXP, SEXP Pr_tlsSEXP, SEXP Pr_ttsSEXP, SEXP AtSEXP, SEXP DtSEXP, SEXP FtSEXP, SEXP HtSEXP, SEXP QtSEXP, SEXP RtSEXP, SEXP Tr_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type B_tlss(B_tlssSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type B_tts(B_ttsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B_tt(B_ttSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::cube> >::type P_tlss(P_tlssSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::cube> >::type P_tts(P_ttsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Pr_tls(Pr_tlsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Pr_tts(Pr_ttsSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type At(AtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Dt(DtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Ft(FtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Ht(HtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Qt(QtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Rt(RtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Tr_mat(Tr_matSEXP);
    rcpp_result_gen = Rcpp::wrap(kim_smoother(B_tlss, B_tts, B_tt, P_tlss, P_tts, Pr_tls, Pr_tts, At, Dt, Ft, Ht, Qt, Rt, Tr_mat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MarkovSwitchingDCF_Rginv", (DL_FUNC) &_MarkovSwitchingDCF_Rginv, 1},
    {"_MarkovSwitchingDCF_gen_inv", (DL_FUNC) &_MarkovSwitchingDCF_gen_inv, 1},
    {"_MarkovSwitchingDCF_ss_prob", (DL_FUNC) &_MarkovSwitchingDCF_ss_prob, 1},
    {"_MarkovSwitchingDCF_kalman_filter", (DL_FUNC) &_MarkovSwitchingDCF_kalman_filter, 9},
    {"_MarkovSwitchingDCF_kalman_smoother", (DL_FUNC) &_MarkovSwitchingDCF_kalman_smoother, 5},
    {"_MarkovSwitchingDCF_kim_filter", (DL_FUNC) &_MarkovSwitchingDCF_kim_filter, 11},
    {"_MarkovSwitchingDCF_kim_smoother", (DL_FUNC) &_MarkovSwitchingDCF_kim_smoother, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_MarkovSwitchingDCF(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
