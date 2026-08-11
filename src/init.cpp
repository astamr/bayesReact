#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" SEXP fast_kmers_counting(SEXP, SEXP);
extern "C" SEXP pwm_prob_dp(SEXP, SEXP, SEXP, SEXP, SEXP);

// register native functions that R may call via .Call()
static const R_CallMethodDef CallEntries[] = {
  {"fast_kmers_counting",
   (DL_FUNC) &fast_kmers_counting, 2},
  {"pwm_prob_dp",
   (DL_FUNC) &pwm_prob_dp, 5},
  {NULL, NULL, 0}
};

// called when the bayesReact package DLL is loaded
extern "C" void R_init_bayesReact(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
