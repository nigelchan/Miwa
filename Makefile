RESULT = main
SOURCES = \
  utils/utils.ml \
  \
  math/mw_math.ml \
  math/mw_stat.ml \
  math/mw_matrix.ml \
  math/mw_integral.ml \
  math/mw_optimize.ml \
  \
  model/model.ml \
  \
  process/process.ml \
  \
  engine/engine.ml \
  engine/payoff.ml \
  engine/analytic/analyticHeston.ml \
  engine/montecarlo/mcTools.ml \
  engine/montecarlo/mcEngine.ml \
  \
  model/calibratedModel.ml \
  main.ml

INCDIRS = math/gsl
LIBS = Bigarray Gsl

OCAMLMAKEFILE = OCamlMakefile
include $(OCAMLMAKEFILE)
