# objective types, used to trigger different code paths in
# tree search and probabilistic inference
abstract type AbstractObjective end
abstract type AbstractTreeSearchSettings end

struct LeastSquares <: AbstractObjective end
Base.string(::LeastSquares) = "LS"

struct KullbackLeibler <: AbstractObjective end
Base.string(::KullbackLeibler) = "KL"
