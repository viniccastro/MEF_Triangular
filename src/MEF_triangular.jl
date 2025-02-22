module MEF_triangular

# bibliotecas
using BenchmarkTools, SparseArrays, Plots, LinearAlgebra

# arquivos
include("erro.jl")
include("base.jl")
include("galerkin.jl")
include("pontos_gauss.jl")
include("gerador_malha.jl")
include("mudanca_variavel.jl")
include("vetor_local_global.jl")
include("matriz_local_global.jl")

# funcoes
export galerkin                                 # galerkin.jl
export calculo_erro                             # erro.jl
export gauss_triangular                         # pontos_gauss.jl
export vetor_local!, vetor_global               # vetor_local_global.jl
export funcao_g, derivada_funcao_g              # mudanca_variavel.jl
export matriz_local!, matriz_global             # matriz_local_global.jl
export phi, derivada_phi, gradiente_phi         # base.jl
export malha_triangular, malha_quadricular      # gerador_malha.jl

end # module MEF_triangular