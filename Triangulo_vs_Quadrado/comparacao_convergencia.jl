include("Elementos Triangulares.jl")
include("Elementos Quadriculares.jl")

# Usar os módulos
using .tri
using .quad
using Plots
using Profile

function problema() :: Tuple
    alpha         = 1.0
    beta          = 1.0
    funcao_f(x) = (2 * π^2 + 1) * sin(π*x[1]) * sin(π*x[2])
    funcao_u(x) = sin(π*x[1]) * sin(π*x[2])       # funcao exata

    return (alpha, beta, funcao_f, funcao_u)

end

function plot_erro_vs(max_h::Vector{Float64}, erros1::Vector{Float64}, erros2::Vector{Float64})
    max_h2 = (x -> x^2).(max_h)

    # grafico
    plot(max_h, max_h2, label="O(h²)", legs =:topleft, seriestype=:path, linestyle=:dash, xaxis=:log, yaxis=:log,)
    plot!(max_h, erros1, label="Erro_trian ", marker=:utriangle)
    plot!(max_h, erros2, label="Erro_quadr ", marker=:square)
    xlabel!("h")
    ylabel!("Erros")
end

gauss_tri = tri.gauss_triangulo(7)
gauss_quad = quad.gauss_quadrado(5)

# maximo = 6
numero_casos = 6
lista_max_h = zeros(numero_casos)
erros_tri = zeros(numero_casos)
erros_quad = zeros(numero_casos)

for i in 2:numero_casos+1
    println("caso ", i, "/", numero_casos+1)

    malha_tri = tri.gerador_malha([2^i, 2^i])
    erros_tri[i-1] = tri.erro( tri.galerkin(problema(), malha_tri) , malha_tri, gauss_tri, problema()[end])
    
    # malha_quad = quad.gerador_malha([2^i, 2^i])
    # erros_quad[i-1] = quad.erro( quad.galerkin(problema(), malha_quad) , malha_quad, gauss_quad, problema()[end])
    
    lista_max_h[i-1] = malha_tri[1]
end
tri.plot_erro(lista_max_h, erros_tri)


# graficos da aproximacao
# elem = 2^7
# malha_tri = tri.gerador_malha([elem, elem])
# vetor_c_tri = tri.galerkin(problema(), malha_tri)
# tri.plot_grafico(1, vetor_c_tri, malha_tri[3], malha_tri[end-2], malha_tri[end])

# malha_quad = quad.gerador_malha([elem, elem])
# vetor_c_quad = quad.galerkin(problema(), malha_quad)
# quad.plot_grafico(1, vetor_c_quad, malha_quad[3], malha_quad[end-2], malha_quad[end])
# Profile.print()
# plot_erro_vs(lista_max_h, erros_tri, erros_quad)

# function plot_grafico(tri_vetor::Vector{Float64}, tri_lg::Matrix{Int}, quad_vetor::Vector{Float64}, quad_lg::Matrix{Int}, pontos::Matrix{Float64}, eq::Vector{Int})
#     tri_numero_p = (Int(sqrt(length(tri_vetor)))+2)^2
#     tri_vetor = [tri_vetor; 0]
#     tri_altura = zeros(tri_numero_p)
#     for i in 1:(tri_numero_p)
#         tri_altura[i] = tri_vetor[eq[i]]
#     end
#     tri_numero_e = size(tri_lg)[2]

#     plt = plot(seriestype=:surface, layout=(1,3), xticks=[0, 1], yticks=[0, 1], zticks=[0, 1], color=:viridis, xlims=(0,1), ylims=(0,1), zlims=(0,1), colorbar=false, legend=false,camera = (30, 30, 1))
#     xlabel!("x")
#     ylabel!("y")
#     for e in 1:tri_numero_e
#         lista_pontos = pontos[:,tri_lg[:,e]]
#         lista_altura = tri_altura[tri_lg[:,e]]
#         lista_pontos = hcat(lista_pontos, lista_pontos[:,end])
#         lista_altura = vcat(lista_altura, lista_altura[end])
#         coordenadas_x = transpose(reshape(lista_pontos[1,:], 2, 2))
#         coordenadas_y = transpose(reshape(lista_pontos[2,:], 2, 2))
#         coordenadas_z = transpose(reshape(lista_altura, 2, 2))
        
#         plot!(coordenadas_x, coordenadas_y, coordenadas_z, color=:viridis, grid=1)
#     end


#     quad_numero_p = (Int(sqrt(length(quad_vetor)))+2)^2
#     quad_vetor = [quad_vetor; 0]
#     quad_altura = zeros(quad_numero_p)
#     for i in 1:(quad_numero_p)
#         quad_altura[i] = quad_vetor[eq[i]]
#     end
    
#     quad_lg[[3, 4], :] .= quad_lg[[4, 3], :]
#     quad_numero_e = size(quad_lg)[2]
    
#     xlabel!("x")
#     ylabel!("y")
#     for e in 1:quad_numero_e
#         lista_pontos = pontos[:,quad_lg[:,e]]
#         lista_altura = quad_altura[quad_lg[:,e]]
#         coordenadas_x = transpose(reshape(lista_pontos[1,:], 2, 2))
#         coordenadas_y = transpose(reshape(lista_pontos[2,:], 2, 2))
#         coordenadas_z = transpose(reshape(lista_altura, 2, 2))
        
#         plot!(coordenadas_x, coordenadas_y, coordenadas_z, color=:viridis, grid=2)
    
#         # coordenadas_x = reshape(lista_pontos[1,:], 2, 2)
#         # coordenadas_y = reshape(lista_pontos[2,:], 2, 2)
#         # coordenadas_z = reshape(lista_altura, 2, 2)
        
#         # plot!(coordenadas_x, coordenadas_y, coordenadas_z, color=:viridis, layout=(1,2))
#     end

#     xs = collect(0:0.01:1)
#     ys = collect(0:0.01:1)
#     X = [x for x = xs for _ = ys]
#     Y = [y for _ = xs for y = ys]
#     Z = ((x,y)->begin
#         sin(π*x) * sin(π*y)
#         end)
#     plot!(X, Y, Z.(X, Y),color=:viridis, grid=(1,3))

#     display(plt)
# end
# plot_grafico(vetor_c_tri, malha_tri[end-2], vetor_c_quad, malha_quad[end-2], malha_tri[3], malha_tri[end])