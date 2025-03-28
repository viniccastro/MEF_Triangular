module tri # begin module

using Plots
using LinearAlgebra
using SparseArrays

export gauss_triangulo, gerador_malha, galerkin, erro, plot_grafico

function problema() :: Tuple
    alpha         = 1.0
    beta          = 1.0
    funcao_f(x) = (2 * π^2 + 1) * sin(π*x[1]) * sin(π*x[2])
    funcao_u(x) = sin(π*x[1]) * sin(π*x[2])       # funcao exata

    return (alpha, beta, funcao_f, funcao_u)

end

function gauss_triangulo(numero_p::Int) :: Tuple{Int, Vector{Vector{Float64}}, Vector{Float64}}
    
    if numero_p == 1

        # coordenadas dos pontos
        pontos = [[1/3, 1/3]]
    
        # pesos dos pontos
        pesos = [1/2]

    elseif numero_p == 3

        # coordenadas dos pontos
        pontos = [[1/6, 1/6], [2/3, 1/6], [1/6, 2/3]]
    
        # pesos dos pontos
        pesos = [1/6, 1/6, 1/6]

    elseif numero_p == 4

        # coordenadas dos pontos
        pontos = [[1/3, 1/3], [1/5, 1/5], [3/5, 1/5], [1/5, 3/5]]
    
        # pesos dos pontos
        pesos = [-27/96, 25/96, 25/96, 25/96]

    elseif numero_p == 6

        # coordenadas dos pontos
        valor1 = 0.445948490915965
        valor2 = 0.091576213509771
        pontos = [[valor1, valor1], [1-2*valor1, valor1], [valor1, 1-2*valor1], [valor2, valor2], [1-2*valor2, valor2], [valor2, 1-2*valor2]]
    
        # pesos dos pontos
        valor3 = 0.111690794839005
        valor4 = 0.054975871827661
        pesos = [valor3, valor3, valor3, valor4, valor4, valor4]

    elseif numero_p == 7

        # coordenadas dos pontos
        valor1 = 0.470142064105115
        valor2 = 0.101286507323456
        pontos = [[1/3, 1/3], [valor1, valor1], [1-2*valor1, valor1], [valor1, 1-2*valor1], [valor2, valor2], [1-2*valor2, valor2], [valor2, 1-2*valor2]]
    
        # pesos dos pontos
        valor3 = 0.0661970763942530
        valor4 = 0.0629695902724135
        pesos = [9/80, valor3, valor3, valor3, valor4, valor4, valor4]

    else
        
        error("Parametro 'numero_p' invalido, informe entre 1, 3, 4, 6 ou 7")

    end

    # retorna o numero de pontos, as coordenadas e os pesos de cada ponto
    return numero_p, pontos, pesos

end

function phi(phi_n::Int, ponto_xi::Vector{Float64}) :: Float64
    if phi_n == 1
        return 1 - ponto_xi[1] - ponto_xi[2]

    elseif phi_n == 2
        return ponto_xi[1]

    elseif phi_n == 3
        return ponto_xi[2]

    else
        error("Parametro ``phi_n`` invalido, argumento esperado é um inteiro de 1 a 3")

    end
end

function derivada_phi(dimensao::Int, phi_n::Int, ponto_xi::Vector{Float64}) :: Float64
    @assert dimensao == 1 || dimensao == 2 "Parametro ``dimensão`` inválido, argumento esperado é 1 ou 2"

    if phi_n == 1
        return -1

    elseif phi_n == 2
        return 2-dimensao

    elseif phi_n == 3
        return dimensao-1

    else
        error("Parametro ``phi_n`` invalido, argumento esperado é um inteiro de 1 a 3")

    end
end

function funcao_g(ponto_xi::Vector{Float64}, pontos_e::Matrix{Float64}) :: Vector{Float64}
    return phi(1, ponto_xi)*pontos_e[:, 1] + phi(2, ponto_xi)*pontos_e[:, 2] + phi(3, ponto_xi)*pontos_e[:, 3]

end

function gradiente_g(ponto_xi::Vector{Float64}, pontos_e::Matrix{Float64}) :: Tuple{Float64, Float64, Float64, Float64}
    grad1 = derivada_phi(1, 1, ponto_xi)*pontos_e[1, 1] + derivada_phi(1, 2, ponto_xi)*pontos_e[1, 2] + derivada_phi(1, 3, ponto_xi)*pontos_e[1, 3]
    grad2 = derivada_phi(1, 1, ponto_xi)*pontos_e[2, 1] + derivada_phi(1, 2, ponto_xi)*pontos_e[2, 2] + derivada_phi(1, 3, ponto_xi)*pontos_e[2, 3]
    grad3 = derivada_phi(2, 1, ponto_xi)*pontos_e[1, 1] + derivada_phi(2, 2, ponto_xi)*pontos_e[1, 2] + derivada_phi(2, 3, ponto_xi)*pontos_e[1, 3]
    grad4 = derivada_phi(2, 1, ponto_xi)*pontos_e[2, 1] + derivada_phi(2, 2, ponto_xi)*pontos_e[2, 2] + derivada_phi(2, 3, ponto_xi)*pontos_e[2, 3]
    
    return grad1, grad2, grad3, grad4
end

function vetor_local!(vector_local!::Vector{Float64}, pontos_e::Matrix{Float64}, gauss::Tuple, funcao_f::Function)
    numero_gauss = gauss[1] :: Int
    pontos_gauss = gauss[2] :: Vector{Vector{Float64}}
    pesos_gauss  = gauss[3] :: Vector{Float64}
    
    vector_local! .= zeros(3)
    
    for i = 1:numero_gauss
        peso = pesos_gauss[i]
        ponto_xi = pontos_gauss[i]

        valor_f = funcao_f(funcao_g(ponto_xi, pontos_e))
        
        jacobiano1, jacobiano2, jacobiano3, jacobiano4 = gradiente_g(ponto_xi, pontos_e)
        determinante_j = jacobiano1*jacobiano4 - jacobiano2*jacobiano3
        @assert determinante_j > 0 "Determinante do jacobiano deve ser positivo"
        
        for a = 1:3
            vector_local![a] += peso * valor_f * phi(a, ponto_xi) * determinante_j
            
        end
    end
end

function matriz_local!(matrix_local!::Matrix{Float64}, pontos_e::Matrix{Float64}, gauss::Tuple, alpha::Float64, beta::Float64)
    numero_gauss = gauss[1] :: Int
    pontos_gauss = gauss[2] :: Vector{Vector{Float64}}
    pesos_gauss  = gauss[3] :: Vector{Float64}
    
    matrix_local! .= zeros(3, 3)

    for i in 1:numero_gauss
        peso = pesos_gauss[i]
        ponto_xi = pontos_gauss[i]

        jacobiano1, jacobiano2, jacobiano3, jacobiano4 = gradiente_g(ponto_xi, pontos_e)
        determinante_j = jacobiano1*jacobiano4 - jacobiano2*jacobiano3
        @assert determinante_j > 0 "Determinante do jacobiano deve ser positivo"

        hh11 =  jacobiano3^2 + jacobiano4^2
        hh21 = -jacobiano1*jacobiano3 - jacobiano2*jacobiano4
        hh22 =  jacobiano1^2 + jacobiano2^2

        for b in 1:3
            phi_b = phi(b, ponto_xi)
            derivada_phi_b1 = derivada_phi(1, b, ponto_xi)
            derivada_phi_b2 = derivada_phi(2, b, ponto_xi)

            for a in 1:3
                derivada_phi_a1 = derivada_phi(1, a, ponto_xi)
                derivada_phi_a2 = derivada_phi(2, a, ponto_xi)

                integral_alpha = alpha * (derivada_phi_b1*(derivada_phi_a1*hh11 + derivada_phi_a2*hh21) + derivada_phi_b2*(derivada_phi_a1*hh21 + derivada_phi_a2*hh22))
                integral_beta  = beta * phi(a, ponto_xi) * phi_b

                matrix_local![a, b] += peso * ((1/determinante_j) * integral_alpha + integral_beta * determinante_j)

            end
        end
    end
end

function monta_global(dados_padrao::Tuple, gauss::Tuple, dados_problema::Tuple)
    alpha    = dados_problema[1] :: Float64
    beta     = dados_problema[2] :: Float64
    funcao_f = dados_problema[3] :: Function

    numero_e = dados_padrao[2] :: Int
    pontos   = dados_padrao[3] :: Matrix{Float64}
    valor_m  = dados_padrao[4] :: Int
    lg       = dados_padrao[5] :: Matrix{Int}
    eqlg     = dados_padrao[6] :: Function
    
    vetor_local  = zeros(3)
    matriz_local = zeros(3,3)

    vetor_global  = zeros(valor_m+1)
    matriz_global = sparse(zeros(valor_m+1, valor_m+1))
   
    for e in 1:numero_e
        pontos_e = pontos[:, lg[:,e]]
        vetor_local!(vetor_local, pontos_e, gauss, funcao_f)
        matriz_local!(matriz_local, pontos_e, gauss, alpha, beta)
        
        for b in 1:3
            j = eqlg(b, e)
            vetor_global[j] += vetor_local[b]

            for a in 1:3
                i = eqlg(a, e)
                matriz_global[i,j] += matriz_local[a,b]
                
            end
        end
    end
    return (matriz_global[1:valor_m, 1:valor_m], vetor_global[1:valor_m])

end

function gerador_malha(numero_e::Vector{Int}) :: Tuple{Float64, Int, Matrix{Float64}, Int, Matrix{Int}, Function, Vector{Int}}

    # verificador
    @assert numero_e[1] > 1 "valor informado deve ser maior iqual a 2"

    # tamaho dos passos espaciais
    valor_h = 1 / numero_e[1]
    max_h   = valor_h*sqrt(2)

    # coordenadas dos pontos da malha
    pontos = transpose([0:1:numero_e[1]  fill(0, numero_e[1]+1)])
    for i in 1:(numero_e[1])
        
        pontos = hcat(pontos, transpose([0:1:numero_e[1]  fill(i, numero_e[1]+1)]))

    end

    # valor do m
    valor_m = (numero_e[1]-1)^2

    # matriz LG
    matriz_lg = zeros(Int,3,0)
    for l in 1:numero_e[1]

        padrao1 = (x -> div(x,2)).(0:(2*numero_e[1]-1))
        padrao2 = (x -> div(x,2)).(1:(2*numero_e[1]))
        padrao3 = (x -> mod(x,2)).(0:(2*numero_e[1]-1))

        linha1 = padrao2 .+ (1 + (l-1)*(numero_e[1]+1))
        linha3 = padrao1 .+ (1 + l*(numero_e[1]+1))

        sequencia1 = padrao1 .+ (linha1[1]+1)
        sequencia2 = padrao1 .+ (linha3[1]+1)
        linha2 = sequencia1 .*(1 .- padrao3) .+ sequencia2 .* padrao3

        matriz_lg = hcat(matriz_lg, transpose([linha1 linha2 linha3]))

    end

    # vetor EQ
    vetor_eq  = zeros(Int,0)
    for i in 1:(numero_e[1]-1)
        
        vetor_eq = [vetor_eq; [valor_m+1 ; (1:(numero_e[1]-1)) .+ ((i-1)*(numero_e[1]-1)) ; valor_m+1]]
        
    end
    tampa = ones(Int,numero_e[1]+1) .* (valor_m+1)
    vetor_eq = [tampa ; vetor_eq ; tampa]

    EQLG(p,e) = vetor_eq[matriz_lg[p,e]]

    # retorna o vetor de equacao
    return max_h, 2*prod(numero_e), valor_h*pontos, valor_m, matriz_lg, EQLG, vetor_eq

end

function galerkin(dados_problema::Tuple, dados_malha::Tuple) :: Vector{Float64}
    matriz_K, vetor_F = monta_global(dados_malha, gauss_triangulo(7), dados_problema)
    return matriz_K \ vetor_F

end

function erro(vetor_C::Vector{Float64}, dados_padrao::Tuple, gauss::Tuple, f_exata::Function) :: Float64
    # dados padroes
    numero_e = dados_padrao[2] :: Int
    pontos   = dados_padrao[3] :: Matrix{Float64}
    lg       = dados_padrao[5] :: Matrix{Int}
    eqlg     = dados_padrao[6] :: Function

    # valores de Gauss
    numero_gauss = gauss[1] :: Int
    pontos_gauss = gauss[2] :: Vector{Vector{Float64}}
    pesos_gauss  = gauss[3] :: Vector{Float64}
    
    # adicao do valor dos extremos
    vetor_C = cat(vetor_C, 0, dims=1)

    # inicializacao do erro
    erro = 0

    for i in 1:numero_gauss
        peso = pesos_gauss[i]
        ponto_xi = pontos_gauss[i]
        
        for e in 1:numero_e
            pontos_e = pontos[:, lg[:,e]]

            jacobiano1, jacobiano2, jacobiano3, jacobiano4 = gradiente_g(ponto_xi, pontos_e)
            determinante_j = jacobiano1*jacobiano4 - jacobiano2*jacobiano3
            @assert determinante_j > 0 "Determinante do jacobiano deve ser positivo"
            
            valor_u = f_exata( funcao_g(ponto_xi, pontos_e))
            valor_e = 0

            for a in 1:3
                valor_e += vetor_C[ eqlg(a, e) ] * phi(a, ponto_xi)

            end
            erro += determinante_j * peso * (valor_u - valor_e)^2

        end
    end
    # retorna a raiz quadrada do erro
    return sqrt(erro)

end

function plot_erro(max_h::Vector{Float64}, erros::Vector{Float64})
    max_h2 = (x -> x^2).(max_h)

    plot(max_h, max_h2, label="O(h²)", legs =:topleft, seriestype=:path, linestyle=:dash, xaxis=:log, yaxis=:log,)
    plot!(max_h, erros, label="Erro_trian ", marker=:utriangle)
    xlabel!("h")
    ylabel!("Erros")
end

function plot_grafico(versao::Int, vetor::Vector{Float64}, pontos::Matrix{Float64}, lg::Matrix{Int}, eq::Vector{Int})
    if versao == 1
        numero_p = (Int(sqrt(length(vetor)))+2)^2
        vetor = [vetor; 0]
        altura = zeros(numero_p)
        for i in 1:(numero_p)
            altura[i] = vetor[eq[i]]
        end
    elseif versao == 2
        altura = vetor
    else
        println("Parametro ``versao`` invalido, informe 1 ou 2")
    end
    
    numero_e = size(lg)[2]

    plt = plot(seriestype=:surface, xticks=[0, 1], yticks=[0, 1], zticks=[0, 1], color=:viridis, xlims=(0,1), ylims=(0,1), zlims=(0,1), colorbar=false, legend=false,camera = (30, 30, 1), size=(400, 400))
    xlabel!("x")
    ylabel!("y")
    for e in 1:numero_e
        lista_pontos = pontos[:,lg[:,e]]
        lista_altura = altura[lg[:,e]]
        lista_pontos = hcat(lista_pontos, lista_pontos[:,end])
        lista_altura = vcat(lista_altura, lista_altura[end])
        coordenadas_x = transpose(reshape(lista_pontos[1,:], 2, 2))
        coordenadas_y = transpose(reshape(lista_pontos[2,:], 2, 2))
        coordenadas_z = transpose(reshape(lista_altura, 2, 2))
        
        plot!(coordenadas_x, coordenadas_y, coordenadas_z, color=:viridis)
    end
    display(plt)
end

# numero_casos = 5
# lista_max_h = zeros(numero_casos)
# lista_erros = zeros(numero_casos)

# for i in 1:numero_casos
#     malha = gerador_malha([2^i, 2^i])
#     lista_erros[i] = erro( galerkin(problema(), malha) , malha, gauss_triangulo(7), problema()[end])
#     lista_max_h[i] = malha[1]

# end

# plot_erro(lista_max_h, lista_erros)

end # end modulo
