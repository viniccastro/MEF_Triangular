module quad # begin module

using Plots
using LinearAlgebra
using SparseArrays
using GaussQuadrature

export gauss_quadrado, gerador_malha, galerkin, erro, plot_grafico

function problema() :: Tuple
    alpha         = 1.0
    beta          = 1.0
    funcao_f(x) = (2 * π^2 + 1) * sin(π*x[1]) * sin(π*x[2])
    funcao_u(x) = sin(π*x[1]) * sin(π*x[2])       # funcao exata

    return (alpha, beta, funcao_f, funcao_u)

end

function gauss_quadrado(numero_p::Int) :: Tuple{Int, Vector{Float64}, Vector{Float64}}
    
    pontos, pesos = legendre(numero_p)

    # retorna o numero de pontos, as coordenadas e os pesos de cada ponto
    return numero_p, pontos, pesos

end

function phi(phi_n::Int, ponto_xi::Vector{Float64}) :: Float64
    if phi_n == 1
        return 0.25*(1 - ponto_xi[1])*(1 - ponto_xi[2])

    elseif phi_n == 2
        return 0.25*(1 + ponto_xi[1])*(1 - ponto_xi[2])

    elseif phi_n == 3
        return 0.25*(1 + ponto_xi[1])*(1 + ponto_xi[2])

    elseif phi_n == 4
        return 0.25*(1 - ponto_xi[1])*(1 + ponto_xi[2])

    else
        error("Parametro ``phi_n`` invalido, argumento esperado é um inteiro de 1 a 4")

    end
end

function derivada_phi(dimensao::Int, phi_n::Int, ponto_xi::Vector{Float64}) :: Float64
    @assert dimensao == 1 || dimensao == 2 "Parametro ``dimensão`` inválido, argumento esperado é 1 ou 2"

    if phi_n == 1
        return 0.25*(ponto_xi[3-dimensao] - 1)

    elseif phi_n == 2
        return (2-dimensao)*0.25*(1 - ponto_xi[2]) + (dimensao-1)*0.25*(-1 - ponto_xi[1])

    elseif phi_n == 3
        return 0.25*(1 + ponto_xi[3-dimensao])

    elseif phi_n == 4
        return (2-dimensao)*0.25*(-1 - ponto_xi[2]) + (dimensao-1)*0.25*(1 - ponto_xi[1])

    else
        error("Parametro ``phi_n`` invalido, argumento esperado é um inteiro de 1 a 4")

    end
end

function funcao_g(ponto_xi::Vector{Float64}, pontos_e::Matrix{Float64}) :: Vector{Float64}
    return phi(1, ponto_xi)*pontos_e[:, 1] + phi(2, ponto_xi)*pontos_e[:, 2] + phi(3, ponto_xi)*pontos_e[:, 3] + phi(4, ponto_xi)*pontos_e[:, 4]

end

function gradiente_g(ponto_xi::Vector{Float64}, pontos_e::Matrix{Float64}) :: Tuple{Float64, Float64, Float64, Float64}
    grad1 = derivada_phi(1, 1, ponto_xi)*pontos_e[1, 1] + derivada_phi(1, 2, ponto_xi)*pontos_e[1, 2] + derivada_phi(1, 3, ponto_xi)*pontos_e[1, 3] + derivada_phi(1, 4, ponto_xi)*pontos_e[1, 4]
    grad2 = derivada_phi(1, 1, ponto_xi)*pontos_e[2, 1] + derivada_phi(1, 2, ponto_xi)*pontos_e[2, 2] + derivada_phi(1, 3, ponto_xi)*pontos_e[2, 3] + derivada_phi(1, 4, ponto_xi)*pontos_e[2, 4]
    grad3 = derivada_phi(2, 1, ponto_xi)*pontos_e[1, 1] + derivada_phi(2, 2, ponto_xi)*pontos_e[1, 2] + derivada_phi(2, 3, ponto_xi)*pontos_e[1, 3] + derivada_phi(2, 4, ponto_xi)*pontos_e[1, 4]
    grad4 = derivada_phi(2, 1, ponto_xi)*pontos_e[2, 1] + derivada_phi(2, 2, ponto_xi)*pontos_e[2, 2] + derivada_phi(2, 3, ponto_xi)*pontos_e[2, 3] + derivada_phi(2, 4, ponto_xi)*pontos_e[2, 4]
    
    return grad1, grad2, grad3, grad4
end

function vetor_local!(vector_local!::Vector{Float64}, pontos_e::Matrix{Float64}, gauss::Tuple, funcao_f::Function)
    numero_gauss = gauss[1] :: Int
    pontos_gauss = gauss[2] :: Vector{Float64}
    pesos_gauss  = gauss[3] :: Vector{Float64}
    
    vector_local! .= zeros(4)
    
    for i = 1:numero_gauss
        peso1 = pesos_gauss[i]

        for j in 1:numero_gauss
            peso2 = pesos_gauss[j]
            ponto_xi = [pontos_gauss[i], pontos_gauss[j]]
    
            valor_f = funcao_f(funcao_g(ponto_xi, pontos_e))
            
            jacobiano1, jacobiano2, jacobiano3, jacobiano4 = gradiente_g(ponto_xi, pontos_e)
        determinante_j = jacobiano1*jacobiano4 - jacobiano2*jacobiano3
            @assert determinante_j > 0 "Determinante do jacobiano deve ser positivo"
            
            for a = 1:4
                vector_local![a] += peso1 * peso2 * valor_f * phi(a, ponto_xi) * determinante_j
                
            end
        end
    end
end

function matriz_local!(matrix_local!::Matrix{Float64}, pontos_e::Matrix{Float64}, gauss::Tuple, alpha::Float64, beta::Float64)
    numero_gauss = gauss[1] :: Int
    pontos_gauss = gauss[2] :: Vector{Float64}
    pesos_gauss  = gauss[3] :: Vector{Float64}
    
    matrix_local! .= zeros(4, 4)

    for i in 1:numero_gauss
        peso1 = pesos_gauss[i]

        for j in 1:numero_gauss
            peso2 = pesos_gauss[j]
            ponto_xi = [pontos_gauss[i], pontos_gauss[j]]

            jacobiano1, jacobiano2, jacobiano3, jacobiano4 = gradiente_g(ponto_xi, pontos_e)
            determinante_j = jacobiano1*jacobiano4 - jacobiano2*jacobiano3
            @assert determinante_j > 0 "Determinante do jacobiano deve ser positivo"

            hh11 =  jacobiano3^2 + jacobiano4^2
            hh21 = -jacobiano1*jacobiano3 - jacobiano2*jacobiano4
            hh22 =  jacobiano1^2 + jacobiano2^2

            for b in 1:4
                phi_b = phi(b, ponto_xi)
                derivada_phi_b1 = derivada_phi(1, b, ponto_xi)
                derivada_phi_b2 = derivada_phi(2, b, ponto_xi)

                for a in 1:4
                    derivada_phi_a1 = derivada_phi(1, a, ponto_xi)
                    derivada_phi_a2 = derivada_phi(2, a, ponto_xi)

                    integral_alpha = alpha * (derivada_phi_b1*(derivada_phi_a1*hh11 + derivada_phi_a2*hh21) + derivada_phi_b2*(derivada_phi_a1*hh21 + derivada_phi_a2*hh22))
                    integral_beta  = beta * phi(a, ponto_xi) * phi_b

                    matrix_local![a, b] += peso1 * peso2 * ((1/determinante_j) * integral_alpha + integral_beta * determinante_j)

                end
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
    
    vetor_local  = zeros(4)
    matriz_local = zeros(4,4)

    vetor_global  = zeros(valor_m+1)
    matriz_global = sparse(zeros(valor_m+1, valor_m+1))
   
    for e in 1:numero_e
        pontos_e = pontos[:, lg[:,e]]
        vetor_local!(vetor_local, pontos_e, gauss, funcao_f)
        matriz_local!(matriz_local, pontos_e, gauss, alpha, beta)
        
        for b in 1:4
            j = eqlg(b, e)
            vetor_global[j] += vetor_local[b]

            for a in 1:4
                i = eqlg(a, e)
                matriz_global[i,j] += matriz_local[a,b]

            end
        end
    end
    return (matriz_global[1:valor_m, 1:valor_m], vetor_global[1:valor_m])

end

function gerador_malha(elementos::Vector{Int}) :: Tuple{Float64, Int, Matrix{Float64}, Int, Matrix{Int}, Function, Vector{Int}}

    # tamaho dos passos espaciais
    h1 = 1 / elementos[1]
    h2 = 1 / elementos[2]
    vetor_h = [h1, h2]
    max_h   = sqrt(vetor_h[1]^2 + vetor_h[2]^2)

    # inicializando a lista de pontos com a primeira fileira de pontos
    pontos = [0:vetor_h[1]:1  fill(0, elementos[1]+1)]'
    
    # adicionando as outras fileiras
    for i in 1:(elementos[2])
        
        pontos = hcat(pontos, [0:vetor_h[1]:1  fill(i*vetor_h[2], elementos[1]+1)]')

    end
    
    
    # contruindo a matriz base
    linha1 = 1:elementos[1]
    linha2 = linha1 .+ 1
    linha3 = linha2 .+ (elementos[1] + 1)
    linha4 = linha1 .+ (elementos[1] + 1)
    matriz_base = [linha1 linha2 linha3 linha4]'
    
    # inicializando a matriz LG
    matriz_LG = matriz_base
    
    # inserindo as proximas linhas de elementos finitos
    for i in 2:elementos[2]
    
        valor = (i-1) * (elementos[1] + 1)
        proxima_linha = matriz_base .+ valor
        matriz_LG = hcat(matriz_LG, proxima_linha)
    
    end

    # calculando o valor de m
    valor_m = prod(elementos .- 1)

    # inicializando o vetor base
    vetor_base = 1:(elementos[1] - 1)
    
    # inicializando o vetor EQ
    vetor_EQ = vcat(valor_m+1, vetor_base, valor_m+1)

    
    # gerando os pontos internos
    for i in 2:(elementos[2]-1)
        
        valor = (i-1) * (elementos[1]-1)
        proximos_pontos = vetor_base .+ valor
        vetor_EQ = vcat(vetor_EQ, valor_m+1, proximos_pontos, valor_m+1)
        
    end

    # montando o vetor de equacao
    tampa = (valor_m + 1) * ones(Int,elementos[1]+1)
    vetor_EQ = vcat(tampa, vetor_EQ, tampa)

    
    EQLG(p,e) = vetor_EQ[matriz_LG[p,e]]

    # retorna o vetor de equacao
    return max_h, prod(elementos), pontos, valor_m, matriz_LG, EQLG, vetor_EQ
    
end

function galerkin(dados_problema::Tuple, dados_malha::Tuple) :: Vector{Float64}

    matriz_K, vetor_F = monta_global(dados_malha, gauss_quadrado(5), dados_problema)

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
    pontos_gauss = gauss[2] :: Vector{Float64}
    pesos_gauss  = gauss[3] :: Vector{Float64}
    
    # adicao do valor dos extremos
    vetor_C = cat(vetor_C, 0, dims=1)

    # inicializacao do erro
    erro = 0

    for i in 1:numero_gauss
        peso1 = pesos_gauss[i]

        for j in 1:numero_gauss
            peso2 = pesos_gauss[j]
            ponto_xi = [pontos_gauss[i], pontos_gauss[j]]
            
            for e in 1:numero_e
                pontos_e = pontos[:, lg[:,e]]
    
                jacobiano1, jacobiano2, jacobiano3, jacobiano4 = gradiente_g(ponto_xi, pontos_e)
                determinante_j = jacobiano1*jacobiano4 - jacobiano2*jacobiano3
                @assert determinante_j > 0 "Determinante do jacobiano deve ser positivo"
                
                valor_u = f_exata( funcao_g(ponto_xi, pontos_e))
                valor_e = 0
    
                for a in 1:4
                    valor_e += vetor_C[ eqlg(a, e) ] * phi(a, ponto_xi)
    
                end
                erro += determinante_j * peso1 * peso2 * (valor_u - valor_e)^2
    
            end
        end
    end
    # retorna a raiz quadrada do erro
    return sqrt(erro)

end

function plot_erro(max_h::Vector{Float64}, erros::Vector{Float64})
    max_h2 = (x -> x^2).(max_h)

    plot(max_h, max_h2, label="O(h²)", legs =:topleft, seriestype=:path, linestyle=:dash, xaxis=:log, yaxis=:log)
    plot!(max_h, erros, label="Erro_quadr ", marker=:square)
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

    lg[[3, 4], :] .= lg[[4, 3], :]
    numero_e = size(lg)[2]

    plt = plot(seriestype=:surface, xticks=[0, 1], yticks=[0, 1], zticks=[0, 1], color=:viridis, xlims=(0,1), ylims=(0,1), zlims=(0,1), colorbar=false, legend=false,camera = (30, 30, 1), size=(400, 400))
    xlabel!("x")
    ylabel!("y")
    for e in 1:numero_e
        lista_pontos = pontos[:,lg[:,e]]
        lista_altura = altura[lg[:,e]]
        coordenadas_x = transpose(reshape(lista_pontos[1,:], 2, 2))
        coordenadas_y = transpose(reshape(lista_pontos[2,:], 2, 2))
        coordenadas_z = transpose(reshape(lista_altura, 2, 2))
        
        plot!(coordenadas_x, coordenadas_y, coordenadas_z, color=:viridis)

        coordenadas_x = reshape(lista_pontos[1,:], 2, 2)
        coordenadas_y = reshape(lista_pontos[2,:], 2, 2)
        coordenadas_z = reshape(lista_altura, 2, 2)
        
        plot!(coordenadas_x, coordenadas_y, coordenadas_z, color=:viridis)
    end
    display(plt)
end

# numero_casos = 5
# lista_max_h = zeros(numero_casos)
# lista_erros = zeros(numero_casos)

# for i in 1:numero_casos
#     malha = gerador_malha([2^i, 2^i])
#     lista_erros[i] = erro( galerkin(problema(), malha) , malha, gauss_quadrado(5), problema()[end])
#     lista_max_h[i] = malha[1]

# end

# plot_erro(lista_max_h, lista_erros)

end # end modulo
