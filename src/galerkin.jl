"""
Resolve numericamente EDPs por meio do metodo de Galerkin no espaco 2D.

# Parâmetros
- `dados_problema::Vector`: Vector com os dados de entrada:
    * `1 alpha::Float64                    `: Constante do problema;
    * `2 beta::Float64                     `: Constante do problema;
    * `3 intervalo_esp::Vector{Vector{Int}}`: Intervalo do espaco analizado;
    * `4 funcao::Function                  `: Funcao f.
- `dados_malha::Vector   `: Vetor com os dados estruturais:
    * `lista_elementos::Vector{Vector{Int}}  `: Lista com a quantidade de elementos;
    * `lista_pontos::Vector{Matrix{Float64}} `: Lista com as coordenadas dos pontos da malha;
    * `lista_m::Vector{Int}                  `: Lista com o valor de m em cada caso;
    * `lista_lg::Vector{Matrix{Int}}         `: Lista com a matriz lg de cada caso;
    * `lista_eq::Vector{Vector{Int}}         `: Lista com o vetor de eq de cada caso.
- `f_exata::Function     `: Funcao exata usada no problema (``Opcional``).

# Retorna
- Se a lista de elementos tem um elemento:
    * O plot da solucao aproximada.
- Se a lista de elementos tem mais de um elemento e a solucao exata é dada:
    * O plot da convergencia do erro.
"""
function galerkin(dados_problema::Vector, dados_malha::Vector, f_exata::Function=identity)
    
    # dados do probema
    alpha         = dados_problema[1] ::Float64
    beta          = dados_problema[2] ::Float64
    funcao_f      = dados_problema[3] ::Function
    # dados da malha
    vetor_h         = dados_malha[1] :: Vector{Float64}
    lista_elementos = dados_malha[2] :: Vector{Int}
    lista_pontos    = dados_malha[3] :: Vector{Matrix{Float64}}
    lista_m         = dados_malha[4] :: Vector{Int}
    lista_lg        = dados_malha[5] :: Vector{Matrix{Int}}
    lista_eq        = dados_malha[6] :: Vector{Vector{Int}}

    # quantidade de casos analizados
    quantidade_casos = length(lista_elementos)
    
    # inicializacao o vetor com todos os erros
    vetor_erro = zeros(quantidade_casos)

    for x in 1:quantidade_casos
        # dados
        numero_e = lista_elementos[x]   # numero de elementos totais
        malha    = lista_pontos[x]      # pontos da malha
        valor_m  = lista_m[x]           # valor de m
        lg       = lista_lg[x]          # matriz localglobal(LG)
        eq       = lista_eq[x]          # matriz equacao(EQ)
        
        # composicao da EQ( LG )
        EqLg(ii,jj) = eq[ lg[ii, jj] ]

        # definindo os pontos e pesos de Gauss
        numero_gauss, pontos_gauss, pesos_gauss = gauss_triangular(7)
        gauss = [numero_gauss, pontos_gauss, pesos_gauss]

        # dados muito usados
        dados_padrao = [malha, numero_e, valor_m, lg, EqLg]
        
        # matriz K(m x m)
        matriz_K = matriz_global(dados_padrao, gauss, alpha, beta)
        
        # vetor F(m x 1)
        vetor_F = vetor_global(dados_padrao, gauss, funcao_f, alpha, beta)

        # calculando C
        vetor_C = matriz_K \ vetor_F

        if quantidade_casos == 1 && f_exata != identity
            # imagem = Plots.surface(malha[1,1:end], malha[2,1:end], zeros(length(malha[1,1:end])), xlabel="X", ylabel="Y", zlabel="Z", title="Superfície 3D", color=:rainbow, colorbar=false)
            # display(imagem)

        elseif quantidade_casos > 1 && f_exata != identity
            # calculando os erros
            vetor_erro[x] = calculo_erro(vetor_C, dados_padrao, gauss, f_exata, alpha, beta)
            
        else
            error("Quantidade de elementos invalido, informe pelo menos 1 elementos")
            
        end
    end

    if quantidade_casos > 1 && f_exata != identity
        # ordem do erro
        vetor_h2 = vetor_h .^2

        # grafico dos erros
        plt2 = plot(vetor_h2, [vetor_h vetor_h2 vetor_erro], label=["O(h)" "O(h^2)" "Erro"], xaxis=:log, yaxis=:log, markershape=:circle)
        display(plt2)

        # retorna o grafico do erro de comvergencia elemento finito dado

    end
end