"""
Gera o vetor local de cada elemento finito baseado em uma funcao dada.

# Parâmetros
- `vector_local!::Vector           `: Vetor que ira armazenar as informacoes.
- `pontos_elemento::Matrix{Float64}`: Os 4 pontos do elemento usado.
- `gauss::Vector                   `: Valores de Gauss;
    * `1 : numero::Int            `: Numero de gauss;
    * `2 : pontos::Vector{Float64}`: Pontos de gauss;
    * `3 : pesos::Vector{Float64} `: Pesos de gauss.
- `funcao::Function                `: Funcao que sera usada.
- `alpha::Float64                  `: Constante do problema.
- `beta::Float64                   `: Constante do problema.

# Retorna
- O vetor local.
# """
function vetor_local!(vector_local!::Vector, pontos_elemento::Matrix{Float64}, gauss::Vector, funcao::Function, alpha::Float64, beta::Float64)
    # valores de Gauss
    numero_gauss = gauss[1] :: Int
    pontos_gauss = gauss[2] :: Vector{Vector{Float64}}
    pesos_gauss  = gauss[3] :: Vector{Float64}

    # limpando o vetor
    fill!(vector_local!, 0)

    # componentes do jacobiano
    x21, y21 = pontos_elemento[:, 2] - pontos_elemento[:, 1]
    x31, y31 = pontos_elemento[:, 3] - pontos_elemento[:, 1]

    # determinante jacobiano
    determinante_j = (x21 * y31) - (x31 * y21)
    @assert determinante_j > 0 "O determinante jacobiano deve ser positivo" # verificador

    # somatorio de aproximacao da integral
    for i in 1:numero_gauss
        peso = pesos_gauss[i]
        ponto_xi = pontos_gauss[i]

        # valor da funcao
        valor_f = funcao( funcao_g(ponto_xi, pontos_elemento) , alpha, beta)

        for a in 1:3
            vector_local![a] += determinante_j * peso * valor_f * phi(a, ponto_xi)

        end
    end
    # retorna o vetor local do elemento finito

end

"""
Monta o vetor global apartir dos vetores locais de cada elemento.

# Parâmetros
- `dados_padrao::Vector`: dados muito usados:
    * `1 : lista_pontos::Matrix{Float64}`: Lista com as coordenadas dos pontos da malha usada;
    * `2 : elementos_totais::Int        `: Numero de elementos finitos totais;
    * `3 : valor_m::Int                 `: Numero de funcoes da base;
    * `4 : Lg::Function                 `: Matriz LG;
    * `5 : EqLg::Function               `: Funcao da composicao da EQ( LG ).
- `gauss::Vector       `: Valores de Gauss;
    * `1 : numero::Int            `: Numero de gauss;
    * `2 : pontos::Vector{Float64}`: Pontos de gauss;
    * `3 : pesos::Vector{Float64} `: Pesos de gauss.
- `funcao::Function    `: Funcao usada no problema.
- `alpha::Float64      `: Constante do problema.
- `beta::Float64       `: Constante do problema.

# Retorna
- O vetor global.
"""
function vetor_global(dados_padrao::Vector, gauss::Vector, funcao::Function, alpha::Float64, beta::Float64) :: Vector{Float64}
    # dados padroes
    lista_pontos     = dados_padrao[1] :: Matrix{Float64}
    elementos_totais = dados_padrao[2] :: Int
    valor_m          = dados_padrao[3] :: Int
    Lg               = dados_padrao[4] :: Matrix{Int}
    EqLg             = dados_padrao[5] :: Function

    # inicializacao do vetor global
    vetor_global = spzeros(valor_m+1)
    
    # inicializacao do vetor local
    vetor_local = zeros(3)
    
    # montando a vetor global
    for e in 1:elementos_totais
        pontos_elemento = lista_pontos[:, Lg[:,e]]
        
        vetor_local!(vetor_local, pontos_elemento, gauss, funcao, alpha, beta)
        
        for a in 1:3
            i = EqLg(a, e)
            vetor_global[i] += vetor_local[a]

        end
    end
    # retorna o vetor global
    return vetor_global[1:valor_m]

end