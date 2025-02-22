"""
Gera a matriz local de cada elemento finito.

# Parâmetros
- `matrix_local!::Matrix           `: Matriz que ira armazenar as informacoes.
- `pontos_elemento::Matrix{Float64}`: Os 3 pontos do elemento usado.
- `gauss::Vector                   `: Valores de Gauss;
    * `1 : numero::Int            `: Numero de gauss;
    * `2 : pontos::Vector{Float64}`: Pontos de gauss;
    * `3 : pesos::Vector{Float64} `: Pesos de gauss.
- `alpha::Float64                  `: Constante do problema.
- `beta::Float64                   `: Constante do problema.

# Retorna
- A matriz local.
"""
function matriz_local!(matrix_local!::Matrix, pontos_elemento::Matrix{Float64}, gauss::Vector, alpha::Float64, beta::Float64)
    # valores de Gauss
    numero_gauss = gauss[1] :: Int
    pontos_gauss = gauss[2] :: Vector{Vector{Float64}}
    pesos_gauss  = gauss[3] :: Vector{Float64}

    # limpando a matriz
    fill!(matrix_local!, 0)

    # componentes do jacobiano
    x21, y21 = pontos_elemento[:, 2] - pontos_elemento[:, 1]
    x31, y31 = pontos_elemento[:, 3] - pontos_elemento[:, 1]
    x23, y23 = pontos_elemento[:, 2] - pontos_elemento[:, 3]
    x32, y32 = pontos_elemento[:, 3] - pontos_elemento[:, 2]

    # determinante jacobiano
    determinante_j = (x21 * y31) - (x31 * y21)
    @assert determinante_j > 0 "O determinante jacobiano deve ser positivo" # verificador

    # matriz constante
    valor11 = +(x32 * x32) + (y23 * y23)
    valor21 = -(x32 * x31) + (y23 * y31)
    valor31 = +(x32 * x21) - (y23 * y21)
    valor22 = +(x31 * x31) + (y31 * y31)
    valor23 = -(x31 * x21) - (y31 * y21)
    valor33 = +(x21 * x21) + (y21 * y21)
    matrix_local! .= alpha / (2 * determinante_j) * [valor11 valor21 valor31 ; valor21 valor22 valor23 ; valor31 valor23 valor33]
    
    # somatorio de aproximacao da integral
    for i in 1:numero_gauss
        peso = pesos_gauss[i]
        ponto_xi = pontos_gauss[i]

        for b in 1:3
            phi_b = phi(b, ponto_xi)

            for a in 1:3
                # atualizando a matriz local
                matrix_local![a, b] += beta * determinante_j * peso * phi_b * phi(a, ponto_xi)

            end
        end
    end
    # atualiza a matriz local do elemento finito

end

"""
Monta a matriz global apartir das matrizes locais de cada elemento.

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
- `alpha::Float64      `: Constante do problema.
- `beta::Float64       `: Constante do problema.

# Retorna
- A matriz global.
"""
function matriz_global(dados_padrao::Vector, gauss::Vector, alpha::Float64, beta::Float64) :: Matrix{Float64}
    # dados padroes
    lista_pontos     = dados_padrao[1] :: Matrix{Float64}
    elementos_totais = dados_padrao[2] :: Int
    valor_m          = dados_padrao[3] :: Int
    Lg               = dados_padrao[4] :: Matrix{Int}
    EqLg             = dados_padrao[5] :: Function
    
    # inicializacao da matriz global
    matriz_global = spzeros(valor_m+1, valor_m+1)

    # inicializacao da matriz local
    matriz_local = zeros(3,3)
    
    for e in 1:elementos_totais
        # defininindo os pontos do elemento
        pontos_elemento = lista_pontos[:, Lg[:,e]]
        
        # calculando e atualizando a matriz local
        matriz_local!(matriz_local, pontos_elemento, gauss, alpha, beta)

        # alocando as locais na global
        for b in 1:3
            j = EqLg(b, e)

            for a in 1:3
                i = EqLg(a, e)
                matriz_global[i, j] += matriz_local[a, b]
                
            end
        end
    end
    # retorna a matriz global
    return matriz_global[1:valor_m, 1:valor_m]
    
end