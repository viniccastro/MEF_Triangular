"""
Calcula o erro entre a solucao exata e a aproximada baseada na integral.

# Parâmetros
- `vetor_C::Vector{Float64}`: Vetor com os valores de C.
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
- `f_exata::Function       `: Funcao exata usada no problema.
- `alpha::Float64          `: Constante do problema.
- `beta::Float64           `: Constante do problema.

# Retorna
- O erro em cada elemento.
"""
function calculo_erro(vetor_C::Vector{Float64}, dados_padrao::Vector, gauss::Vector, f_exata::Function, alpha::Float64, beta::Float64) :: Float64
    # dados padroes
    lista_pontos = dados_padrao[1] :: Matrix{Float64}
    numero_e     = dados_padrao[2] :: Int
    lg           = dados_padrao[4] :: Matrix{Int}
    eqlg         = dados_padrao[5] :: Function

    # adicao do valor dos extremos
    vetor_C = cat(vetor_C, 0, dims=1)
    
    # valores de Gauss
    numero_gauss = gauss[1] :: Int
    pontos_gauss = gauss[2] :: Vector{Vector{Float64}}
    pesos_gauss  = gauss[3] :: Vector{Float64}

    # inicializacao do erro
    erro = 0

    for i in 1:numero_gauss
        # definindo o peso e o ponto xi
        peso = pesos_gauss[i]
        ponto_xi = pontos_gauss[i]
        
        for e in 1:numero_e
            # defininindo os pontos do elemento
            pontos_elemento = lista_pontos[:, lg[:,e]]

            # componentes do jacobiano
            x21, y21 = pontos_elemento[:, 2] - pontos_elemento[:, 1]
            x31, y31 = pontos_elemento[:, 3] - pontos_elemento[:, 1]

            # determinante jacobiano
            determinante_j = (x21 * y31) - (x31 * y21)
            @assert determinante_j > 0 "O determinante jacobiano deve ser positivo" # verificador

            valor_u = f_exata( funcao_g(ponto_xi, pontos_elemento) , alpha, beta)
            valor_c  = 0

            for a in 1:3
                valor_c += vetor_C[ eqlg(a, e) ] * phi(a, ponto_xi)

            end
            # calcula o erro ao quadrado
            erro += determinante_j * peso * (valor_u - valor_c)^2   

        end
    end
    return sqrt(erro)

end