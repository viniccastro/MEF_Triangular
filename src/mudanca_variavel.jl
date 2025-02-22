"""
Calcula o valor da funcao g aplicada em um ponto. A funcao e responsavel por muda a variavel de x para xi

# Parâmetros
- `ponto_xi::Vector{Float64}       `: Ponto xi usado.
- `pontos_elemento::Matrix{Float64}`: Os 3 pontos do elemento usado.

# Retorna
- Um vetor com os valores do ponto aplicado na funcao g.
# """
function funcao_g(ponto_xi::Vector{Float64}, pontos_elemento::Matrix{Float64}) :: Vector{Float64}

    # inicializando as coordenadas
    coordenada1 = 0
    coordenada2 = 0

    # calculando das coordenadas
    for a in 1:3

        coordenada1 += pontos_elemento[1, a] * phi(a, ponto_xi)
        coordenada2 += pontos_elemento[2, a] * phi(a, ponto_xi)

    end

    # retorna o vetor com as coordenadas
    return [coordenada1, coordenada2]

end

"""
Calcula o valor da derivada da funcao g aplicada em um ponto.

# Parâmetros
- `dimencao::Int                   `: Determina em qual dimencao a derivada vai ser aplicada entre 1 à 2:
    * `1`: Primeira dimencao;
    * `2`: Segunda dimencao.
- `ponto_xi::Vector{Float64}       `: Ponto xi usado.
- `pontos_elemento::Matrix{Float64}`: Os 3 pontos do elemento usado.

# Retorna
- Um vetor com os valores do ponto aplicado na derivada da funcao g.
# """
function derivada_funcao_g(dimencao::Int, ponto_xi::Vector{Float64}, pontos_elemento::Matrix{Float64}) :: Vector{Float64}

    # inicializando as coordenadas
    coordenada1 = 0
    coordenada2 = 0

    # calculando das coordenadas
    for a in 1:3

        coordenada1 += pontos_elemento[1,a] * derivada_phi(dimencao, a, ponto_xi)
        coordenada2 += pontos_elemento[2,a] * derivada_phi(dimencao, a, ponto_xi)

    end

    # retorna o vetor com as coordenadas
    return [coordenada1, coordenada2]

end