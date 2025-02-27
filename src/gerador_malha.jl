"""
Gera uma manha regular em formato triangular como apenas elementos triangulares.

# Parametros
- `numero_e::Int`: Quantidade de elementos(valor informado deve ser maior iqual a 3).

# Retorno
- Matriz com as coordenadas dos pontos da malha;
- Valor do m;
- Matriz LG como a numeracao que relaciona a local com global;
- Vetor EQ com a numeracao da fucoes inportantes.
"""
function malha_triangular(numero_e::Int) :: Tuple{Float64, Int, Matrix{Float64}, Int, Matrix{Int}, Vector{Int}}

    # verificador
    @assert numero_e > 2 "valor informado deve ser maior iqual a 3"

    # tamaho dos passos espaciais
    valor_h = 1 / numero_e
    max_h   = valor_h*sqrt(2)

    # coordenadas dos pontos da malha
    pontos = transpose([0 0])
    for i in 1:(numero_e)
        
        pontos = hcat(pontos, transpose([i:-1:0 0:1:i]))

    end

    # valor do m
    valor_m = div((numero_e-2)*(numero_e-1),2)

    # matriz LG
    matriz_lg = [1 2 2 3 ; 2 4 5 5 ; 3 5 3 6]
    for l in 5:2:(2*numero_e-1)

        padrao1 = (x -> div(x,2)).(0:(l-1))
        padrao2 = (x -> div(x,2)).(1:l)
        padrao3 = (x -> mod(x,2)).(0:(l-1))

        linha1 = padrao1 .+ (matriz_lg[1,end]+1)
        linha2 = padrao2 .+ (matriz_lg[2,end]+2)

        sequencia1 = padrao1 .+ (matriz_lg[3,end]+2)
        sequencia2 = padrao1 .+ (matriz_lg[3,end-1]+2)
        linha3 = sequencia1 .*(1 .- padrao3) .+ sequencia2 .* padrao3

        matriz_lg = hcat(matriz_lg, transpose([linha1 linha2 linha3]))

    end

    # vetor EQ
    vetor_eq = ones(Int,3) .* (valor_m+1)
    for i in 1:(numero_e-2)
        
        inicio = div(i*(i-1) + 2,2)
        fim    = div(i*(i+1),2)
        vetor_eq = [vetor_eq; [valor_m+1 ; inicio:fim ; valor_m+1]]
    
    end
    
    vetor_eq = [vetor_eq; ones(Int,numero_e+1) .* (valor_m+1)]

    return max_h, numero_e^2, valor_h*pontos, valor_m, matriz_lg, vetor_eq

end

"""
Gera uma manha regular em formato quadricular como apenas elementos triangulares.

# Parametros
- `numero_e::Int`: Quantidade de elementos(valor informado deve ser maior iqual a 2).

# Retorno
- Matriz com as coordenadas dos pontos da malha;
- Valor do m;
- Matriz LG como a numeracao que relaciona a local com global;
- Vetor EQ com a numeracao da fucoes inportantes.
"""
function malha_quadricular(numero_e::Int) :: Tuple{Float64, Int, Matrix{Float64}, Int, Matrix{Int}, Vector{Int}}

    # verificador
    @assert numero_e > 1 "valor informado deve ser maior iqual a 2"

    # tamaho dos passos espaciais
    valor_h = 1 / numero_e
    max_h   = valor_h*sqrt(2)

    # coordenadas dos pontos da malha
    pontos = transpose([0:1:numero_e  fill(0, numero_e+1)])
    for i in 1:(numero_e)
        
        pontos = hcat(pontos, transpose([0:1:numero_e  fill(i, numero_e+1)]))

    end

    # valor do m
    valor_m = (numero_e-1)^2

    # matriz LG
    matriz_lg = zeros(Int,3,0)
    for l in 1:numero_e

        padrao1 = (x -> div(x,2)).(0:(2*numero_e-1))
        padrao2 = (x -> div(x,2)).(1:(2*numero_e))
        padrao3 = (x -> mod(x,2)).(0:(2*numero_e-1))

        linha1 = padrao2 .+ (1 + (l-1)*(numero_e+1))
        linha3 = padrao1 .+ (1 + l*(numero_e+1))

        sequencia1 = padrao1 .+ (linha1[1]+1)
        sequencia2 = padrao1 .+ (linha3[1]+1)
        linha2 = sequencia1 .*(1 .- padrao3) .+ sequencia2 .* padrao3

        matriz_lg = hcat(matriz_lg, transpose([linha1 linha2 linha3]))

    end

    # vetor EQ
    vetor_eq  = zeros(Int,0)
    for i in 1:(numero_e-1)
        
        vetor_eq = [vetor_eq; [valor_m+1 ; (1:(numero_e-1)) .+ ((i-1)*(numero_e-1)) ; valor_m+1]]
        
    end
    tampa = ones(Int,numero_e+1) .* (valor_m+1)
    vetor_eq = [tampa ; vetor_eq ; tampa]


    return max_h, 2*(numero_e^2), valor_h*pontos, valor_m, matriz_lg, vetor_eq

end