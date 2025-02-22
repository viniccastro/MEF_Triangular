"""
Fornece as coordenadas e os pessos para 7 pontos de Gauss para elementos triangulares.

# Argumentos
- `numero_p::Int`: Numero de pontos de Gauss.

# Retorna
- O numero de pontos, as coordenadas e os pesos de cada ponto.
"""
function gauss_triangular(numero_p::Int) :: Tuple{Int, Vector{Vector{Float64}}, Vector{Float64}}
    
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