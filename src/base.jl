"""
Calcula o valor da funcao phi aplicada em um ponto espacial do elemento baseado em um parametro dado.

# Parâmetros
- `phi_n::Int               `: Determina qual funcao phi vai ser aplicada entre 1 à 3:
    * `1`: Funcao phi_1;
    * `2`: Funcao phi_2;
    * `3`: Funcao phi_3.
- `ponto_xi::Vector{Float64}`: Ponto xi usado.

# Retorna
- O valor do ponto aplicado na funcao phi escolhida.
"""
function phi(phi_n::Int, ponto_xi::Vector{Float64}) :: Float64
    if phi_n == 1
        return 1 - ponto_xi[1] - ponto_xi[2]

    elseif phi_n == 2
        return ponto_xi[1]

    elseif phi_n == 3
        return ponto_xi[2]

    else               
        error("Parametro 'phi_n' invalido, informe uma valor de 1 à 3")

    end
end

"""
Calcula o valor da derivada da funcao phi aplicada em um ponto espacial do elemento baseado em um parametro dado.

# Parâmetros
- `dimencao::Int            `: Determina em qual dimencao a derivada vai ser aplicada entre 1 à 2:
    * `1`: Primeira dimencao;
    * `2`: Segunda dimencao.
- `phi_n::Int               `: Determina qual funcao phi sera derivada entre 1 à 3:
    * `1`: Funcao phi_1;
    * `2`: Funcao phi_2;
    * `3`: Funcao phi_3.
-

# Retorna
- O valor do ponto aplicado na derivada da funcao phi escolhida.
"""
function derivada_phi(dimencao::Int, phi_n::Int) :: Float64
    if dimencao == 1
        # derivada de phi em relacao a primeira dimencao
        if phi_n == 1
            return -1

        elseif phi_n == 2
            return 1

        elseif phi_n == 3
            return 0

        else                          
            error("Parametro 'phi_n' invalido, informe uma valor de 1 à 3")

        end

    elseif dimencao == 2
        # derivada de phi em relacao a segunda dimencao
        if phi_n == 1
            return -1

        elseif phi_n == 2
            return 0

        elseif phi_n == 3
            return 1

        else                        
            error("Parametro 'phi_n' invalido, informe uma valor de 1 à 3")

        end

    else                       
        error("Parametro 'dimencao' invalido, informeu um volor de 1 à 2")

    end 
end

"""
Calcula o valor do gradiente da funcao phi aplicada em um ponto espacial do elemento.

# Parâmetros
- `phi_n::Int               `: Determina qual derivada vai ser aplicada entre 1 à 3:
    * `1`: Derivada de phi_1;
    * `2`: Derivada de phi_2;
    * `3`: Derivada de phi_3.
-

# Retorna
- Um vetor com os valores do ponto aplicado no gradiente da funcao phi escolhida.
"""
function gradiente_phi(phi_n::Int) :: Vector{Float64}
    # retorna o vetor com as coordenadas
    return [derivada_phi(1, phi_n), derivada_phi(2, phi_n)]

end