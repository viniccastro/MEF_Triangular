@testitem "Galerkin, malha triangular" begin
    using MEF_triangular, Test

    # parametros
    elem_1 = 2^2
    elem_2 = 2^3
    elem_3 = 2^4
    elem_4 = 2^5

    max_h1, numero_e1, pontos1, valor_m1, lg1, eq1 = malha_triangular(elem_1)
    max_h2, numero_e2, pontos2, valor_m2, lg2, eq2 = malha_triangular(elem_2)
    max_h3, numero_e3, pontos3, valor_m3, lg3, eq3 = malha_triangular(elem_3)
    max_h4, numero_e4, pontos4, valor_m4, lg4, eq4 = malha_triangular(elem_4)

    todos_max_h  = [max_h1, max_h2, max_h3, max_h4]
    todos_elemet = [numero_e1, numero_e2, numero_e3, numero_e4]
    todos_pontos = [pontos1, pontos2, pontos3, pontos4]
    todos_m      = [valor_m1, valor_m2, valor_m3, valor_m4]
    todas_lg     = [lg1, lg2, lg3, lg4]
    todas_eq     = [eq1, eq2, eq3, eq4]

    malhas = [todos_max_h, todos_elemet, todos_pontos, todos_m, todas_lg, todas_eq]

    function exemplo() :: Vector
        # dados do exemplo
        alpha         = 1.0
        beta          = 1.0
        funcao_u(x, a, b) = sin(π*x[1]) * sin(π*x[2]) * sin(π*(1 - x[1] - x[2]))
        funcao_f(x, a, b) = (4*a*π^2 + b) * sin(π*x[1]) * sin(π*x[2]) * sin(π*(1 - x[1] - x[2])) + 2*a*π^2 * sin(π*(x[1] + x[2])) * cos(π*(1 - x[1] - x[2]))
        
        # retorna os dados
        return [alpha, beta, funcao_f, funcao_u]

    end

    # teste
    galerkin(exemplo(), malhas, exemplo()[end])
end

@testitem "Galerkin, malha quadricular" begin
    using MEF_triangular, Test

    # parametros
    elem_1 = 2^2
    elem_2 = 2^3
    elem_3 = 2^4
    elem_4 = 2^5

    max_h1, numero_e1, pontos1, valor_m1, lg1, eq1 = malha_quadricular(elem_1)
    max_h2, numero_e2, pontos2, valor_m2, lg2, eq2 = malha_quadricular(elem_2)
    max_h3, numero_e3, pontos3, valor_m3, lg3, eq3 = malha_quadricular(elem_3)
    max_h4, numero_e4, pontos4, valor_m4, lg4, eq4 = malha_quadricular(elem_4)

    todos_max_h  = [max_h1, max_h2, max_h3, max_h4]
    todos_elemet = [numero_e1, numero_e2, numero_e3, numero_e4]
    todos_pontos = [pontos1, pontos2, pontos3, pontos4]
    todos_m      = [valor_m1, valor_m2, valor_m3, valor_m4]
    todas_lg     = [lg1, lg2, lg3, lg4]
    todas_eq     = [eq1, eq2, eq3, eq4]

    malhas = [todos_max_h, todos_elemet, todos_pontos, todos_m, todas_lg, todas_eq]

    function exemplo() :: Vector
        # dados do exemplo
        alpha         = 1.0
        beta          = 1.0
        funcao_f(x, a, b) = (2*a * π^2 + b) * sin(π*x[1]) * sin(π*x[2])
        funcao_u(x, a, b) = sin(π*x[1]) * sin(π*x[2])       # funcao exata

        # retorna os dados
        return [alpha, beta, funcao_f, funcao_u]

    end

    # teste
    galerkin(exemplo(), malhas, exemplo()[end])
end