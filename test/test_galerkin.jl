@testitem "Galerkin, malha triangular" begin
    using MEF_triangular, Test

    # parametros
    pontos1, valor_m1, lg1, eq1 = malha_triangular(2^2)
    pontos2, valor_m2, lg2, eq2 = malha_triangular(2^3)
    pontos3, valor_m3, lg3, eq3 = malha_triangular(2^5)
    pontos4, valor_m4, lg4, eq4 = malha_triangular(2^6)

    todos_e      = [[2^2, 2^2], [2^3, 2^3], [2^5, 2^5], [2^6, 2^6]]
    todos_pontos = [pontos1, pontos2, pontos3, pontos4]
    todos_m      = [valor_m1, valor_m2, valor_m3, valor_m4]
    todas_lg     = [lg1, lg2, lg3, lg4]
    todas_eq     = [eq1, eq2, eq3, eq4]

    malhas = [todos_e,todos_pontos, todos_m, todas_lg, todas_eq]

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

@testitem "Galerkin, malha quadricular" begin
    using MEF_triangular, Test

    # parametros
    pontos1, valor_m1, lg1, eq1 = malha_quadricular(5)
    pontos2, valor_m2, lg2, eq2 = malha_quadricular(3)
    pontos3, valor_m3, lg3, eq3 = malha_quadricular(2^4)
    pontos4, valor_m4, lg4, eq4 = malha_quadricular(2^5)

    # todos_e      = [[2^2, 2^2], [2^3, 2^3], [2^4, 2^4], [2^5, 2^5]]
    # todos_pontos = [pontos1, pontos2, pontos3, pontos4]
    # todos_m      = [valor_m1, valor_m2, valor_m3, valor_m4]
    # todas_lg     = [lg1, lg2, lg3, lg4]
    # todas_eq     = [eq1, eq2, eq3, eq4]

    todos_e      = [[5, 5]]
    todos_pontos = [pontos1]
    todos_m      = [valor_m1]
    todas_lg     = [lg1]
    todas_eq     = [eq1]
    # display(lg1)
    # println(eq1)
    # println()

    malhas = [todos_e,todos_pontos, todos_m, todas_lg, todas_eq]
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