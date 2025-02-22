@testitem "malha trangular, erro de entrada" begin
    using MEF_triangular, Test

    # Teste
    @test_throws AssertionError malha_triangular(2)
end

@testitem "malha triangular, tamanho da saida" begin
    using MEF_triangular, Test

    # parametro
    numero_e = 120

    # Teste
    saida = malha_triangular(numero_e)
    @test size(saida[1])   == (2, div((numero_e+1)*(numero_e+2), 2))
    @test size(saida[3])   == (3, numero_e^2)
    @test length(saida[2]) == 1
    @test length(saida[4]) == (div((numero_e+1)*(numero_e+2), 2))
end

@testitem "malha triangular, 3 elementos" begin
    using MEF_triangular, Test
    
    # Retorno
    pontos = [0 1 0 2 1 0 3 2 1 0 ; 0 0 1 0 1 2 0 1 2 3] .* 1/3
    lg     = [1 2 2 3 4 4 5 5 6 ; 2 4 5 5 7 8 8 9 9 ; 3 5 3 6 8 5 9 6 10]
    eq     = [2, 2, 2, 2, 1, 2, 2, 2, 2, 2]
    m      = 1

    # Teste
    saida = malha_triangular(3)
    @test saida[1] == pontos
    @test saida[2] == m
    @test saida[3] == lg
    @test saida[4] == eq
end

@testitem "malha triangular, 4 elementos" begin
    using MEF_triangular, Test

    # Retorno
    pontos = [0 1 0 2 1 0 3 2 1 0 4 3 2 1 0 ; 0 0 1 0 1 2 0 1 2 3 0 1 2 3 4] .* 1/4
    lg     = [1 2 2 3 4 4 5 5 6 7 7 8 8 9 9 10 ; 2 4 5 5 7 8 8 9 9 11 12 12 13 13 14 14 ; 3 5 3 6 8 5 9 6 10 12 8 13 9 14 10 15]
    eq     = [4, 4, 4, 4, 1, 4, 4, 2, 3, 4, 4, 4, 4, 4, 4]
    m      = 3

    # Teste
    saida = malha_triangular(4)
    @test saida[1] == pontos
    @test saida[2] == m
    @test saida[3] == lg
    @test saida[4] == eq
end

@testitem "malha triangular, 5 elementos" begin
    using MEF_triangular, Test

    # Retorno
    pontos = [0 1 0 2 1 0 3 2 1 0 4 3 2 1 0 5 4 3 2 1 0 ; 0 0 1 0 1 2 0 1 2 3 0 1 2 3 4 0 1 2 3 4 5] .* 1/5
    lg     = [1 2 2 3 4 4 5 5 6 7 7 8 8 9 9 10 11 11 12 12 13 13 14 14 15 ; 2 4 5 5 7 8 8 9 9 11 12 12 13 13 14 14 16 17 17 18 18 19 19 20 20 ; 3 5 3 6 8 5 9 6 10 12 8 13 9 14 10 15 17 12 18 13 19 14 20 15 21]
    eq     = [7, 7, 7, 7, 1, 7, 7, 2, 3, 7, 7, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7]
    m      = 6

    # Teste
    saida = malha_triangular(5)
    @test isapprox(saida[1],pontos)
    @test saida[2] == m
    @test saida[3] == lg
    @test saida[4] == eq
end

@testitem "malha quadricular, erro de entrada" begin
    using MEF_triangular, Test

    # Teste
    @test_throws AssertionError malha_quadricular(1)
end

@testitem "malha quadricular, tamanho da saida" begin
    using MEF_triangular, Test

    # parametro
    numero_e = 2

    # Teste
    saida = malha_quadricular(numero_e)
    @test size(saida[1])   == (2, (numero_e+1)^2)
    @test size(saida[3])   == (3, 2*numero_e^2)
    @test length(saida[2]) == 1
    @test length(saida[4]) == ((numero_e+1)^2)
end

@testitem "malha quadricular, 3 elementos" begin
    using MEF_triangular, Test
    
    # Retorno
    pontos = [0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 ; 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3] .* 1/3
    lg     = [1  2  2  3  3  4  5  6  6  7  7  8  9 10 10 11 11 12 ; 2  6  3  7  4  8  6 10  7 11  8 12 10 14 11 15 12 16 ; 5  5  6  6  7  7  9  9 10 10 11 11 13 13 14 14 15 15]
    eq     = [5, 5, 5, 5, 5, 1, 2, 5, 5, 3, 4, 5, 5, 5, 5, 5]
    m      = 4

    # Teste
    saida = malha_quadricular(3)
    @test saida[1] == pontos
    @test saida[2] == m
    @test saida[3] == lg
    @test saida[4] == eq
end

@testitem "malha quadricular, 4 elementos" begin
    using MEF_triangular, Test

    # Retorno
    pontos = [0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 ; 0 0 0 0 0 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4] .* 1/4
    lg     = [1  2  2  3  3  4  4  5  6  7  7  8  8  9  9 10 11 12 12 13 13 14 14 15 16 17 17 18 18 19 19 20 ; 2  7  3  8  4  9  5 10  7 12  8 13  9 14 10 15 12 17 13 18 14 19 15 20 17 22 18 23 19 24 20 25 ; 6  6  7  7  8  8  9  9 11 11 12 12 13 13 14 14 16 16 17 17 18 18 19 19 21 21 22 22 23 23 24 24]
    eq     = [10, 10, 10, 10, 10, 10,  1,  2,  3, 10, 10,  4,  5,  6, 10, 10,  7,  8,  9, 10, 10, 10, 10, 10, 10]
    m      = 9

    # Teste
    saida = malha_quadricular(4)
    @test saida[1] == pontos
    @test saida[2] == m
    @test saida[3] == lg
    @test saida[4] == eq
end

@testitem "malha quadricular, 5 elementos" begin
    using MEF_triangular, Test

    # Retorno
    pontos = [0 1 2 3 4 5 0 1 2 3 4 5 0 1 2 3 4 5 0 1 2 3 4 5 0 1 2 3 4 5 0 1 2 3 4 5 ; 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5] .* 1/5
    lg     = [1  2  2  3  3  4  4  5  5  6  7  8  8  9  9 10 10 11 11 12 13 14 14 15 15 16 16 17 17 18 19 20 20 21 21 22 22 23 23 24 25 26 26 27 27 28 28 29 29 30 ; 2  8  3  9  4 10  5 11  6 12  8 14  9 15 10 16 11 17 12 18 14 20 15 21 16 22 17 23 18 24 20 26 21 27 22 28 23 29 24 30 26 32 27 33 28 34 29 35 30 36 ; 7  7  8  8  9  9 10 10 11 11 13 13 14 14 15 15 16 16 17 17 19 19 20 20 21 21 22 22 23 23 25 25 26 26 27 27 28 28 29 29 31 31 32 32 33 33 34 34 35 35]
    eq     = [17, 17, 17, 17, 17, 17, 17,  1,  2,  3,  4, 17, 17,  5,  6,  7,  8, 17, 17,  9, 10, 11, 12, 17, 17, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17]
    m      = 16

    # Teste
    saida = malha_quadricular(5)
    @test isapprox(saida[1],pontos)
    @test saida[2] == m
    @test saida[3] == lg
    @test saida[4] == eq
end