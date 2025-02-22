@testitem "matriz local, malha com 3 elementos" begin
    using MEF_triangular, Test
    
    # parametros
    numero_g, pontos_g, pesos_g = gauss_triangular(7)
    gauss   = [numero_g, pontos_g, pesos_g]
    matriz  = zeros(3,3)
    pontos1 = [0 1 1/3 ; 0 0 1/3]
    pontos2 = [1/3 1 0 ; 1/3 0 1]
    pontos3 = [0 1/3 0 ; 0 1/3 1]

    # testes
    # integral do alpha
    matriz_local!(matriz, pontos1, gauss, 6.0, 0.0)
    @test matriz ≈ [ 5.0  1.0 -6.0 ;  1.0  2.0 -3.0 ; -6.0 -3.0  9.0]
    matriz_local!(matriz, pontos2, gauss, 6.0, 0.0)
    @test matriz ≈ [18.0 -9.0 -9.0 ; -9.0  5.0  4.0 ; -9.0  4.0  5.0]
    matriz_local!(matriz, pontos3, gauss, 6.0, 0.0)
    @test matriz ≈ [ 5.0 -6.0  1.0 ; -6.0  9.0 -3.0 ;  1.0 -3.0  2.0]

    # integral do beta
    matriz_local!(matriz, pontos1, gauss, 0.0, 72.0)
    @test matriz ≈ [2.0  1.0  1.0 ; 1.0  2.0  1.0 ; 1.0  1.0  2.0]
    matriz_local!(matriz, pontos2, gauss, 0.0, 72.0)
    @test matriz ≈ [2.0  1.0  1.0 ; 1.0  2.0  1.0 ; 1.0  1.0  2.0]
    matriz_local!(matriz, pontos3, gauss, 0.0, 72.0)
    @test matriz ≈ [2.0  1.0  1.0 ; 1.0  2.0  1.0 ; 1.0  1.0  2.0]

    # integral do alpha + beta
    matriz_local!(matriz, pontos1, gauss, 2.0, 24.0)
    @test matriz ≈ [2.33333  0.666667 -1.666667 ;  0.666667  1.33333 -0.666667 ; -1.666667 -0.666667  3.66667] atol=0.0001
    matriz_local!(matriz, pontos2, gauss, 2.0, 24.0)
    @test matriz ≈ [6.66667 -2.666667 -2.666667 ; -2.666667  2.33333  1.666667 ; -2.666667  1.666667  2.33333] atol=0.0001
    matriz_local!(matriz, pontos3, gauss, 2.0, 24.0)
    @test matriz ≈ [2.33333 -1.666667  0.666667 ; -1.666667  3.66667 -0.666667 ;  0.666667 -0.666667  1.33333] atol=0.0001
end

@testitem "matriz global, malha com 3 elementos" begin
    using MEF_triangular, Test

    # parametros
    numero_g, pontos_g, pesos_g = gauss_triangular(7)
    gauss    = [numero_g, pontos_g, pesos_g]
    l_pontos = [0 1 0 1/3 ; 0 0 1 1/3]
    lg       = [1 4 1 ; 2 2 4 ; 4 3 3]
    eq       = [2, 2, 2, 1]
    e_totais = 3
    valor_m  = 1
    eqlg(pt, el) = eq[ lg[pt, el] ]
    dados    = [l_pontos, e_totais, valor_m, lg, eqlg]

    # testes
    @test matriz_global(dados, gauss, 6.0,  0.0) ≈ [36] 
    @test matriz_global(dados, gauss, 0.0, 72.0) ≈ [6]
    @test matriz_global(dados, gauss, 2.0, 24.0) ≈ [42/3] 
end

@testitem "matriz local, malha com 9 elementos" begin
    using MEF_triangular, Test
    
    # parametros
    numero_g, pontos_g, pesos_g = gauss_triangular(7)
    gauss   = [numero_g, pontos_g, pesos_g]
    matriz  = zeros(3,3)
    pontos1 = [1/3 2/3 1/3 ; 0/3 0/3 1/3]
    pontos2 = [1/3 1/3 0/3 ; 0/3 1/3 1/3]
    pontos3 = [0/3 1/3 0/3 ; 1/3 1/3 2/3]
    pontos4 = [2/3 2/3 1/3 ; 0/3 1/3 1/3]
    pontos5 = [1/3 2/3 1/3 ; 1/3 1/3 2/3]
    pontos6 = [1/3 1/3 0/3 ; 1/3 2/3 2/3]

    # testes
    # integral do alpha
    matriz_local!(matriz, pontos1, gauss, 6.0, 0.0)
    @test matriz ≈ [6.0 -3.0 -3.0 ; -3.0  3.0  0.0 ; -3.0  0.0  3.0]
    matriz_local!(matriz, pontos2, gauss, 6.0, 0.0)
    @test matriz ≈ [3.0 -3.0  0.0 ; -3.0  6.0 -3.0 ;  0.0 -3.0  3.0]
    matriz_local!(matriz, pontos3, gauss, 6.0, 0.0)
    @test matriz ≈ [6.0 -3.0 -3.0 ; -3.0  3.0  0.0 ; -3.0  0.0  3.0]
    matriz_local!(matriz, pontos4, gauss, 6.0, 0.0)
    @test matriz ≈ [3.0 -3.0  0.0 ; -3.0  6.0 -3.0 ;  0.0 -3.0  3.0]
    matriz_local!(matriz, pontos5, gauss, 6.0, 0.0)
    @test matriz ≈ [6.0 -3.0 -3.0 ; -3.0  3.0  0.0 ; -3.0  0.0  3.0]
    matriz_local!(matriz, pontos6, gauss, 6.0, 0.0)
    @test matriz ≈ [3.0 -3.0  0.0 ; -3.0  6.0 -3.0 ;  0.0 -3.0  3.0]

    # integral do beta
    matriz_local!(matriz, pontos1, gauss, 0.0, 72.0)
    @test matriz ≈ [0.666667  0.333333  0.333333 ; 0.333333  0.666667  0.333333 ; 0.333333  0.333333  0.666667] atol=0.00001
    matriz_local!(matriz, pontos2, gauss, 0.0, 72.0)
    @test matriz ≈ [0.666667  0.333333  0.333333 ; 0.333333  0.666667  0.333333 ; 0.333333  0.333333  0.666667] atol=0.00001
    matriz_local!(matriz, pontos3, gauss, 0.0, 72.0)
    @test matriz ≈ [0.666667  0.333333  0.333333 ; 0.333333  0.666667  0.333333 ; 0.333333  0.333333  0.666667] atol=0.00001
    matriz_local!(matriz, pontos4, gauss, 0.0, 72.0)
    @test matriz ≈ [0.666667  0.333333  0.333333 ; 0.333333  0.666667  0.333333 ; 0.333333  0.333333  0.666667] atol=0.00001
    matriz_local!(matriz, pontos5, gauss, 0.0, 72.0)
    @test matriz ≈ [0.666667  0.333333  0.333333 ; 0.333333  0.666667  0.333333 ; 0.333333  0.333333  0.666667] atol=0.00001
    matriz_local!(matriz, pontos6, gauss, 0.0, 72.0)
    @test matriz ≈ [0.666667  0.333333  0.333333 ; 0.333333  0.666667  0.333333 ; 0.333333  0.333333  0.666667] atol=0.00001

    # integral do alpha + beta
    matriz_local!(matriz, pontos1, gauss, 2.0, 24.0)
    @test matriz ≈ [2.22222 -0.888889 -0.888889 ; -0.888889  1.22222  0.111111 ; -0.888889  0.111111  1.22222] atol=0.00001
    matriz_local!(matriz, pontos2, gauss, 2.0, 24.0)
    @test matriz ≈ [1.22222 -0.888889  0.111111 ; -0.888889  2.22222 -0.888889 ;  0.111111 -0.888889  1.22222] atol=0.00001
    matriz_local!(matriz, pontos3, gauss, 2.0, 24.0)
    @test matriz ≈ [2.22222 -0.888889 -0.888889 ; -0.888889  1.22222  0.111111 ; -0.888889  0.111111  1.22222] atol=0.00001
    matriz_local!(matriz, pontos4, gauss, 2.0, 24.0)
    @test matriz ≈ [1.22222 -0.888889  0.111111 ; -0.888889  2.22222 -0.888889 ;  0.111111 -0.888889  1.22222] atol=0.00001
    matriz_local!(matriz, pontos5, gauss, 2.0, 24.0)
    @test matriz ≈ [2.22222 -0.888889 -0.888889 ; -0.888889  1.22222  0.111111 ; -0.888889  0.111111  1.22222] atol=0.00001
    matriz_local!(matriz, pontos6, gauss, 2.0, 24.0)
    @test matriz ≈ [1.22222 -0.888889  0.111111 ; -0.888889  2.22222 -0.888889 ;  0.111111 -0.888889  1.22222] atol=0.00001
end

@testitem "matriz global, malha com 9 elementos" begin
    using MEF_triangular, Test

    # parametros
    numero_g, pontos_g, pesos_g = gauss_triangular(7)
    gauss    = [numero_g, pontos_g, pesos_g]
    l_pontos = [0 1 0 2 1 0 3 2 1 0 ; 0 0 1 0 1 2 0 1 2 3] .* 1/3
    lg       = [1 2 2 3 4 4 5 5 6 ; 2 4 5 5 7 8 8 9 9 ; 3 5 3 6 8 5 9 6 10]
    eq       = [2, 2, 2, 2, 1, 2, 2, 2, 2, 2]
    e_totais = 9
    valor_m  = 1
    eqlg(pt, el) = eq[ lg[pt, el] ]
    dados    = [l_pontos, e_totais, valor_m, lg, eqlg]

    # testes
    @test matriz_global(dados, gauss, 6.0, 0.0) ≈ [24]
    @test matriz_global(dados, gauss, 0.0, 72.0) ≈ [4]
    @test matriz_global(dados, gauss, 2.0, 24.0) ≈ [28/3]
end