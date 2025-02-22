@testitem "vetor local, malha com 3 elementos" begin
    using MEF_triangular, Test
    
    # parametros
    numero_g, pontos_g, pesos_g = gauss_triangular(4)
    gauss   = [numero_g, pontos_g, pesos_g]
    vetor   = zeros(3)
    pontos1 = [0 1 1/3 ; 0 0 1/3]
    pontos2 = [1/3 1 0 ; 1/3 0 1]
    pontos3 = [0 1/3 0 ; 0 1/3 1]
    funcao_teste1(x, a, b) = 1.0
    funcao_teste2(x, a, b) = x[1] + x[2]

    # testes
    # funcao f constante
    vetor_local!(vetor, pontos1, gauss, funcao_teste1, 0.0, 0.0)
    @test vetor ≈ [1/18, 1/18, 1/18]
    vetor_local!(vetor, pontos2, gauss, funcao_teste1, 0.0, 0.0)
    @test vetor ≈ [1/18, 1/18, 1/18]
    vetor_local!(vetor, pontos3, gauss, funcao_teste1, 0.0, 0.0)
    @test vetor ≈ [1/18, 1/18, 1/18]

    # funcao f linear
    vetor_local!(vetor, pontos1, gauss, funcao_teste2, 0.0, 0.0)
    @test vetor ≈ [0.02314814814814815, 0.03703703703703703, 0.0324074074074074]
    vetor_local!(vetor, pontos2, gauss, funcao_teste2, 0.0, 0.0)
    @test vetor ≈ [0.046296296296296294, 0.050925925925925944, 0.05092592592592596]
    vetor_local!(vetor, pontos3, gauss, funcao_teste2, 0.0, 0.0)
    @test vetor ≈ [0.023148148148148147, 0.032407407407407406, 0.037037037037037035]
end

@testitem "vetor global, malha com 3 elementos" begin
    using MEF_triangular, Test

    # parametros
    numero_g, pontos_g, pesos_g = gauss_triangular(4)
    gauss    = [numero_g, pontos_g, pesos_g]
    l_pontos = [0 1 0 1/3 ; 0 0 1 1/3]
    lg       = [1 4 1 ; 2 2 4 ; 4 3 3]
    eq       = [2, 2, 2, 1]
    e_totais = 3
    valor_m  = 1
    eqlg(pt, el) = eq[ lg[pt, el] ]
    dados    = [l_pontos, e_totais, valor_m, lg, eqlg]
    funcao_teste1(x, a, b) = 1.0
    funcao_teste2(x, a, b) = x[1] + x[2]
    
    @test vetor_global(dados, gauss, funcao_teste1, 0.0, 0.0) ≈ [1/6]
    @test vetor_global(dados, gauss, funcao_teste2, 0.0, 0.0) ≈ [1/9]
end

@testitem "vetor local, malha com 9 elementos" begin
    using MEF_triangular, Test
    
    # parametros
    numero_g, pontos_g, pesos_g = gauss_triangular(4)
    gauss   = [numero_g, pontos_g, pesos_g]
    vetor   = zeros(3)
    pontos1 = [1/3 2/3 1/3 ; 0/3 0/3 1/3]
    pontos2 = [1/3 1/3 0/3 ; 0/3 1/3 1/3]
    pontos3 = [0/3 1/3 0/3 ; 1/3 1/3 2/3]
    pontos4 = [2/3 2/3 1/3 ; 0/3 1/3 1/3]
    pontos5 = [1/3 2/3 1/3 ; 1/3 1/3 2/3]
    pontos6 = [1/3 1/3 0/3 ; 1/3 2/3 2/3]
    funcao_teste1(x, a, b) = 1.0
    funcao_teste2(x, a, b) = x[1] + x[2]

    # testes
    # funcao f constante
    vetor_local!(vetor, pontos1, gauss, funcao_teste1, 0.0, 0.0)
    @test vetor ≈ [1/54, 1/54, 1/54]
    vetor_local!(vetor, pontos2, gauss, funcao_teste1, 0.0, 0.0)
    @test vetor ≈ [1/54, 1/54, 1/54]
    vetor_local!(vetor, pontos3, gauss, funcao_teste1, 0.0, 0.0)
    @test vetor ≈ [1/54, 1/54, 1/54]
    vetor_local!(vetor, pontos4, gauss, funcao_teste1, 0.0, 0.0)
    @test vetor ≈ [1/54, 1/54, 1/54]
    vetor_local!(vetor, pontos5, gauss, funcao_teste1, 0.0, 0.0)
    @test vetor ≈ [1/54, 1/54, 1/54]
    vetor_local!(vetor, pontos6, gauss, funcao_teste1, 0.0, 0.0)
    @test vetor ≈ [1/54, 1/54, 1/54]

    # funcao f linear
    vetor_local!(vetor, pontos1, gauss, funcao_teste2, 0.0, 0.0)
    @test vetor ≈ [0.00925925925925926, 0.010802469135802469, 0.010802469135802465]
    vetor_local!(vetor, pontos2, gauss, funcao_teste2, 0.0, 0.0)
    @test vetor ≈ [0.007716049382716049, 0.009259259259259259, 0.007716049382716049]
    vetor_local!(vetor, pontos3, gauss, funcao_teste2, 0.0, 0.0)
    @test vetor ≈ [0.00925925925925926, 0.010802469135802469, 0.010802469135802465]
    vetor_local!(vetor, pontos4, gauss, funcao_teste2, 0.0, 0.0)
    @test vetor ≈ [0.013888888888888885, 0.015432098765432101, 0.013888888888888888]
    vetor_local!(vetor, pontos5, gauss, funcao_teste2, 0.0, 0.0)
    @test vetor ≈ [0.015432098765432093, 0.016975308641975308, 0.016975308641975315]
    vetor_local!(vetor, pontos6, gauss, funcao_teste2, 0.0, 0.0)
    @test vetor ≈ [0.013888888888888885, 0.015432098765432101, 0.013888888888888888]
end

@testitem "vetor global, malha regular com 9 elementos" begin
    using .MEF_triangular

    # parametros
    numero_g, pontos_g, pesos_g = gauss_triangular(4)
    gauss    = [numero_g, pontos_g, pesos_g]
    l_pontos = [0 1 0 2 1 0 3 2 1 0 ; 0 0 1 0 1 2 0 1 2 3] .* 1/3
    lg       = [1 2 2 3 4 4 5 5 6 ; 2 4 5 5 7 8 8 9 9 ; 3 5 3 6 8 5 9 6 10]
    eq       = [2, 2, 2, 2, 1, 2, 2, 2, 2, 2]
    e_totais = 9
    valor_m  = 1
    eqlg(pt, el) = eq[ lg[pt, el] ]
    dados    = [l_pontos, e_totais, valor_m, lg, eqlg]
    funcao_teste1(x, a, b) = 1.0
    funcao_teste2(x, a, b) = x[1] + x[2]

    # testes
    @test vetor_global(dados, gauss, funcao_teste1, 0.0, 0.0) ≈ [1/9]
    @test vetor_global(dados, gauss, funcao_teste2, 0.0, 0.0) ≈ [0.07407407407407406]
end