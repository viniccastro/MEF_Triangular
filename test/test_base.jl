@testitem "phi tests" begin
    using MEF_triangular, Test
    
    @test phi(1, [0.0, 0.00]) ==  1.00
    @test phi(1, [0.5, 0.25]) ==  0.25
    @test phi(1, [1.0, 1.00]) == -1.00
    
    @test phi(2, [0.0, 0.00]) ==  0.00
    @test phi(2, [0.5, 0.25]) ==  0.50
    @test phi(2, [1.0, 1.00]) ==  1.00
    
    @test phi(3, [0.0, 0.00]) ==  0.00
    @test phi(3, [0.5, 0.25]) ==  0.25
    @test phi(3, [1.0, 1.00]) ==  1.00

    @test_throws ErrorException phi(4, [0.0, 0.0])
end

@testitem "phi X phi_vs2" begin
    using MEF_triangular, Test
    using BenchmarkTools
    
    @btime phi(1, [0.5, 0.25])
    @btime phi(2, [0.5, 0.25])
    @btime phi(3, [0.5, 0.25])

    @btime phi_vs2(1, [0.5, 0.25])
    @btime phi_vs2(2, [0.5, 0.25])
    @btime phi_vs2(3, [0.5, 0.25])
end

#######################################################################################

@testitem "derivada_phi tests" begin
    using MEF_triangular, Test
    
    @test derivada_phi(1, 1, [0.0, 0.0]) == -1.00
    @test derivada_phi(1, 2, [0.0, 0.0]) ==  1.00
    @test derivada_phi(1, 3, [0.0, 0.0]) ==  0.00
    
    @test derivada_phi(2, 1, [0.0, 0.0]) == -1.00
    @test derivada_phi(2, 2, [0.0, 0.0]) ==  0.00
    @test derivada_phi(2, 3, [0.0, 0.0]) ==  1.00

    @test_throws ErrorException derivada_phi(1, 4, [0.0, 0.0])
    @test_throws ErrorException derivada_phi(3, 1, [0.0, 0.0])
end

@testitem "derivada_phi time" begin
    using MEF_triangular, Test
    using BenchmarkTools
    
    @btime derivada_phi(1, 1, [0.5, 0.25])
    @btime derivada_phi(1, 2, [0.5, 0.25])
    @btime derivada_phi(1, 3, [0.5, 0.25])
end

#######################################################################################

@testitem "gradiente_phi tests" begin
    using MEF_triangular, Test
    
    @test gradiente_phi(1, [0.0, 0.00]) == [-1.0, -1.0]
    @test gradiente_phi(2, [0.0, 0.00]) == [ 1.0,  0.0]
    @test gradiente_phi(3, [0.0, 0.00]) == [ 0.0,  1.0]
    
    @test gradiente_phi(1, [0.5, 0.25]) == [-1.0, -1.0]
    @test gradiente_phi(2, [0.5, 0.25]) == [ 1.0,  0.0]
    @test gradiente_phi(3, [0.5, 0.25]) == [ 0.0,  1.0]

    @test_throws ErrorException gradiente_phi(4, [0.0, 0.0])
end

@testitem "gradiente_phi X gradiente_phi_vs2" begin
    using MEF_triangular, Test
    using BenchmarkTools
    
    @btime gradiente_phi(1, [0.5, 0.25])
    @btime gradiente_phi(2, [0.5, 0.25])
    @btime gradiente_phi(3, [0.5, 0.25])

    @btime gradiente_phi_vs2(1, [0.5, 0.25])
    @btime gradiente_phi_vs2(2, [0.5, 0.25])
    @btime gradiente_phi_vs2(3, [0.5, 0.25])
end