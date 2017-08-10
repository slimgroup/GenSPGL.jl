@testset "NormL1 Tests" begin

    x = Float64.([1; 2; 3; 4; 5])
    xc = complex(Float64.([1; 2; 3; 4; 5]))
    weights_scalar = 2.0
    weights_vector = weights_scalar.*ones(length(x))
    params = Dict{String, Any}()
    tau = 3.14

    @testset "NormL1_primal Tests" begin
        
        @test NormL1_primal(x, weights_scalar, params) == 30.0
        @test NormL1_primal(x, weights_vector, params) == 30.0

        @test NormL1_primal(xc, weights_scalar, params) == 30.0
        @test NormL1_primal(xc, weights_vector, params) == 30.0
            
    end

    @testset "NormL1_dual Tests" begin

        @test NormL1_dual(x, weights_scalar, params) == 2.5
        @test NormL1_dual(x, weights_vector, params) == 2.5

        @test NormL1_dual(xc, weights_scalar, params) == 2.5
        @test NormL1_dual(xc, weights_vector, params) == 2.5

    end

    @testset "NormL1_project Tests" begin

        @test isapprox(NormL1_project(x, tau, weights_scalar,
                                        spgOptions())[1], [0.0; 0.0; 0.0; 0.285; 1.285])
        @test isapprox(NormL1_project(x, tau, weights_vector,
                                        spgOptions())[1], [0.0; 0.0; 0.0; 0.285; 1.285])

        @test isapprox(NormL1_project(xc, tau, weights_scalar,
                                        spgOptions())[1], complex.([0.0; 0.0; 0.0; 0.285; 1.285]))
        @test isapprox(NormL1_project(xc, tau, weights_vector,
                                        spgOptions())[1], complex.([0.0; 0.0; 0.0; 0.285; 1.285]))

    end
end
