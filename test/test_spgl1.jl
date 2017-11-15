@testset "spgl1 tests" begin

    DTs = (Float32, Float64, Complex{Float32}, Complex{Float64})

    for DTA in DTs, DTx in DTs
        @testset "SPGL1 - $DTA $DTx" begin
            A = rand(DTA, 100, 10)
            x = rand(DTx, 10)
            b = A*x
            x_out = spgl1(A,b, options = spgOptions(verbosity = 0))[1]
            @test snr(x, x_out) > 80
        end
    end
end

