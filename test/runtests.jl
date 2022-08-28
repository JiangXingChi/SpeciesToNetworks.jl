using SpeciesToNetworks
using Test

@testset "SpeciesToNetworks" begin
    @test SpeciesCor([5,4,3,2,1],[1,2,3,4,5],"pearson") == -1.0
end