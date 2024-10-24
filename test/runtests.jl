using CQDBasics
using Test
using LinearAlgebra

@testset "CQDBase functions" begin
    @test CQDBase.transform_vectors([[1, 2, 3], [4, 5, 6]]) == [[1, 4], [2, 5], [3, 6]]
    @test CQDBase.transform_vectors([1, 2]) == [1, 2]
    @test collect(CQDBase.spherical_to_cartesian(1, π / 2, π / 2)) ≈ [0, 1, 0]
    @test collect(CQDBase.cartesian_to_spherical(0, 1, 0)) ≈ [1, π / 2, π / 2]
    @test norm(CQDBase.get_perpendicular_norm_vector([5, π, -3])) ≈ 1
    @test CQDBase.get_perpendicular_norm_vector([2, 3, -1]) ⋅ [2, 3, -1] ≈ 0
    @test CQDBase.wrap(1.2π) ≈ 0.8π
    @test CQDBase.latex_exponential(10^3) == "10^{3}"
end