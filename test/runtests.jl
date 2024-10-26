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
    @test CQDBase.sample_atom_once((1, 1), "CQD", "up", [0, 0, 0, 1])[[1, 2, 5, 6, 7]] ≈ [CQDBase.δθ / 2, CQDBase.δθ / 2, 1.1884177310293015e-5, 0.05580626719338844, 3/2]
    @test CQDBase.sample_atom_once((0.2, 2.3), "quantum", 1/2, [0, 1, 0, 0])[[1, 2, 5, 6, 7]] ≈ [π - CQDBase.δθ / 2, π - CQDBase.δθ / 2, 0.2 * 12.36e-3 / 3, 2.3 * 58.12, -1/2]
    @test collect(get_external_magnetic_fields((1.2, -0.2), 0, 1e-2, [5.2, -1, 0], "exact", "off")) ≈ [5.2, -1, 0]
    @test CQDBase.apply_branching([0.2π, 0.6π], [0.1π, 0.7π], [0, π], [1.2π, 5.5π], 3/2, "BE", "up", "HS", "B₀ dominant") == [true, false]
    @test CQDBase.is_flipped([[0.2π, 0.1π, 0, 1.2π], [0.6π, 0.7π, π, 5.5π]], 3/2, (0, 0, 1), (0, 0, 1), "BE", "up", "HS", "B₀ dominant") == [true, false]
    @test CQDBase.is_flipped([[0.2π, 0.1π, 0, 1.2π], [0.6π, 0.7π, π, 5.5π]], 3/2, (0, 0, 1), (0, 0, -1), "BE", "up", "HS", "B₀ dominant") == [false, true]
end