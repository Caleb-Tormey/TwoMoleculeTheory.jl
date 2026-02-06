using TwoMoleculeTheory
using Test
using StaticArrays

@testset "TwoMoleculeTheory.jl Tests" begin
    
    @testset "Data Structures" begin
        # 1. Grid
        grid = RadialGrid(2048, 0.1)
        @test length(grid.r) == 2048
        @test grid.dk * grid.dr * grid.N ≈ π

        # 2. Topology (Polyethylene Case)
        ch2 = SiteType(:CH2, 3.93, 0.0739, 14.0)
        pe_seq = ones(Int, 24)
        pe_bonds = [SVector(i, i+1) for i in 1:23]
        pe_topo = MolecularTopology([ch2], pe_seq, pe_bonds)

        @test count_site_types(pe_topo) == 1
        @test count_monomers(pe_topo) == 24

        # 3. Topology (Trimer Case)
        site_E = SiteType(:End, 1.0, 1.0, 1.0)
        site_M = SiteType(:Mid, 1.0, 1.0, 1.0)
        trimer_topo = MolecularTopology([site_E, site_M], [1, 2, 1], [SVector(1,2), SVector(2,3)])

        @test count_site_types(trimer_topo) == 2
    end
end