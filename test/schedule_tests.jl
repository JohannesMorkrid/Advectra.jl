# Written by Claude
using Test

@testset "Every" begin
    @testset "construction" begin
        @test Every(10) isa Every
        @test_throws ArgumentError Every(0)
        @test_throws ArgumentError Every(-1)
    end

    @testset "fresh start" begin
        s = Every(10)
        setup!(s, 0.1, 0.0, 0.0)
        @test !isdue(s, 0, 0.0)
        @test !isdue(s, 9, 0.9)
        @test isdue(s, 10, 1.0)
        advance_schedule!(s)
        @test !isdue(s, 11, 1.1)
        @test isdue(s, 20, 2.0)
    end

    @testset "resume mid-stride" begin
        s = Every(10)
        setup!(s, 0.1, 0.0, 0.6)
        @test !isdue(s, 6, 0.6)
        @test !isdue(s, 9, 0.9)
        @test isdue(s, 10, 1.0)
    end

    @testset "resume on stride boundary" begin
        s = Every(10)
        setup!(s, 0.1, 0.0, 1.0)
        @test !isdue(s, 10, 1.0)
        @test isdue(s, 20, 2.0)
    end

    @testset "resume changed Δt warns" begin
        s = Every(10)
        @test_warn "Every" setup!(s, 0.2, 0.0, 0.6; Δt_prev=0.1)
    end

    @testset "resume same Δt no warn" begin
        s = Every(10)
        @test_nowarn setup!(s, 0.1, 0.0, 0.6; Δt_prev=0.1)
    end

    @testset "nonzero t0" begin
        s = Every(10)
        setup!(s, 0.1, 1.0, 1.0)
        @test !isdue(s, 0, 1.0)
        @test isdue(s, 10, 2.0)
    end

    @testset "stride of 1" begin
        s = Every(1)
        setup!(s, 0.1, 0.0, 0.0)
        @test !isdue(s, 0, 0.0)
        @test isdue(s, 1, 0.1)
        advance_schedule!(s)
        @test isdue(s, 2, 0.2)
    end
end

@testset "EveryT" begin
    @testset "construction" begin
        @test EveryT(1.0) isa EveryT
        @test_throws ArgumentError EveryT(0.0)
        @test_throws ArgumentError EveryT(1.0, -0.1)
        @test_throws ArgumentError EveryT(-1.0, 0.1)
    end

    @testset "fresh start" begin
        s = EveryT(1.0)
        setup!(s, 0.1, 0.0, 0.0)
        @test !isdue(s, 0, 0.0)
        @test !isdue(s, 9, 0.9)
        @test isdue(s, 10, 1.0)
        advance_schedule!(s)
        @test !isdue(s, 11, 1.1)
        @test isdue(s, 20, 2.0)
    end

    @testset "snapping with warning" begin
        s = EveryT(0.35)
        @test_warn "EveryT" setup!(s, 0.1, 0.0, 0.0)
        @test s.dt ≈ 0.4
    end

    @testset "snapping without warning" begin
        s = EveryT(0.35, 0.1)
        @test s.dt ≈ 0.4
        @test_nowarn setup!(s, 0.1, 0.0, 0.0)
        @test s.dt ≈ 0.4
    end

    @testset "snapping nearest 10 rule" begin
        s = EveryT(0.45, 0.1)
        @test s.dt ≈ 0.5
        s = EveryT(0.35, 0.1)
        @test s.dt ≈ 0.4
    end

    @testset "resume mid-interval" begin
        s = EveryT(1.0)
        setup!(s, 0.1, 0.0, 0.6)
        @test !isdue(s, 6, 0.6)
        @test isdue(s, 10, 1.0)
    end

    @testset "resume on boundary" begin
        s = EveryT(1.0)
        setup!(s, 0.1, 0.0, 1.0)
        @test !isdue(s, 10, 1.0)
        @test isdue(s, 20, 2.0)
    end

    @testset "nonzero t0" begin
        s = EveryT(1.0)
        setup!(s, 0.1, 1.0, 1.0)
        @test !isdue(s, 0, 1.0)
        @test isdue(s, 10, 2.0)
    end

    @testset "mismatched sign errors" begin
        s = EveryT(1.0)
        @test_throws ArgumentError setup!(s, -0.1, 0.0, 0.0)
    end

    @testset "no drift over many steps" begin
        s = EveryT(1.0)
        setup!(s, 0.1, 0.0, 0.0)
        for i in 1:1000
            if isdue(s, i, i * 0.1)
                advance_schedule!(s)
            end
        end
        n_firings = round(Int, 1000 * 0.1 / s.dt)
        @test s.occurrence ≈ (n_firings + 1) * s.dt
    end

    @testset "backward integration" begin
        s = EveryT(-1.0)
        setup!(s, -0.1, 10.0, 10.0)
        @test !isdue(s, 0, 10.0)
        @test isdue(s, 10, 9.0)
        advance_schedule!(s)
        @test !isdue(s, 11, 8.9)
        @test isdue(s, 20, 8.0)
    end
end

@testset "At" begin
    @testset "construction" begin
        @test_throws ArgumentError At(Float64[])
        @test_throws ArgumentError At([1.0, 0.5, 2.0])
        @test_throws ArgumentError At([1.0, 1.0, 2.0])
        @test_throws ArgumentError At([2.0, 0.5, 1.0])
        @test At([0.5, 1.0, 2.0]) isa At
        @test At([2.0, 1.0, 0.5]) isa At
        @test At(0.5) isa At
        @test At(0.5).occurrences == [0.5]
    end

    @testset "fresh start" begin
        s = At([0.5, 1.0, 2.0])
        setup!(s, 0.1, 0.0, 0.0)
        @test !isdue(s, 0, 0.0)
        @test !isdue(s, 4, 0.4)
        @test isdue(s, 5, 0.5)
        advance_schedule!(s)
        @test !isdue(s, 6, 0.6)
        @test isdue(s, 10, 1.0)
        advance_schedule!(s)
        @test isdue(s, 20, 2.0)
        advance_schedule!(s)
        @test !isdue(s, 21, 2.1)
    end

    @testset "single value constructor" begin
        s = At(0.5)
        setup!(s, 0.1, 0.0, 0.0)
        @test !isdue(s, 0, 0.0)
        @test isdue(s, 5, 0.5)
        advance_schedule!(s)
        @test !isdue(s, 6, 0.6)
    end

    @testset "resume skips past occurrences" begin
        s = At([0.5, 1.0, 2.0])
        setup!(s, 0.1, 0.0, 0.8)
        @test !isdue(s, 8, 0.8)
        @test isdue(s, 10, 1.0)
    end

    @testset "resume on occurrence" begin
        s = At([0.5, 1.0, 2.0])
        setup!(s, 0.1, 0.0, 0.5)
        @test !isdue(s, 5, 0.5)
        @test isdue(s, 10, 1.0)
    end

    @testset "skipped due to coarse Δt warns" begin
        s = At([0.5, 0.75, 1.0])
        @test_warn "At" setup!(s, 0.1, 0.0, 0.0)
    end

    @testset "no overflow past end forward" begin
        s = At([0.5, 1.0])
        setup!(s, 0.1, 0.0, 0.0)
        isdue(s, 5, 0.5) && advance_schedule!(s)
        isdue(s, 10, 1.0) && advance_schedule!(s)
        for _ in 1:10
            advance_schedule!(s)
        end
        @test !isdue(s, 100, 10.0)
        @test s.idx == length(s.occurrences) + 1
    end

    @testset "nonzero t0" begin
        s = At([1.5, 2.0, 3.0])
        setup!(s, 0.1, 1.0, 1.0)
        @test !isdue(s, 0, 1.0)
        @test isdue(s, 5, 1.5)
    end

    @testset "direction mismatch errors" begin
        s = At([2.0, 1.0, 0.5])
        @test_throws ArgumentError setup!(s, 0.1, 2.0, 2.0)
        s = At([0.5, 1.0, 2.0])
        @test_throws ArgumentError setup!(s, -0.1, 2.0, 2.0)
    end

    @testset "backward integration" begin
        s = At([2.0, 1.0, 0.5])
        setup!(s, -0.1, 3.0, 3.0)
        @test !isdue(s, 0, 3.0)
        @test isdue(s, 10, 2.0)
        advance_schedule!(s)
        @test !isdue(s, 15, 1.5)
        @test isdue(s, 20, 1.0)
        advance_schedule!(s)
        @test isdue(s, 25, 0.5)
        advance_schedule!(s)
        @test !isdue(s, 30, 0.0)
    end

    @testset "no overflow past end backward" begin
        s = At([2.0, 1.0])
        setup!(s, -0.1, 3.0, 3.0)
        isdue(s, 10, 2.0) && advance_schedule!(s)
        isdue(s, 20, 1.0) && advance_schedule!(s)
        for _ in 1:10
            advance_schedule!(s)
        end
        @test !isdue(s, 100, 0.0)
        @test s.idx == length(s.occurrences) + 1
    end
end

@testset "build_schedule" begin
    @testset "Every" begin
        s = build_schedule(every=10)
        @test s isa Every
        @test s.stride == 10
    end

    @testset "EveryT" begin
        s = build_schedule(dt=0.1)
        @test s isa EveryT
        @test s.dt ≈ 0.1
    end

    @testset "At vector" begin
        s = build_schedule(at=[0.5, 1.0])
        @test s isa At
        @test s.occurrences == [0.5, 1.0]
    end

    @testset "At scalar" begin
        s = build_schedule(at=0.5)
        @test s isa At
        @test s.occurrences == [0.5]
    end

    @testset "no schedule kwarg returns nothing" begin
        @test isnothing(build_schedule())
    end

    @testset "multiple schedule kwargs errors" begin
        @test_throws ArgumentError build_schedule(every=10, dt=0.1)
        @test_throws ArgumentError build_schedule(every=10, at=0.5)
        @test_throws ArgumentError build_schedule(dt=0.1, at=0.5)
    end

    @testset "ignores non-schedule kwargs" begin
        s = build_schedule(every=10, some_other_kwarg=42)
        @test s isa Every
    end
end