using MultiLine
using Base.Test
using TestSetExtensions


@testset ExtendedTestSet "Running MultiLine tests" begin
    

    @testset "setup" begin
        x = collect(0:0.1:1)
        y = rand(11)
        m = Line(x,y)
        @test m.n == 11
        @test size(m) == (11,)
    end

    @testset "Interpolations" begin
        @testset "interpolate at log(1)" begin
            x = collect(0:0.1:1)
            y = log.(1.+x)
            m = Line(x,y)
            i = interp(m,0.0)
            @test i == 0.0
        end
    end

    @testset "getter and setters" begin
        x = [1,3,4,10]
        y = [6,0,8,1]
        l = Line(x,y)
        @test size(l) == (4,)

        @test l[1] == (1,6)
        @test isa(l[1:2],Line{Int64})
        @test l[1:2].x == [1,3]
        @test l[1:2].y == [6,0]

        setindex!(l,-1,-1,1)
        @test size(l) == (4,)
        @test l[1] == (-1,-1)

        setindex!(l,[10,11],[90,91],3:4)
        @test l[3:4].x == [10,11]
        @test l[3:4].y == [90,91]
        
    end

    @testset "Modifying methods" begin
        
    end

end


