using MultiLine
using Base.Test
using TestSetExtensions


@testset ExtendedTestSet "Running MultiLine tests" begin
    

    @testset "setup" begin
        x = collect(0:0.1:1)
        y = Vector{Float64}[]
        push!(y,log.(1.+x),2*log.(1.+x),3*log.(1.+x))
        m = Mline(x,y)
        @test m.n == 11
        @test m.m == 3
        @test m.nm == 33
    end

    @testset "interpolate at log(1)" begin
        x = collect(0:0.1:1)
        y = Vector{Float64}[]
        push!(y,log.(1.+x),2*log.(1.+x),3*log.(1.+x))
        m = Mline(x,y)
        i = interpolate(m,1.0)
        @test i == zeros(3)
    end
end


