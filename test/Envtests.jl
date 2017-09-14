

@testset "Envelope" begin
    @testset "Constructors" begin
        n = 15
        x1 = collect(linspace(0,10,n))
        L1 = Line(x1,x1)

        en = Envelope(L1)
        @test isa(en,Envelope)
        @test length(en.L)==0
        @test length(en.env)==n
        @test length(gets(en))==0
        @test length(getr(en))==1
        @test length(getr(en)[1])==0
        @test eltype(getr(en)[1])==Point{Float64}

        en = Envelope([L1,L1])
        @test isa(en,Envelope)
        @test length(en.L)==2
        @test length(en.env)==1
        @test length(gets(en))==0
        @test length(getr(en))==1
        @test length(getr(en)[1])==0
        @test eltype(getr(en)[1])==Point{Float64}
        @test isa(getr(en)[1],Vector{Point{Float64}})
    end
end

@testset "Testing Envelopes over 2 lines" begin 
    @testset "upper_env test 1" begin
        n = 15
        x1 = collect(linspace(0,10,n))
        x2 = collect(linspace(-1,9,n))
        L1 = Line(x1,x1)
        L2 = Line(x2,ones(n)*5)
        e = Envelope([L1,L2])
        upper_env!(e)
        @test issorted(getx(e))
        @test getx(e) == sort(unique(vcat(x1,x2)))
        @test gets(e)[1].x ≈ 5.0
        @test gets(e)[1].y ≈ 5.0
    end
    @testset "upper_env test 2" begin
        n = 15
        x1 = collect(linspace(0,10,n))
        L1 = Line(x1,x1)
        L2 = Line(x1,ones(n)*5)
        e = Envelope([L1,L2])
        upper_env!(e)
        @test getx(e) == x1
        @test issorted(getx(e))
        @test gets(e)[1].x ≈ 5.0
        @test gets(e)[1].y ≈ 5.0
    end
    @testset "upper_env test 3" begin
        n = 15
        x1 = collect(linspace(0,10,n))
        L1 = Line(x1,x1)
        L2 = Line(x1,5.0+0.3*x1)
        e = Envelope([L1,L2])
        upper_env!(e)
        @test getx(e) == x1
        @test issorted(getx(e))
        @test gets(e)[1].x ≈ 5.0/0.7
        @test gets(e)[1].y ≈ 5.0/0.7
    end
    @testset "upper_env test 4" begin
        n = 15
        x1 = collect(linspace(0,10,n))
        x2 = collect(linspace(-1,9,n))
        L1 = Line(x1,x1)
        L2 = Line(x2,5.0+0.3*x2)
        e = Envelope([L1,L2])
        upper_env!(e)
        @test issorted(getx(e))
        @test getx(e) == sort(unique(vcat(x1,x2)))
        @test gets(e)[1].x ≈ 5.0/0.7
        @test gets(e)[1].y ≈ 5.0/0.7
    end
    @testset "upper_env test: decreasing " begin
        n = 15
        x1 = collect(linspace(0,10,n))
        x2 = collect(linspace(-1,9,n))
        L1 = Line(x1,x1[end:-1:1])
        L2 = Line(x2,ones(n)*5)
        e = Envelope([L1,L2])
        upper_env!(e)
        @test issorted(getx(e))
        @test getx(e) == sort(unique(vcat(x1,x2)))
        @test gets(e)[1].x ≈ 5.0
        @test gets(e)[1].y ≈ 5.0
    end
end
@testset "Testing Envelopes over more lines" begin 
    @testset "upper_env test 1" begin
        n = 15
        x1 = collect(linspace(0,10,n))
        x2 = collect(linspace(-1,9,n))
        x3 = collect(linspace(-0.1,10,n))
        L1 = Line(x1,x1)
        L2 = Line(x2,ones(n)*5)
        L3 = Line(x3,(x3.^2)/8)
        e = Envelope([L1,L2,L3])
        upper_env!(e)
        @test issorted(getx(e))
        @test getx(e) == sort(unique(vcat(x1,x2,x3,gets(e)[2].x)))
        @test length(gets(e)) == 2
        @test gets(e)[1].x ≈ 5.0
        @test gets(e)[1].y ≈ 5.0
        @test isapprox(gets(e)[2].x,8.0,rtol=0.01)
        @test isapprox(gets(e)[2].y,8.0,rtol=0.01)
    end
    @testset "upper_env test 2" begin
        n = 10
        x1 = collect(linspace(1,10,n))
        L1 = Line(x1,x1.*(x1.<6))
        L2 = Line(x1,ones(n)*5.*(x1.>4))
        L3 = Line(x1,x1-3)
        e = Envelope([L1,L2,L3])
        upper_env!(e)
        @test issorted(getx(e))
        @test getx(e) == sort(unique(vcat(x1,[gets(e)[i].x for i in 1:2])))
        @test length(gets(e)) == 2
        @test gets(e)[1].x == 5.0
        @test gets(e)[1].y == 5.0
        @test gets(e)[2].x == 8.0
        @test gets(e)[2].y == 5.0
    end
end

@testset "create_envelope" begin
    @testset "test 1: simple" begin
        x = [1,2,3,1.5,2.1,2.9]
        y = [1,1.5,1.7,1.2,1.8,2.1]
        L1 = Line(x,y)
        e = create_envelope(L1)
        @test isa(e,Envelope)
        @test length(getx(e))==1
        @test length(gety(e))==1
        @test length(gets(e))==0
        @test length(getr(e))==1
        @test length(e.L) == 3
        @test eltype(e.L)== Line{Float64}

        @test extrema(e.L[1].x) == (1.0,3.0)
        @test extrema(e.L[2].x) == (1.5,3.0)
        @test extrema(e.L[3].x) == (1.5,2.9)
    end
    @testset "test 2: zig-zag" begin
        x = [1,2,3,4,5.0,0,2,3,4,5]
        y = [2,1,2,1,2.0,1,2,1,2,1]
        L = Line(x,y)
        e = create_envelope(L)
        @test isa(e,Envelope)
        @test length(getx(e))==1
        @test length(gety(e))==1
        @test length(gets(e))==0
        @test length(getr(e))==1
        @test length(e.L) == 3

        upper_env!(e)
        @test issorted(getx(e))
        @test gety(e) == [3,2,3-4/3,2,1.5,2,1.5,2,1.5,2]
        @test getx(e) == [0,1,4/3,collect(2:0.5:5)...]
        @test gets(e)[1].x == 4/3

    end
    @testset "4 lines" begin
        f1(x) = ones(length(x))
        f2(x) = 0.5x
        f3(x) = x-2
        f4(x) = 2x-8 
        x1 = [1.0,2.0]
        x2 = collect(linspace(1,7,15))
        L = Line([x1...,vcat([x2 for i in 1:3]...)...],vcat(f1(x1),f2(x2),f3(x2),f4(x2)))
        en = create_envelope(L)
        @test isa(en,Envelope)
    end

end



