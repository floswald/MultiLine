    @testset "4 lines - fails!" begin
        f1(x) = ones(length(x))
        f2(x) = 0.5x
        f3(x) = x-2
        f4(x) = 2x-8 
        x1 = collect(linspace(0.9,1.9,14))
        x2 = collect(linspace(1,7,19))
        x3 = collect(linspace(1,7,15))
        x4 = collect(linspace(1,8,25))
        X = [x1...,x2...,x3...,x4...]
        L = Line([x1...,x2...,x3...,x4...],vcat(f1(x1),f2(x2),f3(x3),f4(x4)))
        en = create_envelope(L)
 
        @test isa(en,Envelope)
        @test length(getx(en))==1
        @test length(gety(en))==1
        @test length(gets(en))==0
        @test length(getr(en))==1
        @test length(en.L) == 7-3

        upper_env!(en)
        @test issorted(getx(en))
        @test gets(en) == [Point(2.0,1.0),Point(4.0,2.0),Point(6.0,4.0)]
        xx = unique(sort(vcat(X,vcat([gets(en)[i].x for i in 1:3]...))))
        @test getx(en) == xx
        yy = reshape(vcat([[f1(ix) f2(ix) f3(ix) f4(ix)] for ix in xx]...),length(xx),4)
        @test isapprox(gety(en)[:],findmax(yy,2)[1][:])
        println(en.removed)


    end
