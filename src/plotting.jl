
@recipe function f(x::Envelope; removed=false)

    # defaults
    grid --> true
    xticks := true
    legend := false

    # if line array exists, plot
    if length(x.L) > 0
        for l in x.L
            @series begin
                # subplot := 1
                linetype := :line 
                linecolor := :black
                linewidth := 1
                markershape := :circle
                markerstrokecolor := :black
                markercolor := :white
                markersize := 2
                (l.x,l.y)
            end
        end
    end
    # plot envelope, if exists
    if x.env_set
        @series begin
            # subplot := 1
            linetype := :line 
            linecolor --> :red
            linewidth --> 4
            markershape := :circle
            markercolor := :white
            # markeralpha := 0.5
            markerstrokecolor := :black
            markersize := 3
            (getx(x),gety(x))
        end
        if removed
            for ir in x.removed
                if length(ir) > 0
                    @series begin
                        seriestype = :scatter
                        markershape := :rect
                        markersize := 3
                        markerstrokecolor := :black
                        markercolor := :white
                        markeralpha := 0.5
                        [ir[i].x for i in 1:length(ir)],[ir[i].y for i in 1:length(ir)]
                    end
                end
            end
        end
    end
end

function tplot1()
    n = 15
    x1 = collect(linspace(0,10,n))
    x2 = collect(linspace(-1,9,n))
    L1 = Line(x1,x1)
    L2 = Line(x2,ones(n)*5)
    e = Envelope([L1,L2])
    plot(e)
end

function tplot2()
    n = 15
    x1 = collect(linspace(0,10,n))
    x2 = collect(linspace(-1,9,n))
    L1 = Line(x1,x1)
    L2 = Line(x2,ones(n)*5)
    e = Envelope([L1,L2])
    upper_env!(e)
    plot(e)
end

function f3a()

    f1(x) = ones(length(x))
    f2(x) = 0.5x
    f3(x) = x-2
    f4(x) = 2x-8 
    x1 = collect(linspace(-1,0.9,6))
    x2 = collect(linspace(1,7,19))
    x3 = collect(linspace(2,7,15))
    x4 = collect(linspace(4,8,25))
    X = [x1...,x2...,x3...,x4...]
    L = Line([x1...,x2...,x3...,x4...],vcat(f1(x1),f2(x2),f3(x3),f4(x4)))
    en = create_envelope(L)
    return en
end
function f3b()

    f1(x) = ones(length(x))
    f2(x) = 0.5x
    f3(x) = x-2
    f4(x) = 2x-8 
    x1 = collect(linspace(-1,3.1,6))
    x2 = collect(linspace(1,7,19))
    x3 = collect(linspace(2,7,15))
    x4 = collect(linspace(4,8,25))
    X = [x1...,x2...,x3...,x4...]
    L = Line([x1...,x2...,x3...,x4...],vcat(f1(x1),f2(x2),f3(x3),f4(x4)))
    en = create_envelope(L)
    return en
end

function tplot3a()

    en = f3a()

    p1 = plot(en)

    upper_env!(en)
    p2 = plot(en,removed=true)
    plot(p1,p2)

end

function tplot3()

    en = f3b()

    p1 = plot(en)

    upper_env!(en)
    p2 = plot(en,removed=true)

    plot(p1,p2)

end

function tplot4()

    f1(x) = ones(length(x))
    f2(x) = 0.5x
    f3(x) = x-2
    f4(x) = 2x-8 
    x1 = collect(linspace(0.1,1.5,5))
    x2 = collect(linspace(1,7,19))
    x3 = collect(linspace(2,7,15))
    x4 = collect(linspace(4,8,25))
    X = [x1...,x2...,x3...,x4...]
    L = Line([x1...,x2...,x3...,x4...],vcat(f1(x1),f2(x2),f3(x3),f4(x4)))
    en = create_envelope(L)

    p1 = plot(en)

    upper_env!(en)

    p2 = plot(en,removed=true)

    plot(p1,p2)
end

