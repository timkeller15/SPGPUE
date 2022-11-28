function draw_wf(par::Params,cn::CuArray{Complex{Float32}},plotgrid::Int64,posmax::Float64)

    u = variable(Polynomial{Rational{Int}})
    HermiteH = [basis(Hermite, i)(u) for i in 0:par.Ngrid]
    dx = 2*posmax/plotgrid
    x = collect(1:plotgrid).*dx .- posmax
    plotbasis = Float32.(zeros(plotgrid,par.cutoff+1))
    for n = 0:par.cutoff
        plotbasis[:,n+1] = @. HermiteH[n+1](x)*exp(-0.5*x^2)/sqrt(2^n*factorial(big(n))*sqrt(pi))
    end
    wf = Complex{Float32}.(zeros(plotgrid^3,1))
    for i = 1:par.modes
        wf .+= cn[i].*kron(plotbasis[:,par.combinations[i,3]+1],ones(plotgrid^2,1)).*repeat(kron(plotbasis[:,par.combinations[i,2]+1],ones(plotgrid,1)),plotgrid,1).*kron(ones(plotgrid^2,1),plotbasis[:,par.combinations[i,1]+1])
    end
    wf = reshape(wf,plotgrid,plotgrid,plotgrid)
    return x, wf
end



function draw_current(par::Params,cn::CuArray{Complex{Float32}},plotgrid::Int64,posmax::Float64,dim::Int64)

    u = variable(Polynomial{Rational{Int}})
    HermiteH = [basis(Hermite, i)(u) for i in 0:par.Ngrid]
    dx = 2*posmax/plotgrid
    x = collect(1:plotgrid).*dx .- posmax
    plotbasis = Float32.(zeros(plotgrid,par.cutoff+1))
    for n = 0:par.cutoff
        plotbasis[:,n+1] = @. HermiteH[n+1](x)*exp(-0.5*x^2)/sqrt(2^n*factorial(big(n))*sqrt(pi))
    end

    # how to handle a^dagger for last mode?
    dn = Tridiagonal(-sqrt.(par.combinations[1:end-1,dim]) , zeros(par.modes), sqrt.(par.combinations[1:end-1,dim] .+ 1))*cn
    wf, current = Complex{Float32}.(zeros(plotgrid^3,1)), Complex{Float32}.(zeros(plotgrid^3,1))
    for i = 1:par.modes
        state = kron(plotbasis[:,par.combinations[i,3]+1],ones(plotgrid^2,1)).*repeat(kron(plotbasis[:,par.combinations[i,2]+1],ones(plotgrid,1)),plotgrid,1).*kron(ones(plotgrid^2,1),plotbasis[:,par.combinations[i,1]+1])
        wf .+= cn[i].*state
        current .+= dn[i].*state
    end

    current = @. real(1im*wf*current/sqrt(2))
    current = reshape(current,plotgrid,plotgrid,plotgrid)
    return x, current
end
