mutable struct Operators

    xn::Array{Float32}
    weights1D::Array{Float32}
    weights3D::CuArray{Float32}
    basis1D::Array{Float32}
    basis3D::CuArray{Float32}
    basis3Dc::CuArray{Complex{Float32}}

    cn::CuArray{Complex{Float32}}
    en::CuArray{Float32}
    an::CuArray{Complex{Float32}}
    dAn::CuArray{Complex{Float32}}
    Gn::CuArray{Complex{Float32}}
    Q::CuArray{Complex{Float32}}
end


function Operators(par::Params,cn_ini::Array{Complex{Float32}})

    u = variable(Polynomial{Rational{Int}})
    HermiteH = [basis(Hermite, i)(u) for i in 0:par.Ngrid]
    xn = sort(roots(HermiteH[par.Ngrid+1])) # Ngrid + 1 for HermiteH of degree Ngrid
    weights1D = @. 2^(par.Ngrid-1)*factorial(big(par.Ngrid))*sqrt(pi)/(par.Ngrid^2*HermiteH[par.Ngrid].(xn)^2) # Ngrid for HermiteH of degree Ngrid - 1
    weights3D = Float32.(kron(weights1D,ones(par.Ngrid^2,1)).*repeat(kron(weights1D,ones(par.Ngrid,1)),par.Ngrid,1).*kron(ones(par.Ngrid^2,1),weights1D))

    basis1D = zeros(par.Ngrid,par.cutoff+1)
    for n = 0:par.cutoff
        basis1D[:,n+1] = @. HermiteH[n+1](xn/sqrt(2))/sqrt(2^n*factorial(big(n))*sqrt(pi)); # prefactor from rescaling integral variable x = y/sqrt(2)
    end

    basis3D = zeros(par.Ngrid^3,par.modes)
    for i = 1:par.modes
        basis3D[:,i] = kron(basis1D[:,par.combinations[i,3]+1],ones(par.Ngrid^2,1)).*repeat(kron(basis1D[:,par.combinations[i,2]+1],ones(par.Ngrid,1)),par.Ngrid,1).*kron(ones(par.Ngrid^2,1),basis1D[:,par.combinations[i,1]+1])
    end

    en = sum(par.combinations,dims=2) .+ 3/2
    an = CuArray{Complex{Float32}}(undef,par.modes,1)
    dAn = CuArray{Complex{Float32}}(undef,par.modes,1)
    Gn = CuArray{Complex{Float32}}(undef,par.modes,1)
    Q = CuArray{Complex{Float32}}(undef,par.Ngrid^3,1)

  return Operators(xn, weights1D, CuArray(weights3D), basis1D, CuArray(basis3D), CuArray(complex.(basis3D)), CuArray(cn_ini), en, an, dAn, Gn, Q)
end
