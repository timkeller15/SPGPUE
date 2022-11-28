mutable struct Params

    NBEC::Integer
    g::Float32
    gamma::Float32
    T::Float32
    mu::Float32
    muTF::Float32
    ETF::Float32

    Tf::Float32
    dt::Float32
    steps::Integer

    cutoff::Integer
    Ngrid::Integer
    modes::Integer
    combinations::Array{Integer}
end

function Params(;NBEC = 3e3, g = 2, gamma = 5e-1, T = 0, Tf = pi, cutoff = 20)

    muTF = (15*NBEC*g/(16*sqrt(2)*pi))^(2/5)
    ETF = 5*NBEC*muTF/7
    mu = muTF

    dt = 2*pi/3200
    steps = Int(round(Tf/dt))

    Ngrid = 2*cutoff+1
    modes = Int(sum(@. 0.5*((0:cutoff) + 1)*((0:cutoff) + 2)))

    combinations = hcat(kron(0:cutoff,ones((cutoff+1)^2,1)),repeat(kron(0:cutoff,ones(cutoff+1,1)),cutoff+1,1),kron(ones((cutoff+1)^2,1),0:cutoff))
    combinations = transpose(hcat([Int.(combinations[i,:]) for i=1:size(combinations,1) if sum(combinations[i,:])<=cutoff]...))

  return Params(NBEC, g, gamma, T, mu, muTF, ETF, Tf, dt, steps, cutoff, Ngrid, modes, combinations)
end
