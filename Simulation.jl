"""
# Universal Evolutionary Model for Periodical Organisms

Accompanying Julia code for the paper "Universal evolutionary model for periodical organisms"

```
Eric Goles, Ivan Slapničar and Marco A. Lardies: Universal evolutionary model for periodical organisms, arXiv:2010.00940, submitted
```

The manuscript can be found [here](https://arxiv.org/abs/2010.00940).

The code in this file (also in the Pluto notebook `Simulation_p.jl`) was used to run all simulations and produce all figures in the paper.
"""

using LinearAlgebra, Plots, Plots.PlotMeasures, Random, Primes
Random.seed!(123)
import Plots.spy
spy(A)=heatmap(A, yflip=true, legend=false, c=cgrad([:white,:gray,:red]),
	   aspectratio=1,clim=(0.0,1.0))

mkpath("figures")

# Reference images
# gcd=1
n=20
P=zeros(Float64,n,n)
for i=2:n
	for j=2:n
	    if gcd(i,j)==1
	        P[i,j]=0.5
	    end
	end
end
figG=spy(P)
xlabel!("c1")
ylabel!("c2")

# Primes
P=zeros(Float64,n,n)
for i=2:n
    for j=2:n
        if isprime(i)
            P[i,j]=0.5
        end
    end
end
figP=spy(P)
xlabel!("c1")
ylabel!("c2")

# Divisors
P=zeros(Float64,n,n)
for i=2:n
    for j=2:n
        if mod(i,j)==0
            P[i,j]=0.5
        end
    end
end
figD=spy(P)
xlabel!("c1")
ylabel!("c2")

fᵣ=plot(figG,figP,figD,layout=(1,3),titlefontsize=10,xguidefontsize=10,
    yguidefontsize=10,left_margin=[0mm 0mm],bottom_margin=[10mm 10mm],
    top_margin=[0mm 0mm],right_margin=[0mm 0mm],size=(645,200),grid=false)
savefig(fᵣ,"./figures/Reference.png")

# Local fitness functions
f₁(χ₁,χ₂,v::Tuple)= (χ₁,χ₂)==(1,0) ? v[1] : (χ₁,χ₂)==(1,1) ? v[3] : 0
f₂(χ₁,χ₂,v)= (χ₁,χ₂)==(0,1) ? v[2] : (χ₁,χ₂)==(1,1) ? v[4] : 0
# This is special case for long-time masting bamboos
α=0.8
β=1/α
f₁(χ₁,χ₂,v::Vector)= (χ₁,χ₂)==(1,0) ? v[1] : (χ₁,χ₂)==(1,1) ? α : (χ₁,χ₂)==(0,1) ? β : 0

# Global fitness functions
# v is instance, c₁ and c₂ are starting periods of species C₁ and C₂
function F₁(c₁,c₂,v)
    T=c₁*c₂
    χ₁(t)=map(Int,iszero(mod(t,c₁)))
    χ₂(t)=map(Int,iszero(mod(t,c₂)))
    fitness=0
    for t=1:T
        fitness+=f₁(χ₁(t),χ₂(t),v)
    end
    fitness/c₂
end

function F₂(c₁,c₂,v)
    T=c₁*c₂
    χ₁(t)=map(Int,iszero(mod(t,c₁)))
    χ₂(t)=map(Int,iszero(mod(t,c₂)))
    fitness=0
    for t=1:T
        fitness+=f₂(χ₁(t),χ₂(t),v)
    end
    fitness/c₁
end

# Evolution dynamics algorithm
function Evolution(c₁,c₂,v,p)
    # v is instance, c₁ and c₂ are starting periods of species C₁ and C₂,
    # p is the tuple of possible period changes
    # Initial fitnesses
    F1=F₁(c₁,c₂,v)
    F2=F₂(c₁,c₂,v)
    F1new=copy(F1)
    F2new=copy(F2)
    E=Vector{Tuple{Int32,Int32,Number,Number}}(undef,0)
    push!(E,(c₁,c₂,F1,F2))
    change=true
    steps=0
    maxsteps=30
    while change || steps<maxsteps
        # Random change of species C₁
        d₁=c₁+rand(p)
        change=false
        if d₁>1
            F1new=F₁(d₁,c₂,v)
            # Switch to new period if fitness is better
            if F1new > F1
                c₁=d₁
                F1=F1new
                F2=F₂(c₁,c₂,v)
                push!(E,(c₁,c₂,F1,F2))
                change=true
            end
        end
        # Random change of species C₂
        d₂=c₂+rand(p)
        if d₂>1
            F2new=F₂(c₁,d₂,v)
            # Switch to new period if fitness is better
            if F2new > F2
                c₂=d₂
                F2=F2new
                F1=F₁(c₁,c₂,v)
                push!(E,(c₁,c₂,F1,F2))
                change=true
            end
        end
        steps+=1
    end
    E
    return (E[end][1],E[end][2])
end

# Define fitness tuples v:
v=[(-1, 0, 1, 1), (-1, 1, 1, 1), (0, 1, 1, 1),
(0, 0, 1, 1), (-1, -1, 1, 1), (1, 1, -1, -1),
(1, 1, -1, 0), (1, 1, 1, -1), (-1, 1, 1, -1),
(-1, 1, 1, 0), (1, 0, -1, 0), (0, 1, 1, -1),
(0, 1, 1, 0), (1, 1, 1, 0), (0, -1, 1, 0),
(-1, 0, 1, 0), (1, -1, -1, 0), (1, -1, 1, 0)]

# Plots of attractors
function myPlots(n::Int, k::Int, v::Vector, quantitative::Bool=false)
    # n is size of plot
    # k defines the evolution neighbourhood p
	# v vector of tuples
    # quantitative==false means to use qualitative model (-1,0,1)

    # Define neighbourhood p
    q=collect(-k:k)
    deleteat!(q,k+1)
    p=Tuple(1q)
    # String for title
    stringp=", p=[-"*string(k)*","*string(k)*"]"

    Ev=Array{Tuple{Int32,Int32}}(undef,n,n)
    for i=1:n
        Ev[1,i]=(0,0)
        Ev[i,1]=(0,0)
    end
    # Floating-point array to plot
    FpPlot=zeros(Float64,n,n)
    # Array for figures
    fig=Array{Any}(undef,length(v))
    l=1
    for i=1:length(v)
        FpPlot=zeros(Float64,n,n)
        ν=v[i]
        if quantitative==true
			μ=Tuple([rand(1:10) for j=1:length(ν)])
            ν=μ.*v[i]
        end
        for c₁=2:n
            for c₂=2:n
                Evo=Evolution(c₁,c₂,ν,p)
                Ev[c₁,c₂]=Evo
                if Evo==(c₁,c₂)
                    if F₁(c₁,c₂,ν)>=0 && F₂(c₁,c₂,ν)>=0
                        # Both species survive
                        FpPlot[c₂,c₁]=1.0
                    else
                        # Someone dies
                        if F₁(c₁,c₂,ν)<0
                            FpPlot[c₂,c₁]=0.3
                        else
                            FpPlot[c₂,c₁]=0.3
                        end
                    end
                end
            end
        end
        # Plot
        fig[l]=spy(FpPlot)
        xlabel!("c1")
        ylabel!("c2")
        title!("v["*string(i)*"]="*string(ν)*stringp)
        l+=1
    end
    return fig
end

# Attractors for neigbourhood [-1,1]
fig1=myPlots(20,1,v)

# Attractors for neigbourhood [-4,4]
fig4=myPlots(20,4,v)

# Quantitative plots for neigbourhood [-1,1]
figQ=myPlots(20,1,v,true)

# Figures for the paper
for k=1:length(v)
    f1=fig1[k]
    f4=fig4[k]
    fQ=figQ[k]
    f=plot(f1,f4,fQ,layout=(1,3),titlefontsize=10,xguidefontsize=10,yguidefontsize=10,
        left_margin=[0mm 0mm],bottom_margin=[10mm 10mm],top_margin=[0mm 0mm],
        right_margin=[0mm 0mm],size=(645,220),grid=false)
    savefig(f,"./figures/S"*string(k)*".png")
end


# Long-period masting bamboos
# Redefine fitness vector
v₀=[[-1, 1, α, 1, β, 0, 0, 0]]
# Time period
# n=40
# k defines the evolution neighbourhood p
# k=6
# Attractors
figB=myPlots(40,6,v₀)
# Figure
fb=plot(figB[1],titlefontsize=10,xguidefontsize=10,yguidefontsize=10,
    left_margin=[0mm 0mm],bottom_margin=[0mm 0mm],top_margin=[0mm 0mm],
    right_margin=[0mm 0mm],size=(400,400),grid=:show)
savefig(fb,"./figures/Bamboos.png")
