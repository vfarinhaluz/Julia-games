# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()

module NormalForms

using LazySets, Optim, Plots, Polyhedra
import LinearAlgebra: I
# import GLMakie # Needs to be activated to run 3D plots. Doesn't run with binder

export NormalForm, randomNormalForm, plotall2, plotFeasible3, miniMaxProfile, plotFeasible2!, plotFeasible2, plotIR_Set2!

# STRUCTURE AND CONSTRUCTORS:

struct NormalForm
    nPlayers::Int64
    nMoves::Tuple{Vararg{Int64}}    
    payoffMatList::Tuple{ Vararg{ Array{<:Real} } }

    function NormalForm(players, nmoves, payofflist)
        if size(payofflist[1])==nmoves && players==ndims(payofflist[1])
            new(players, nmoves, payofflist)
        else
            error("Violates dimension requirements for NormalForm.")
        end
    end

end

function NormalForm(normalFormTable::Array{<:Tuple})
    nplayers=ndims(normalFormTable)
    nmoves=size(normalFormTable)

    payoffmatlist=ntuple(
        i -> [normalFormTable[k][i] for k in CartesianIndices(normalFormTable)],
        nplayers
    )

    return NormalForm(nplayers, nmoves, payoffmatlist)
end

function randomNormalForm(nMovesList::Tuple{Vararg{<:Real}})
    nplayers=length(nMovesList)
    payoffmatlist= Tuple(rand(Float64, nMovesList) for i in 1:nplayers)
    nform=NormalForm(nplayers, nMovesList, payoffmatlist);
    display("The payoff table(s) for the random game is:")
    display(
    [
    round.(
    Tuple(nform.payoffMatList[k][i] for k in 1:nform.nPlayers),
    digits=2)
    for i in CartesianIndices(nform.payoffMatList[1])
    ],
    )
    return nform;
end

#CALCULATING VALUES

function pureMinMax2(normalform::NormalForm)
    pure_mm1= minimum(
        maximum(
            normalform.payoffMatList[1];
            dims=1
        ),
        dims=2
    )
    pure_mm2= minimum(
        maximum(
            normalform.payoffMatList[2];
            dims=2
        ),
        dims=1
    )
    pureMMprofile=[pure_mm1[1], pure_mm2[1]]
    println("The pure minimax profile of this game is: ", round.(pureMMprofile,digits=2))
    return pureMMprofile
end

function brCorrelatMix(normalform::NormalForm,player::Int64,corrMix::Vector)
    payoff=normalform.payoffMatList[player]
    permutation=Tuple(
        unique( [player; 1:normalform.nPlayers] )
    )
    permutedPayoff=permutedims(payoff,permutation)
    permutedPayoff= reshape(permutedPayoff, ( size(permutedPayoff,1) , :))
    return maximum(permutedPayoff*corrMix)
end

function miniMax_i(normalform::NormalForm, player::Int64)
    uBR_i(μ) = brCorrelatMix(normalform, player, μ)
    parametrize_mix(z)=exp.(z)/ sum(exp.(z))
    obj(z) = uBR_i(parametrize_mix(z))

    z0=zeros( prod(normalform.nMoves) ÷ normalform.nMoves[player])
    solution=optimize(obj, z0)
    return solution.minimum
end

function miniMaxProfile(normalform::NormalForm, verbose::Bool =false)
    mmProfile=[miniMax_i(normalform, i) for i=1:normalform.nPlayers]
    if verbose
    println("The minimax profile of this game is: ",
            round.(mmProfile,digits=2))
    end
    return mmProfile
end

function greaterThanSet(minPoint::Vector{<:Real})::HPolyhedron
    HPolyhedron(
        -Matrix{Float64}(I,length(minPoint), length(minPoint)),
        -minPoint
    )
end

#PLOTTING

function newplot2()
    p=plot(
        xlabel="Payoff player 1",
        ylabel="Payoff player 2",
        title="Normal Form Payoffs",
        legend = :outerright
    )
end

function plotFeasible2!(p,normalform::NormalForm)
    payoff1=normalform.payoffMatList[1]
    payoff2=normalform.payoffMatList[2]
    
    payoff_points=[
        Vector{Float64}( [payoff1[i], payoff2[i]] )
        for i in eachindex(payoff1)
    ]
    feasible = VPolygon(payoff_points)
    Plots.plot!(p,
        feasible,
        label="Feasible set"
    );
    
    Plots.scatter!(p,
        [payoff_points[i][1] for i in eachindex(payoff_points)],
        [payoff_points[i][2] for i in eachindex(payoff_points)],
        label="Pure payoffs"
    );
end

function plotFeasible2(normalform::NormalForm)
    p=newplot2()
    plotFeasible2!(p,normalform)
    display(p)
    return p
end

function plotIR_Set2!(p::Plots.Plot, normalform::NormalForm)
    plot!(
        p,
        greaterThanSet( miniMaxProfile(normalform) ),
        label="IR Set",
        legend=:bottomleft
    )
end

function plotMinimax2!(p::Plots.Plot, normalform::NormalForm)
    mm=miniMaxProfile(normalform)
    Plots.scatter!(p,
        [mm[1]],
        [mm[2]],
        label="Minimax"
    );
end

function plotall2(normalform::NormalForm)
    p=newplot2()
    plotFeasible2!(p,normalform)
    plotIR_Set2!(p,normalform)
    plotMinimax2!(p,normalform)
    # display(p)
    return p
end

function plotFeasible3(normalform::NormalForm)
    payoff1=normalform.payoffMatList[1]
    mm=miniMaxProfile(normalform)
    nplayers=normalform.nPlayers
    
    payoff_points=[
        [ normalform.payoffMatList[k][i] for k in 1:nplayers]
        for i in eachindex(payoff1)
    ]
    feasible = VPolytope(payoff_points)
    
    fig=GLMakie.Figure()
    fig=LazySets.plot3d(feasible, alpha=0.)
    display(fig)
    GLMakie.meshscatter!(
        [payoff_points[i][1] for i in eachindex(payoff1)],
        [payoff_points[i][2] for i in eachindex(payoff1)],
        [payoff_points[i][3] for i in eachindex(payoff1)],
        color=:red,
        markersize=.02
    )

    GLMakie.meshscatter!(
        [mm[1]],
        [mm[2]],
        [mm[3]],
        color=:black,
        markersize=.05
    )
    return fig
end;

# Correlated equilibria calculations
include("CorrelatedEqPolytopes.jl")

end