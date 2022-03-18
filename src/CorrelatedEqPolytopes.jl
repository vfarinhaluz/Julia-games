# Uses LazySets to calculate correlated equilibria payoffs (to be added to main module)

export plotCorrelated, CorrEquilibriumPayoffs, plotCorrelated!

struct SetofActionDistributions
    S::LazySet{N} where N
    normalform::NormalForm

    function SetofActionDistributions(S::LazySet{N} where N, normalform::NormalForm)
        LazySets.dim(S)==prod(normalform.nMoves) ? new(S, normalform) : error("Inputs with incompatible dimensions.")
    end

end

# Polytope of probabilities
function PolyhedronProbConstraints(normalform::NormalForm)::SetofActionDistributions
    dimension=prod(normalform.nMoves)
    P0=VPolytope( Matrix{Float64}(I,dimension, dimension) )
    return SetofActionDistributions(P0, normalform)
end

# Polytope of incentive compatibility for player 1

function matrixCalcPayoff1(recommended::Int64, deviation::Int64, normalform::NormalForm)
    payoff1=normalform.payoffMatList[1]
    swapedMat = [
                            (i[1]==recommended).*payoff1[deviation,i[2]] for i in CartesianIndices(payoff1)
    ]
    return swapedMat[:]
end

function polyhedronIC_1(normalform::NormalForm)
    arrayOfHalfspacesDeviation = [
        matrixCalcPayoff1(i, j, normalform) - matrixCalcPayoff1(i, i, normalform)
        for i in 1:normalform.nMoves[1], j in 1:normalform.nMoves[1]
    ]

    halfSpaceVector=LazySets.HalfSpace.(
        arrayOfHalfspacesDeviation[ findall(!iszero, arrayOfHalfspacesDeviation) ],
        0.
    )
    return SetofActionDistributions(HPolytope(halfSpaceVector), normalform) 
end

# Polytope of incentive compatibility for player 2

function matrixCalcPayoff2(recommended::Int64, deviation::Int64, normalform::NormalForm)
    payoff2=normalform.payoffMatList[2]
    swapedMat = [
                            (i[2]==recommended).*payoff2[i[1],deviation] for i in CartesianIndices(payoff2)
    ]
    return swapedMat[:]
end

function polyhedronIC_2(normalform::NormalForm)
    arrayOfHalfspacesDeviation = [
        matrixCalcPayoff2(i, j, normalform) - matrixCalcPayoff2(i, i, normalform)
        for i in 1:normalform.nMoves[2], j in 1:normalform.nMoves[2]
    ]

    halfSpaceVector=LazySets.HalfSpace.(
        arrayOfHalfspacesDeviation[ findall(!iszero, arrayOfHalfspacesDeviation) ],
        0.
    )
    return SetofActionDistributions(HPolytope(halfSpaceVector), normalform)  
end

### Correlated Equilibrium polytopes (Intersection of three above)

function CorrEqPolytope(normalform::NormalForm)
    probSet=PolyhedronProbConstraints(normalform).S
    IC1_Set=polyhedronIC_1(normalform).S
    IC2_Set=polyhedronIC_2(normalform).S

    return SetofActionDistributions(probSet ∩ IC1_Set ∩ IC2_Set, normalform ) 
end

function CorrEquilibriumPayoffs(normalform::NormalForm)
    payoffMat=[
        normalform.payoffMatList[1][:]';
        normalform.payoffMatList[2][:]'
    ]
    return payoffMat * CorrEqPolytope(normalform::NormalForm).S
end

# PLOTTING

function plotCorrelated!(p,normalform::NormalForm)
    p=plot!(
        p,
        CorrEquilibriumPayoffs(normalform),
        label="Correlated eq."
    )
end

function plotCorrelated(normalform::NormalForm)
    p=newplot2(normalform)
    p=plotCorrelated!(p,normalform)
end
