
#= To do:
- Polyhedron with prob restrictions
- Polyhedron with incentive A
- Polyhedron with incentive Bool
- Intersection of three
- linear projection
- plot
 =#
 using LazySets, Plots, LinearAlgebra

struct SetofActionDistributions
    S::LazySet{N} where N
    normalform::NormalForm

    function SetofActionDistributions(S::LazySet{N} where N, normalform::NormalForm)
        dim(S)==prod(normalform.nMoves) ? new(S, normalform) : error("Inputs with incompatible dimensions.")
    end

end

function PolyhedronProbConstraints(normalform::NormalForm)::SetofActionDistributions
    dimension=prod(normalform.nMoves)
    P0=VPolytope( Matrix{Float64}(I,dimension, dimension) )
    return SetofActionDistributions(P0, normalform)
end

#    KEY MATRICES
#= This section construct matrices that can be multiplied by vector of probabilities to calculate payoffs from deviations =#

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

    halfSpaceVector=HalfSpace.(
        arrayOfHalfspacesDeviation[ findall(!iszero, arrayOfHalfspacesDeviation) ],
        0.
    )
    return SetofActionDistributions(HPolytope(halfSpaceVector), normalform) 
end

### IC Player 2

function matrixCalcPayoff2(recommended::Int64, deviation::Int64, normalform::NormalForm)
    payoff2=normalform.payoffMatList[2]
    swapedMat = [
                            (i[2]==recommended).*payoff2[deviation,i[1]] for i in CartesianIndices(payoff2)
    ]
    return swapedMat[:]
end

function polyhedronIC_2(normalform::NormalForm)
    arrayOfHalfspacesDeviation = [
        matrixCalcPayoff2(i, j, normalform) - matrixCalcPayoff2(i, i, normalform)
        for i in 1:normalform.nMoves[2], j in 1:normalform.nMoves[2]
    ]

    halfSpaceVector=HalfSpace.(
        arrayOfHalfspacesDeviation[ findall(!iszero, arrayOfHalfspacesDeviation) ],
        0.
    )
    return SetofActionDistributions(HPolytope(halfSpaceVector), normalform)  
end


### Correlated Equilibrium polytope

function CorrEqPolytope(normalform::NormalForm)
    probSet=PolyhedronProbConstraints(normalform).S
    IC1_Set=polyhedronIC_1(normalform).S
    IC2_Set=polyhedronIC_2(normalform).S

    return SetofActionDistributions(probSet ∩ IC1_Set ∩ IC2_Set, normalform ) 
end



#= TO DO:
- Complete matrix to check deviations of player 1
- construct matrices to check deviation for player 2
- Construct respective polyhedra
- Build intersection
- projection in 2 dim
- Plot =#
