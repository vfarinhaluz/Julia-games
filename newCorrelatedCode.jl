
#= To do:
- Polyhedron with prob restrictions
- Polyhedron with incentive A
- Polyhedron with incentive Bool
- Intersection of three
- linear projection
- plot
 =#

struct SetofActionDistributions
    S::LazySet{N} where N
    normalform::NormalForm

    function SetofActionDistributions(S, normalform)
        dim(S)==prod(normalform.nMoves) ? new(S, normalform) : error("Incompatible inputs")
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


#= TO DO:
- Complete matrix to check deviations of player 1
- construct matrices to check deviation for player 2
- Construct respective polyhedra
- Build intersection
- projection in 2 dim
- Plot =#


function Polyhedron_IC_1(normalform::NormalForm)
    payoff1=normalform.payoffMatList[1]
    return nothing
end
