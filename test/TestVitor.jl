
includet("../src/Games_src.jl")

using .NormalForms, Plots

normalform=randomNormalForm((2,2))
plotall2(normalform) #Code that plots all relevant values

p=plot()
plotIR_Set2!(p, normalform) #Code that plots all relevant values
plotFeasible2!(p, normalform) #Code that plots all relevant values

plotFeasible2(normalform)



plotall2(normalform)


p=plotIR_Set2(normalform)
