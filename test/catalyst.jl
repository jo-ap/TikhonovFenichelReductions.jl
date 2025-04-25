using Catalyst 

# define system as reaction network
rn = @reaction_network begin
  @parameters e₀ k₁ k₋₁ k₂
  @species S(t) C(t) 
  k₁*e₀, S --> C 
  k₁, S+C --> 2S
  k₋₁, C --> S 
  k₂, C --> 0
end

# check if ReductionProblem can be constructed correctly
problem_catalyst = ReductionProblem(rn, 1; idx_slow_fast = [true, true, false, true])
@test all(problem_catalyst.f .- problem.f .== 0)



