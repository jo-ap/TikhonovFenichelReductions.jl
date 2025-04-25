# ReductionProblem 
@test isnothing(print_tfpvs(problem, tfpvs))
@test isnothing(print_tfpvs(problem, tfpvs; idx=[1,2,4]))
@test isnothing(print_tfpvs(problem, tfpvs; latex=true, idx=Bool[1,1,0,1,0,0]))

@test isnothing(print_varieties(varieties))
@test isnothing(print_varieties(varieties; idx=[1,2,4]))
@test isnothing(print_varieties(varieties; latex=true, idx=Bool[1,1,0,1,0,0]))

@test isnothing(print_results(problem, tfpvs, varieties))
@test isnothing(print_results(problem, tfpvs, varieties; idx=[1,2,4]))
@test isnothing(print_results(problem, tfpvs, varieties; idx=Bool[1,1,0,1,0,0]))

# Reduction

@test isnothing(print_system(problem))
@test isnothing(print_system(problem; latex=true))

@test isnothing(print_reduced_system(red_i))
@test isnothing(print_reduced_system(red_i; latex=true))
@test isnothing(print_reduced_system(red_i; rewrite=true))
@test isnothing(print_reduced_system(red_i; rewrite=true, latex=true))

@test isnothing(print_slow_fast(red_i))

