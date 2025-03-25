# ReductionProblem 

@test isnothing(print_tfpvs(prob, tfpvs))
@test isnothing(print_tfpvs(prob, tfpvs; idx=[1,2,4]))
@test isnothing(print_tfpvs(prob, tfpvs; latex=true, idx=Bool[1,1,0,1,0,0]))

@test isnothing(print_slow_manifolds(manifolds))
@test isnothing(print_slow_manifolds(manifolds; idx=[1,2,4]))
@test isnothing(print_slow_manifolds(manifolds; latex=true, idx=Bool[1,1,0,1,0,0]))

@test isnothing(print_results(prob, tfpvs, manifolds))
@test isnothing(print_results(prob, tfpvs, manifolds; idx=[1,2,4]))
@test isnothing(print_results(prob, tfpvs, manifolds; idx=Bool[1,1,0,1,0,0]))

# Reduction

@test isnothing(print_system(prob))
@test isnothing(print_system(prob; latex=true))

@test isnothing(print_reduced_system(red_i))
@test isnothing(print_reduced_system(red_i; latex=true))
@test isnothing(print_reduced_system(red_i; rewrite=true))
@test isnothing(print_reduced_system(red_i; rewrite=true, latex=true))

@test isnothing(print_slow_fast(red_i))

