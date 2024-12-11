
@test isnothing(print_tfpv(prob, sf_separations))
@test isnothing(print_tfpv(prob, sf_separations; idx=[1,2,4]))
@test isnothing(print_tfpv(prob, sf_separations; latex=true, idx=Bool[1,1,0,1,0,0]))

@test isnothing(print_varieties(V, dim_V))
@test isnothing(print_varieties(V, dim_V; idx=[1,2,4]))
@test isnothing(print_varieties(V, dim_V; latex=true, idx=Bool[1,1,0,1,0,0]))

@test isnothing(print_results(prob, sf_separations, V, dim_V))
@test isnothing(print_results(prob, sf_separations, V, dim_V; idx=[1,2,4]))
@test isnothing(print_results(prob, sf_separations, V, dim_V; idx=Bool[1,1,0,1,0,0]))
