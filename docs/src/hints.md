# Usage Hints

## Catalyst
`TikhonovFenichelReductions.jl` can use reaction networks defined with
`Catalyst.jl` as input as of version `v0.2.7`. 
Note, that this is implemented as an extension and is therefore only supported
by Julia `1.9+` (this prevents unnecessary precompilation when `Catalyst.jl` is
not loaded).

Here is an example:
```@example 2
using Catalyst 
using TikhonovFenichelReductions

# define enzyme kinetics as reaction network
rn = @reaction_network begin
  @parameters e₀ k₁ k₋₁ k₂
  @species S(t) C(t) 
  k₁*e₀, S --> C 
  k₁, S+C --> 2S
  k₋₁, C --> S 
  k₂, C --> 0
end

# Make `ReductionProblem`
problem_catalyst = ReductionProblem(rn, 1)
```

## ``\LaTeX`` 
All of the `print_...` functions can be used with the keyword `latex=true` to
produce LaTeX code (this relies on `Latexify.jl` internally).
The output strings can also be written to arbitrary `IO` streams, which allows
to write directly to a file.
```@example 2 
print_system(problem_catalyst; latex=true)
```
