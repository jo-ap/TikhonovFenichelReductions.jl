using Catalyst
using Documenter
using DocumenterCitations
using TikhonovFenichelReductions

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric
)

makedocs(
  sitename="TikhonovFenichelReductions.jl",
  repo=Remotes.GitHub("jo-ap", "TikhonovFenichelReductions.jl"),
  pages=["index.md", "gettingstarted.md", "catalyst.md", "api.md", "references.md"],
  clean=true;
  plugins=[bib]
)

deploydocs(
    repo="github.com//jo-ap/TikhonovFenichelReductions.jl.git",
)
