#!/usr/bin/julia
using ColorSchemes

cs = ColorSchemes.inferno
L = length(cs)
@info "Selected colorscheme" cs typeof(cs) L
l = 100
r = floor(Int64,L/l)
@info r
