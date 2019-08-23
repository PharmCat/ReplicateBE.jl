# ReplicateBE
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
# Licence: GNU General Public License v3.0

using Test, DataFrames, CSV

include("testdata.jl")

@testset "  Basic mixed model test       " begin
    df = CSV.read(IOBuffer(minibe)) |> DataFrame
    be = ReplicateBE.rbe(df, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence);
    @test be.β[6]  ≈  -0.0791666 atol=1E-6
    @test be.se[6] ≈   0.0903709 atol=1E-6
    @test be.reml  ≈  10.0652386 atol=1E-6
end

@testset "  Bioequivalence 2x2x4 DS test " begin
    df = CSV.read(IOBuffer(be1)) |> DataFrame
    df[:, :var] = log.(df[:, :var])
    be = ReplicateBE.rbe(df, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence);
    ci = confint(be, 0.1, inv=true)
    @test be.reml  ≈  660.0465401 atol=1E-6
    @test ci[5][1] ≈    0.7792800 atol=1E-6
    @test ci[5][2] ≈    0.9810167 atol=1E-6
end
