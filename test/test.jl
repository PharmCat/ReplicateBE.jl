# ReplicateBE
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
# Licence: GNU General Public License v3.0

using Test, DataFrames, CSV

include("testdata.jl")

@testset "  Basic mixed model test       " begin
    df = CSV.read(IOBuffer(minibe)) |> DataFrame
    be = ReplicateBE.rbe(df, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence);
    @test be.β[6]  ≈  -0.0791666 atol=1E-5
    @test be.se[6] ≈   0.09037378448083119 atol=1E-5
    @test be.reml  ≈  10.065238638105903 atol=1E-5
end

@testset "  Bioequivalence 2x2x4 DS test " begin
    df = CSV.read(IOBuffer(be1)) |> DataFrame
    df[:, :var] = log.(df[:, :var])
    be = ReplicateBE.rbe(df, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence);
    ci = ReplicateBE.confint(be, 0.1, expci = true, inv = true)
    @test be.reml  ≈  660.0465401 atol=1E-5
    @test ci[5][1] ≈    0.7792777889433989 atol=1E-5
    @test ci[5][2] ≈    0.9810195569153635 atol=1E-5
end
