using Test, DataFrames, CSV
include("testdata.jl")
df4 = CSV.read(IOBuffer(be4)) |> DataFrame
df5 = CSV.read(IOBuffer(be5)) |> DataFrame
df6 = CSV.read(IOBuffer(be6)) |> DataFrame
