# ReplicateBE
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
# Licence: GNU General Public License v3.0

using Test, CSV, DataFrames

include("testdata.jl")

@testset "  Basic mixed model test                        " begin
    #df = CSV.read(IOBuffer(minibe)) |> DataFrame
    df[!,:var] = float.(df[!,:var])
    be = ReplicateBE.rbe(df, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10);
    e1 = be.fixed.est[6]
    @test be.fixed.est[6]  ≈  -0.0791666 atol=1E-5
    @test be.fixed.se[6] ≈   0.09037378448083119 atol=1E-5
    @test be.reml  ≈  10.065238638105903 atol=1E-5
    ci = ReplicateBE.confint(be, 0.1, expci = false, inv = false)
    @test ci[5][1] ≈  -0.25791330363201714 atol=1E-5
    @test ci[5][2] ≈   0.09957997029868393 atol=1E-5

    be = ReplicateBE.rbe(df, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10, memopt = false);
    @test be.fixed.est[6] == e1

    io = IOBuffer();
    Base.show(io, be)
    Base.show(io, ci)
    Base.show(io, be.design)
    Base.show(io, be.fixed)
    Base.show(io, be.typeiii)
    Base.show(io, ReplicateBE.estimate(be, [0 0 0 0 0 1]))

end

@testset "  #4 QA 1 Bioequivalence 2x2x4, UnB, NC Dataset " begin
    #REML 530.14451303
    #SE 0.04650
    #DF 208
    #df = CSV.read(IOBuffer(be4)) |> DataFrame
    df4[!,:var1] = float.(df4[!,:var1])
    be = ReplicateBE.rbe(df4, dvar = :var1, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10);
    ci = ReplicateBE.confint(be, 0.1, expci = true, inv = true)
    @test be.reml  ≈  530.1445137281626  atol=1E-5
    @test be.fixed.se[6] ≈    0.04650123700721 atol=1E-5
    @test be.fixed.f[6]  ≈    9.78552229238432 atol=1E-5
    @test be.fixed.df[6] ≈  208.08115303672898 atol=1E-5
    @test ci[5][1] ≈    1.07104135588792 atol=1E-5
    @test ci[5][2] ≈    1.24894237034602 atol=1E-5
end

#Patterson SD, Jones B. Viewpoint: observations on scaled average bioequivalence. Pharm Stat. 2012; 11(1): 1–7. doi:10.1002/pst.498
@testset "  #5 Pub Bioequivalence Dataset                 " begin
    #REML 321.44995530 - SAS STOP!
    #df = CSV.read(IOBuffer(be5)) |> DataFrame
    df5[!,:var1] = float.(df5[!,:var1])
    be = ReplicateBE.rbe(df5, dvar = :var1, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10);
    ci = ReplicateBE.confint(be, 0.1, expci = true, inv = true)
    @test be.reml  ≈  314.2217688405106 atol=1E-5
    @test ci[5][1] ≈    1.1875472284034538 atol=1E-5
    @test ci[5][2] ≈    1.5854215760408064 atol=1E-5
    #119-159%
end

#Shumaker RC, Metzler CM. The Phenytoin Trial is a Case Study of ‘Individual’ Bioequivalence. Drug Inf J. 1998; 32(4): 1063–72
@testset "  #6 Pub Bioequivalence TTRR/RRTT Dataset       " begin
    #REML 329.25749378
    #SE 0.04153
    #DF 62
    #F 2.40
    #df = CSV.read(IOBuffer(be6)) |> DataFrame
    df6[!,:var1] = float.(df6[!,:var1])
    be = ReplicateBE.rbe(df6, dvar = :var1, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10);
    ci = ReplicateBE.confint(be, 0.1, expci = true, inv = true)
    @test be.reml  ≈  329.25749377843033 atol=1E-5
    @test be.fixed.f[6]  ≈  2.399661661708039 atol=1E-5
    @test ci[5][1] ≈    0.8754960202413755 atol=1E-5
    @test ci[5][2] ≈    1.0042930817939983 atol=1E-5
end

@testset "  #  Utils test                                 " begin

    #df = CSV.read(IOBuffer(be6)) |> DataFrame
    be = ReplicateBE.rbe(df6, dvar = :var1, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10);
    @test ReplicateBE.contrast(be, [0 0 0 0 0 1]).f[1]  ≈ 2.3996616631488368 atol=1E-5
    @test ReplicateBE.estimate(be, [0 0 0 0 0 1]).est[1] ≈ 0.06434039007812514 atol=1E-5

    lsmean = ReplicateBE.lsm(be, [0 0 0 0 0 1])
    @test lsmean[1][1] ≈ 0.0643403 atol=1E-5
    @test lsmean[2][1] ≈ 0.0415345 atol=1E-5
    lsm = ReplicateBE.emm(be, [1 1 1 1 1 0], [0 0 0 0 0 0])
    @test lsm[1][1]    ≈ 4.616254407007809     atol=1E-5
    @test lsm[2][1]    ≈ 0.08217365963420642   atol=1E-5
    @test ReplicateBE.reml2(be, [0.1, 0.2, 0.3, 0.4, 1.0]) ≈ 357.238054967491   atol=1E-5

    @test ReplicateBE.coef(be)[6]      ≈ 0.06434039007812514    atol=1E-5
    @test ReplicateBE.coefse(be)[6]    ≈ 0.041534470947138996   atol=1E-5
    @test ReplicateBE.coefnum(be)      == 6
    @test ReplicateBE.reml2(be)        ≈ 329.25749377843044     atol=1E-5
    @test ReplicateBE.design(be).obs   == 256
    @test ReplicateBE.fixed(be)[1,2]   ≈ 4.4215751512542125     atol=1E-5
    @test ReplicateBE.typeiii(be)[1,2] ≈ 4.968210074464397      atol=1E-5
end

@testset "  #  Random DataSet test                        " begin
    rds = ReplicateBE.randrbeds()
    @test size(rds)[1] == 96
    @test size(rds)[2] == 5
    @test rds[1, :sequence] == "TRTR"
    @test rds[1:4, :formulation] == ["T", "R", "T", "R"]
    @test rds[93:96, :formulation] == ["R", "T", "R", "T"]
    @test rds[5:8, :period] == ["1", "2", "3", "4"]
    rds = ReplicateBE.randrbeds(dropobs = 6)
    @test size(rds)[1] == 90

    #TRTR/RTRT
    #rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0])
    #TTRR/RRTT
    #rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "T" "R" "R"; "R" "R" "T" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0])
    #TRT/RTR
    #rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "R" "T"; "R" "T" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0])
    #TRTR/RTRT/TTRR/RRTT
    #rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1,1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T" ; "T" "T" "R" "R"; "R" "R" "T" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0])
    #TRR/RTR/RRT
    #rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1,1], design = ["T" "R" "R"; "R" "T" "R"; "R" "R" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0])
end

@testset "  #  Validation with generated datasets         " begin
    #1
    #TRTR/RTRT
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10001)
    be  = ReplicateBE.rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)         ≈ 164.61336006747743   atol=1E-5
    @test ReplicateBE.coefse(be)[6]     ≈ 0.120162             atol=1E-5
    #2
    #TRRT/RTTR
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "R" "R" "T"; "R" "T" "T" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10002)
    be  = ReplicateBE.rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)         ≈ 197.20037077620236   atol=1E-5
    @test ReplicateBE.coefse(be)[6]     ≈ 0.1067239435821625   atol=1E-5
    #3
    #TTRR/RRTT
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "T" "R" "R"; "R" "R" "T" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10003)
    be  = ReplicateBE.rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)         ≈ 149.25493525435516   atol=1E-5
    @test ReplicateBE.coefse(be)[6]     ≈ 0.0875493081845683   atol=1E-5
    #4
    #TRTR/RTRT/TRRT/RTTR
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1,1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T" ; "T" "R" "R" "T"; "R" "T" "T" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10004)
    be  = ReplicateBE.rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)         ≈ 178.85707596709256   atol=1E-5
    @test ReplicateBE.coefse(be)[8]     ≈ 0.1126566088472447   atol=1E-5

    #5
    #TRRT/RTTR/TTRR/RRTT
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1,1,1], design = ["T" "R" "R" "T"; "R" "T" "T" "R" ; "T" "T" "R" "R"; "R" "R" "T" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10005)
    be  = ReplicateBE.rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)         ≈ 173.3987026483377   atol=1E-5
    @test ReplicateBE.coefse(be)[8]     ≈ 0.126898491641026    atol=1E-5
    #6
    #TRTR/RTRT/TTRR/RRTT
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1,1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T" ; "T" "T" "R" "R"; "R" "R" "T" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10006)
    be  = ReplicateBE.rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)         ≈ 149.2780196484186   atol=1E-5
    @test ReplicateBE.coefse(be)[8]     ≈ 0.109788602721634   atol=1E-5
    #7
    #TRT/RTR
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "R" "T"; "R" "T" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10007)
    be  = ReplicateBE.rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)         ≈ 138.0448457610444   atol=1E-5
    @test ReplicateBE.coefse(be)[5]     ≈ 0.125891677227639   atol=1E-5
    #8
    #TRR/RTT
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "R" "R"; "R" "T" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10008)
    be  = ReplicateBE.rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)         ≈ 126.8935514618381   atol=1E-5
    @test ReplicateBE.coefse(be)[5]     ≈ 0.099838635201318   atol=1E-5
    #9
    #TR/RT/TT/RR
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1,1,1], design = ["T" "R"; "R" "T"; "T" "T"; "R" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0, 0.0], periodcoef = [0.0, 0.0], formcoef = [0.0, 0.0], seed = 10009)
    be  = ReplicateBE.rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)         ≈ 85.85655641168267   atol=1E-5
    @test ReplicateBE.coefse(be)[6]     ≈ 0.201099945972322   atol=1E-5
    #10
    #TRR/RTR/RRT
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1,1], design = ["T" "R" "R"; "R" "T" "R"; "R" "R" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10010)
    be  = ReplicateBE.rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)         ≈ 131.24074835974176   atol=1E-5
    @test ReplicateBE.coefse(be)[6]     ≈ 0.10385830218446117  atol=1E-5
    #11
    #TRR/RTR
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "R" "R"; "R" "T" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10011)
    be  = ReplicateBE.rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)         ≈ 129.06198795781413   atol=1E-5
    @test ReplicateBE.coefse(be)[5]     ≈ 0.1096846988112254   atol=1E-5

end
