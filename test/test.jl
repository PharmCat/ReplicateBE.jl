# ReplicateBE
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
# Licence: GNU General Public License v3.0

using Test, CSV, DataFrames, StatsBase

path    = dirname(@__FILE__)
println("Load data...")
include("testdata.jl")
println("Start test...")
@testset "  Basic mixed model test                         " begin
    be = ReplicateBE.rbe!(df0, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10);
    e1 = ReplicateBE.fixed(be).est[6]
    @test ReplicateBE.fixed(be).est[6]       ≈  -0.0791666 atol=1E-5
    @test ReplicateBE.fixed(be).se[6]        ≈   0.09037378448083119 atol=1E-5
    @test ReplicateBE.reml2(be) ≈   10.065238638105903 atol=1E-5
    @test ReplicateBE.reml2(be, ReplicateBE.theta(be), memopt = false) ≈ ReplicateBE.reml2(be) atol=1E-8
    ci = confint(be, 0.1, expci = false, inv = false)
    @test ci[end][1]            ≈  -0.25791330363201714 atol=1E-5
    @test ci[end][2]            ≈   0.09957997029868393 atol=1E-5
    ci = confint(be, 0.1; expci = true, inv = true, df = :df2)
    @test ci[end][1]            ≈  0.8750195990056241 atol=1E-5
    @test ci[end][2]            ≈  1.338891894355988 atol=1E-5
    ci = confint(be, 0.1; expci = true, inv = true, df = :contw)
    @test ci[end][1]            ≈  0.9080640550377563 atol=1E-5
    @test ci[end][2]            ≈  1.2901696108459495 atol=1E-5
    ci = confint(be; level = 0.90)
    @test ci[end][1]            ≈  -0.25791330363201714 atol=1E-5
    @test ci[end][2]            ≈   0.09957997029868393 atol=1E-5
    ci = confint(be, 0.1, df = :cont)
    @test ci[end][1]            ≈  -0.24721576677454637 atol=1E-5
    @test ci[end][2]            ≈   0.088882433441213 atol=1E-5
    @test dof(be)[end]          ≈   5.463110799437906 atol=1E-5

    be = ReplicateBE.rbe!(df0, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10, memopt = false);
    @test ReplicateBE.fixed(be).est[6]       ≈ e1 atol=1E-10
    @test ReplicateBE.estimate(be, [0 0 0 0 0 1], df = :cont, name = "Formulation")[1,4] == 8
    ci0 = confint(be)[end]
    io = IOBuffer();
    Base.show(io, be)
    Base.show(io, ci)
    Base.show(io, be.design)
    Base.show(io, ReplicateBE.fixed(be))
    #Base.show(io, be.typeiii)
    Base.show(io, ReplicateBE.estimate(be, [0 0 0 0 0 1]))
    #StatsBase.coeftable
    Base.show(io, coeftable(be))
    #StatsBase.vcov
    vcm = vcov(be)
    #POSTOPT+
    #Not recommended
    be = ReplicateBE.rbe!(df0, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10, postopt = true)
    #CONTRAST MULTIDIM+
    L = [0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0]
    t =  ReplicateBE.typeiii(be)
    c =  ReplicateBE.contrast(be, L)
    @test t[2, 5] ≈ c[1, 5]
    Base.show(io, t)

    be  = ReplicateBE.rbe!(df0, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, rholink = :sigmoid)
    ci1 = confint(be)[end]
    @test ci0[1] ≈ ci1[1] atol=1E-5
    @test ci0[2] ≈ ci1[2] atol=1E-5
    #Experimental
    be  = ReplicateBE.rbe!(df0, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, rholink = :arctgsigmoid, singlim = 1e-4)
    ci2 = confint(be)[end]

end

@testset "  #1 Patterson and Jones 2017 E 4.3              " begin
    be = ReplicateBE.rbe!(df1, dvar = :logCmax, subject = :id, formulation = :formulation, period = :period, sequence = :sequence)
    ci = confint(be, 0.1; expci = true)[end]
    @test ci[1]                             ≈  0.831853  atol=1E-4
    @test ci[2]                             ≈  1.000517  atol=1E-4
    #SPSS show -2REML = 15.013102
    @test ReplicateBE.reml2(be)             ≈ 13.441805425068509   atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.05486050426851888   atol=1E-5
    be = ReplicateBE.rbe!(df1, dvar = :logAUC, subject = :id, formulation = :formulation, period = :period, sequence = :sequence)
    ci = confint(be, 0.1; expci = true)[end]
    @test ci[1]                             ≈ 0.993140 atol=1E-5
    @test ci[2]                             ≈ 1.081396 atol=1E-5
    @test ReplicateBE.reml2(be)             ≈ -49.71856515650643   atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.024308218717669812   atol=1E-5
end

@testset "  #2 Shein-Chung Chow, Jen-pei Liu 2009 E 9.4.1  " begin
    be = ReplicateBE.rbe!(df2, dvar = :logAUC, subject = :id, formulation = :formulation, period = :period, sequence = :sequence)
    #SPSS -2REML 45.242642
    @test ReplicateBE.reml2(be)             ≈ 45.242642227970045 atol=1E-5
    ci = confint(be, 0.1; expci = true)[end]
    @test ci[1]                             ≈  0.9819173624303978 atol=1E-5
    @test ci[2]                             ≈  1.38703413072702 atol=1E-5
end

@testset "  #3 Patterson and Jones 2017 E 4.4              " begin
    be = ReplicateBE.rbe!(df3, dvar = :logCmax, subject = :id, formulation = :formulation, period = :period, sequence = :sequence)
    ci = confint(be, 0.1; expci = true)[end]
    @test ReplicateBE.reml2(be)             ≈ 433.84147581860884
    @test ci[1]                             ≈   1.322569 atol=1E-5
    @test ci[2]                             ≈   1.723750  atol=1E-5
    be = ReplicateBE.rbe!(df3, dvar = :logAUC, subject = :id, formulation = :formulation, period = :period, sequence = :sequence)
    ci = confint(be, 0.1; expci = true)[end]
    #SPSS not positive Hessian
    @test ReplicateBE.reml2(be)             ≈ 245.65265598596116
    @test ci[1]                             ≈  1.0290230615220375 atol=1E-5  #0.841779
    @test ci[2]                             ≈  1.1880048531419543 atol=1E-5  #0.971759
end

@testset "  #4 QA 1 Bioequivalence 2x2x4, UnB, NC Dataset  " begin
    #REML 530.14451303
    #SE 0.04650
    #DF 208
    #df4[!,:var1] = float.(df4[!,:var1])
    be = ReplicateBE.rbe!(df4, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10);
    ci = ReplicateBE.confint(be, 0.1, expci = true, inv = true)
    @test ReplicateBE.reml2(be)  ≈  530.1445137281626  atol=1E-5
    @test ReplicateBE.fixed(be).se[6] ≈    0.04650123700721 atol=1E-5
    @test ReplicateBE.fixed(be).f[6]  ≈    9.78552229238432 atol=1E-4
    @test ReplicateBE.fixed(be).df[6] ≈  208.08115303672898 atol=1E-2
    @test ci[end][1]     ≈    1.071047105  atol=1E-5
    @test ci[end][2]     ≈    1.248935873  atol=1E-5
end

#Patterson SD, Jones B. Viewpoint: observations on scaled average bioequivalence. Pharm Stat. 2012; 11(1): 1–7. doi:10.1002/pst.498
@testset "  #5 Patterson 2012 doi:10.1002/pst.498 AUC      " begin
    #SAS  REML 321.44995530 - SAS STOP!
    #SPSS REML 314.221769
    df5[!,:var] = float.(df5[!,:var])
    be = ReplicateBE.rbe!(df5, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10);
    ci = ReplicateBE.confint(be, 0.1, expci = true)
    @test ReplicateBE.reml2(be)                           ≈  314.2217688405106 atol=1E-5
    #1.0.4
    @test ci[end][1]                        ≈  0.6307479743996646 atol=1E-4 #1.187496 SPSS
    @test ci[end][2]                        ≈  0.8420705538500828 atol=1E-4 #1.585490 SPSS
    #1.0.5
    #@test ci[end][1]                        ≈  0.6307285753672432 atol=1E-5 #1.187496 SPSS
    #@test ci[end][2]                        ≈  0.8420964530317809 atol=1E-5 #1.585490 SPSS
    #0.8420964530317809
    #119-159%
end

#Shumaker RC, Metzler CM. The Phenytoin Trial is a Case Study of ‘Individual’ Bioequivalence. Drug Inf J. 1998; 32(4): 1063–72
@testset "  #6 Shumaker 1998 doi:10.1177/009286159803200426" begin
    #REML 329.25749378
    #SE 0.04153
    #DF 62
    #F 2.40
    #df6[!,:var1] = float.(df6[!,:var1])
    be = ReplicateBE.rbe!(df6, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10);
    ci = confint(be, 0.1, expci = true)
    @test ReplicateBE.reml2(be)  ≈  329.25749377843033 atol=1E-5
    @test ReplicateBE.fixed(be).f[end]        ≈  2.399661661708039 atol=1E-5
    @test ci[end][1]             ≈  0.994999  atol=1E-5
    @test ci[end][2]             ≈  1.143044  atol=1E-5
end

@testset "  #  Utils test                                  " begin
    be = ReplicateBE.rbe!(df6, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10);
    @test ReplicateBE.contrast(be, [0 0 0 0 0 1]).f[1]   ≈ 2.3996616631488368 atol=1E-5
    estt = ReplicateBE.estimate(be, [0 0 0 0 0 1])
    @test estt.est[1]                                    ≈ 0.06434039007812514 atol=1E-5
    @test estt.df[1]                                     ≈ 62.0                atol=1E-5
    @test estt.ul[1]                                     ≈ 0.14736661447650518 atol=1E-5

    @test ReplicateBE.contrast(be, [0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0]).df[1] ≈ 144.97812496691395 atol=1E-4

    #=
    lsmean = ReplicateBE.lsm(be, [0 0 0 0 0 1])
    @test lsmean[1][1] ≈ 0.0643403 atol=1E-5
    @test lsmean[2][1] ≈ 0.0415345 atol=1E-5
    lsm = ReplicateBE.emm(be, [1 1 1 1 1 0], [0 0 0 0 0 0])
    @test lsm[1][1]    ≈ 4.616254407007809     atol=1E-5
    @test lsm[2][1]    ≈ 0.08217365963420642   atol=1E-5
    =#
    @test ReplicateBE.reml2(be, [0.1, 0.2, 0.3, 0.4, 1.0]) ≈ 357.238054967491   atol=1E-5
    @test ReplicateBE.coefnum(be)      == 6
    @test ReplicateBE.reml2(be)        ≈ 329.25749377843044     atol=1E-5
    @test ReplicateBE.design(be).obs   == 256
    @test ReplicateBE.fixed(be)[1,2]   ≈ 4.4215751512542125     atol=1E-5
    @test ReplicateBE.typeiii(be)[1,2] ≈ 4.968210074464397      atol=1E-5
    @test ReplicateBE.optstat(be)

    @test nobs(be)             == 64
    @test coef(be)[end]                ≈ 0.06434039007812514    atol=1E-5
    @test stderror(be)[end]            ≈ 0.041534470947138996   atol=1E-5

end

@testset "  #  Random DataSet test                         " begin
    rds = ReplicateBE.randrbeds()
    @test size(rds)[1] == 96
    @test size(rds)[2] == 5
    @test rds[1, :sequence] == "TRTR"
    @test rds[1:4, :formulation] == ["T", "R", "T", "R"]
    @test rds[93:96, :formulation] == ["R", "T", "R", "T"]
    @test rds[5:8, :period] == [1, 2, 3, 4]
    rds = ReplicateBE.randrbeds(dropobs = 6)
    @test size(rds)[1] == 90
end

@testset "  #  Validation with generated datasets          " begin
    using CSV
    path    = dirname(@__FILE__)
    #println("File path: " , path)
    #println("File rds12.csv ($(path*"/csv/rds12.csv")) exist: ", isfile(path*"/csv/rds12.csv"))
    #1
    #TRTR/RTRT
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10001)
    be  = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, store_trace = true, extended_trace = true)
    @test ReplicateBE.reml2(be)             ≈ 164.61336006747743   atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.120162             atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈  0.800182 atol=1E-5
    @test ci[end][2]                        ≈  1.208955 atol=1E-5
    print("[.")
    #DF contain 46
    #DF contain form 44
    #2
    #TRRT/RTTR
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "R" "R" "T"; "R" "T" "T" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10002)
    be  = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 197.20037077620236   atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.1067239435821625   atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.877858             atol=1E-5
    @test ci[end][2]                        ≈ 1.266491             atol=1E-5
    print(".")
    #3
    #TTRR/RRTT
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "T" "R" "R"; "R" "R" "T" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10003)
    be  = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 149.25493525435516   atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.0875493081845683   atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.838738             atol=1E-5
    @test ci[end][2]                        ≈ 1.132936             atol=1E-5
    print(".")
    #4
    #TRTR/RTRT/TRRT/RTTR
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1,1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T" ; "T" "R" "R" "T"; "R" "T" "T" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10004)
    be  = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 178.85707596709256   atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.1126566088472447   atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.901150701434849    atol=1E-5  #SPSS 0.900982
    @test ci[end][2]                        ≈ 1.3260902090001896   atol=1E-5  #     1.326338
    print(".")
    #5
    #TRRT/RTTR/TTRR/RRTT
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1,1,1], design = ["T" "R" "R" "T"; "R" "T" "T" "R" ; "T" "T" "R" "R"; "R" "R" "T" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10005)
    be  = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)         ≈ 173.3987026483377   atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.126898491641026    atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.8697647630450958    atol=1E-4
    @test ci[end][2]                        ≈ 1.3438521036812334    atol=1E-4
    print(".")
    #6
    #TRTR/RTRT/TTRR/RRTT
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1,1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T" ; "T" "T" "R" "R"; "R" "R" "T" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10006)
    be  = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 149.2780196484186   atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.109788602721634   atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.7280915876704905    atol=1E-5
    @test ci[end][2]                        ≈ 1.0609077632653205   atol=1E-5
    print(".")
    #7
    #TRT/RTR
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "R" "T"; "R" "T" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10007)
    be  = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 138.0448457610444   atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.125891677227639   atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.824266             atol=1E-5
    @test ci[end][2]                        ≈ 1.269283             atol=1E-5
    #DF contain 23
    #DF contain form 44
    print(".")
    #8
    #TRR/RTT
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "R" "R"; "R" "T" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10008)
    be  = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 126.8935514618381   atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.099838635201318   atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.8957103074455708             atol=1E-5
    @test ci[end][2]                        ≈ 1.2534507598650542             atol=1E-5
    print(".")
    #9
    #TR/RT/TT/RR
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1,1,1], design = ["T" "R"; "R" "T"; "T" "T"; "R" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0, 0.0], periodcoef = [0.0, 0.0], formcoef = [0.0, 0.0], seed = 10009)
    be  = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 85.85655641168267   atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.201099945972322   atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.5581809471131624  atol=1E-5
    @test ci[end][2]                        ≈ 1.1201862625078982  atol=1E-5
    print(".")
    #10
    #TRR/RTR/RRT
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1,1], design = ["T" "R" "R"; "R" "T" "R"; "R" "R" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10010)
    be  = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 131.24074835974176   atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.10385830218446117  atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.8192235933294734  atol=1E-5
    @test ci[end][2]                        ≈ 1.1698775266746582  atol=1E-5
    print(".")
    #11
    #TRR/RTR
    #SPSS REML 129.655513635398
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "R" "R"; "R" "T" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10011)
    be  = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 129.06198795781413   atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.1096846988112254   atol=1E-5
    ci = confint(be, 0.1, expci = true)

    #1.0.4
    #@test ci[end][1]                        ≈ 0.8719537380504332   atol=1E-5
    #@test ci[end][2]                        ≈ 1.2621358539962937   atol=1E-5
    #1.0.5
    #@test ci[end][1]                        ≈ 0.8719023753585434   atol=1E-3
    #@test ci[end][2]                        ≈ 1.2622102048603623   atol=1E-3
    #1.0.11
    #@test ci[end][1]                        ≈ 0.8719313261195931   atol=1E-5
    #@test ci[end][2]                        ≈ 1.2621682956584102   atol=1E-5
    #J1.4.0
    @test ci[end][1]                        ≈ 0.8719537380504332   atol=1E-4
    @test ci[end][2]                        ≈ 1.2621358539962937   atol=1E-4

    #DF contain 74
    #DF contain form 92
    print(".")
    #Unbalanced + dropouts
    #12
    #TRTR/RTRT
    #SPSS DF 42.3382785451983
    rds = CSV.File(path*"/csv/rds12.csv") |> DataFrame
    #rds = ReplicateBE.randrbeds(;n=48, sequence=[1,2], design = ["T" "R" "T" "R"; "R" "T" "R" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], dropobs = 20, seed = 10012)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 329.76454300481595    atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.07412076135308149   atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.9223372675334348   atol=1E-5
    @test ci[end][2]                        ≈ 1.183433775497828   atol=1E-5
    print(".")
    #13
    #TRRT/RTTR
    rds = CSV.File(path*"/csv/rds13.csv") |> DataFrame
    #rds = ReplicateBE.randrbeds(;n=48, sequence=[1,2], design = ["T" "R" "R" "T"; "R" "T" "T" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], dropobs = 20, seed = 100013)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 305.21958862123165    atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.07880499833587161   atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.9109673193677009   atol=1E-5
    @test ci[end][2]                        ≈ 1.1871890135539473   atol=1E-5
    print(".")
    #14
    #TTRR/RRTT
    rds = CSV.File(path*"/csv/rds14.csv") |> DataFrame
    #rds = ReplicateBE.randrbeds(;n=48, sequence=[1,2], design = ["T" "T" "R" "R"; "R" "R" "T" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], dropobs = 20, seed = 100014)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 277.976236251267      atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.07503459156463281   atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.9372828683843416   atol=1E-5
    @test ci[end][2]                        ≈ 1.2057922163897947   atol=1E-5
    print(".")
    #15
    #TRTR/RTRT/TRRT/RTTR
    rds = CSV.File(path*"/csv/rds15.csv") |> DataFrame
    #rds = ReplicateBE.randrbeds(;n=48, sequence=[1,2,3,4], design = ["T" "R" "T" "R"; "R" "T" "R" "T" ; "T" "R" "R" "T"; "R" "T" "T" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], dropobs = 20, seed = 10015)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 331.77741574 atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.08337      atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.9431875747076479   atol=1E-5
    @test ci[end][2]                        ≈ 1.2478691227393788   atol=1E-5
    print(".")
    #16
    #TRRT/RTTR/TTRR/RRTT
    rds = CSV.File(path*"/csv/rds16.csv") |> DataFrame
    #rds = ReplicateBE.randrbeds(;n=48, sequence=[1,2,3,4], design = ["T" "R" "R" "T"; "R" "T" "T" "R" ; "T" "T" "R" "R"; "R" "R" "T" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], dropobs = 20, seed = 10016)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 285.69277340 atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.07008      atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.9088567332123243   atol=1E-5
    @test ci[end][2]                        ≈ 1.1476219506765315   atol=1E-5
    print(".")
    #17
    #TRTR/RTRT/TTRR/RRTT
    rds = CSV.File(path*"/csv/rds17.csv") |> DataFrame
    #rds = ReplicateBE.randrbeds(;n=48, sequence=[1,2,3,4], design = ["T" "R" "T" "R"; "R" "T" "R" "T" ; "T" "T" "R" "R"; "R" "R" "T" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], dropobs = 20, seed = 10017)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 292.49505051 atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.07570    atol=1E-5
    #@test ReplicateBE.coef(be)[end]        ≈ 0.05136    atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.92723952053598   atol=1E-5
    @test ci[end][2]                        ≈ 1.1951418270978003   atol=1E-5
    print(".")
    #18
    #TRT/RTR
    rds = CSV.File(path*"/csv/rds18.csv") |> DataFrame
    #rds = ReplicateBE.randrbeds(;n=48, sequence=[1,2], design = ["T" "R" "T"; "R" "T" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0], dropobs = 20, seed = 10018)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 237.09185442 atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.09209 atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.878355738329916   atol=1E-5
    @test ci[end][2]                        ≈ 1.1975977027952331   atol=1E-5
    #DF contain 34
    print(".")
    #19
    #TRR/RTT
    rds = CSV.File(path*"/csv/rds19.csv") |> DataFrame
    #rds = ReplicateBE.randrbeds(;n=48, sequence=[1,2], design = ["T" "R" "R"; "R" "T" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0], dropobs = 20, seed = 10019)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 255.99536281 atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.10125043978567916 atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.7755005422832147             atol=1E-5
    @test ci[end][2]                        ≈ 1.0916738532751222             atol=1E-5
    print(".")
    #20
    #TR/RT/TT/RR
    #SPSS REML 151.783195849874
    #SPSS SE 0.214518365508227
    #SPSS DF 10.0167597858775
    rds = CSV.File(path*"/csv/rds20.csv") |> DataFrame
    #rds = ReplicateBE.randrbeds(;n=48, sequence=[1,2,3,4], design = ["T" "R"; "R" "T"; "T" "T"; "R" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0, 0.0], periodcoef = [0.0, 0.0], formcoef = [0.0, 0.0], dropobs = 20, seed = 10020)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    @test ReplicateBE.reml2(be)             ≈ 151.78319585 atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.2145183655081034 atol=1E-5
    ci = confint(be, 0.1, expci = true)
    @test ci[end][1]                        ≈ 0.7865283528654503             atol=1E-5
    @test ci[end][2]                        ≈ 1.6924035882963178             atol=1E-5
    #DF contain 20
    print(".")
    #21
    #TRR/RTR/RRT
    rds = CSV.File(path*"/csv/rds21.csv") |> DataFrame
    #rds = ReplicateBE.randrbeds(;n=48, sequence=[1,2,3], design = ["T" "R" "R"; "R" "T" "R"; "R" "R" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0], dropobs = 20, seed = 10021)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    #SPSS REML 237.076723026247
    @test ReplicateBE.reml2(be)             ≈ 237.07672185 atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.09665    atol=1E-5
    ci = confint(be, 0.1, expci = true)
    #v1.0.4
    @test ci[end][1]                        ≈ 0.8296448401676422           atol=1E-4
    @test ci[end][2]                        ≈ 1.14767886003498             atol=1E-4
    #v1.0.5
    #@test ci[end][1]                        ≈ 0.8293880338983111             atol=1E-3
    #@test ci[end][2]                        ≈ 1.1480342197874598             atol=1E-3
    #1.4.0

    #DF contain 35
    print(".")
    #22
    #TRR/RTR
    rds = CSV.File(path*"/csv/rds22.csv") |> DataFrame
    #rds = ReplicateBE.randrbeds(;n=48, sequence=[1,2], design = ["T" "R" "R"; "R" "T" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0], dropobs = 20, seed = 10022)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    #SPSS REML  234.103074
    @test ReplicateBE.reml2(be)             ≈ 232.55722966 atol=1E-5
    @test ReplicateBE.stderror(be)[end]     ≈ 0.1090 atol=1E-5
    ci = confint(be, 0.1, expci = true)
    #v1.0.4
    @test ci[end][1]                        ≈ 0.90656475483436             atol=1E-5
    @test ci[end][2]                        ≈ 1.3053853980517822           atol=1E-5
    #v1.0.5
    #@test ci[end][1]                        ≈ 0.9064394831113237             atol=1E-3
    #@test ci[end][2]                        ≈ 1.3055658048846954             atol=1E-3
    #1.4.0

    #DF contain 32
    print(".")
    # Unbalanced by sequences
    #23
    #TRTR/RTRT
    rds = CSV.File(path*"/csv/rds23.csv") |> DataFrame
    #rds = ReplicateBE.randrbeds(;n=36, sequence=[1,2], design = ["T" "R" "T" "R"; "R" "T" "R" "T"], inter=[0.5, 0.4, 0.1], intra=[0.1, 0.15], intercept = 1.0, seqcoef = [1.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10023)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    ci = confint(be, 0.1, expci = true)
    @test ReplicateBE.reml2(be)             ≈ 252.0449059441563    atol=1E-5
    @test ci[end][1]                        ≈ 0.638039             atol=1E-5
    @test ci[end][2]                        ≈ 1.135006             atol=1E-5
    print(".")
    #24
    #TRT/RTR
    rds = ReplicateBE.randrbeds(;n=36, sequence=[1,2], design = ["T" "R" "T"; "R" "T" "R"], inter=[0.4, 0.3, 0.2], intra=[0.05, 0.02], intercept = 1.0, seqcoef = [0.1, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.1, 0.0], seed = 10024)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    ci = confint(be, 0.1, expci = true)
    @test ReplicateBE.reml2(be)             ≈ 140.10714908682417    atol=1E-5
    @test ci[end][1]                        ≈ 0.854985             atol=1E-5
    @test ci[end][2]                        ≈ 1.293459             atol=1E-5
    print(".")
    #25
    #TRR/RTT
    rds = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "R" "R"; "R" "T" "T"], inter=[0.5, 0.4, 0.5], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 10008)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    ci = confint(be, 0.1, expci = true)
    @test ReplicateBE.reml2(be)             ≈ 140.37671499958694    atol=1E-5
    @test ci[end][1]                        ≈ 0.8208303872086634    atol=1E-5
    @test ci[end][2]                        ≈ 1.3228068920153602    atol=1E-5
    print(".")
    # Unbalanced by sequences
    #26
    #TRTR/RTRT
    rds = ReplicateBE.randrbeds(;n=128, sequence=[1,2], design = ["T" "R" "T" "R"; "R" "T" "R" "T"], inter=[0.5, 0.4, 0.5], intra=[0.1, 0.15], intercept = 1.0, seqcoef = [1.0, 0.0], periodcoef = [0.0, 1.0, 0.0, 0.0], formcoef = [0.0, 0.2], seed = 10026)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    ci = confint(be, 0.1, expci = true)
    @test ReplicateBE.reml2(be)             ≈ 899.044671814185    atol=1E-5
    @test ci[end][1]                        ≈ 0.6189149073430663    atol=1E-5
    @test ci[end][2]                        ≈ 0.7726263704382558    atol=1E-5
    print(".")
    #27
    #TRT/RTR
    rds = ReplicateBE.randrbeds(;n=128, sequence=[1,2], design = ["T" "R" "T"; "R" "T" "R"], inter=[0.4, 0.3, 0.7], intra=[0.05, 0.2], intercept = 1.0, seqcoef = [1.1, 0.0], periodcoef = [0.0, 1.0, 0.0], formcoef = [0.2, 0.0], seed = 10027)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    ci = confint(be, 0.1, expci = true)
    @test ReplicateBE.reml2(be)             ≈ 614.3407387528628    atol=1E-5
    @test ci[end][1]                        ≈ 1.18696161597237    atol=1E-5
    @test ci[end][2]                        ≈ 1.43349025480276    atol=1E-5
    print(".")
    # Unbalanced by sequences

    #28
    #TRTR/RTRT
    rds = ReplicateBE.randrbeds(;n=512, sequence=[1,2], design = ["T" "R" "T" "R"; "R" "T" "R" "T"], inter=[0.5, 0.4, 0.7], intra=[0.1, 0.15], intercept = 1.0, seqcoef = [1.0, 0.0], periodcoef = [0.0, 1.0, 0.0, 0.0], formcoef = [0.0, 0.3], seed = 10028)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    ci = confint(be, 0.1, expci = true)
    @test ReplicateBE.reml2(be)             ≈ 3495.580039200916    atol=1E-5
    @test ci[end][1]                        ≈ 0.6975351860484861    atol=1E-5
    @test ci[end][2]                        ≈ 0.7698493312540238    atol=1E-5
    print(".")
    #29
    #TRT/RTR
    rds = ReplicateBE.randrbeds(;n=512, sequence=[1,2], design = ["T" "R" "T"; "R" "T" "R"], inter=[0.4, 0.3, 0.7], intra=[0.05, 0.2], intercept = 1.0, seqcoef = [1.1, 0.0], periodcoef = [0.0, 0.0, 1.0], formcoef = [0.4, 0.0], seed = 10029)
    be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
    ci = confint(be, 0.1, expci = true)
    @test ReplicateBE.reml2(be)             ≈ 2540.296180266908    atol=1E-5
    @test ci[end][1]                        ≈ 1.361487533935354    atol=1E-5
    @test ci[end][2]                        ≈ 1.498934137926506    atol=1E-5

    println(".]")
end

@testset "  #  Simulation                                  " begin
    #Power simulation by task
    io = IOBuffer()
    task = ReplicateBE.randrbetask(
        ;
        n = 12,
        sequence = [1, 1],
        design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
        inter = [0.05, 0.04, 0.6],
        intra = [0.02, 0.02],
        intercept = 1.0,
        seqcoef = [0.0, 0.0],
        periodcoef = [0.0, 0.0, 0.0, 0.0],
        formcoef = [0.0, log(0.8)],
        seed = 10001,
        dropobs = 2,
    )
    pow = ReplicateBE.simulation(task; io = io, num = 1007, seed = 1234, verbose = true)
    @test pow.result ≈ 0.06 atol=1E-2

    pow = ReplicateBE.simulation(task; io = io, num = 100, seed = 1234, verbose = true, rsabe = true)
    #@test pow.result ≈ 0.06 atol=1E-2

    #Custom simulation

    task = ReplicateBE.randrbetask(;n=24,
    sequence=[1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
    inter=[0.04, 0.04, 1.0], intra=[0.02, 0.02],
    intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0],
    formcoef = [0.0, log(0.8)], seed = 0, dropobs = 0)
    out = zeros(Float64, 0)
    function simfunc!(out, be)
        push!(out, ReplicateBE.theta(be)[5])
    end
    result             = ReplicateBE.simulation!(task, out, simfunc!; num = 10, seed = 10, verbose = false)
    Base.show(io, result)
    @test mean(result) ≈ 0.8904498245891423 atol=1E-2
end
