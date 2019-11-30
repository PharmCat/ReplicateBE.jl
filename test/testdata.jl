# ReplicateBE
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
# Licence: GNU General Public License v3.0

#Simple dataset
df0         = CSV.File(path*"/csv/df0.csv") |> DataFrame
df1         = CSV.File(path*"/csv/df1.csv") |> DataFrame
df2         = CSV.File(path*"/csv/df2.csv") |> DataFrame
df3         = CSV.File(path*"/csv/df3.csv") |> DataFrame
#QA 1 Bioequivalence 2x2x4, UnB, NC Dataset
df4         = CSV.File(path*"/csv/df4.csv") |> DataFrame
#Patterson SD, Jones B. Viewpoint: observations on scaled average bioequivalence. Pharm Stat. 2012; 11(1): 1–7. doi:10.1002/pst.498
df5         = CSV.File(path*"/csv/df5.csv") |> DataFrame
#Shumaker RC, Metzler CM. The Phenytoin Trial is a Case Study of ‘Individual’ Bioequivalence. Drug Inf J. 1998; 32(4): 1063–72
df6         = CSV.File(path*"/csv/df6.csv") |> DataFrame
