# ReplicateBE
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
# Licence: GNU General Public License v3.0

println("df0...")
#Simple dataset
df0         = CSV.File(path*"/csv/df0.csv") |> DataFrame
println("df1...")
#Scott D. Patterson, Byron Jones, Bioequivalence and Statistics in Clinical Pharmacology 2nd Edition - ISBN 9781466585201, example 4.3, t 4.29, p. 104
df1         = CSV.File(path*"/csv/df1.csv") |> DataFrame
println("df2...")
#Shein-Chung Chow, Jen-pei Liu, Design and Analysis of Bioavailability and Bioequivalence Studies, 3rd Edition, ISBN 9781584886686, example 9.4.1, t. 9.4.3, p. 284
df2         = CSV.File(path*"/csv/df2.csv") |> DataFrame
println("df3...")
#Scott D. Patterson, Byron Jones, Bioequivalence and Statistics in Clinical Pharmacology 2nd Edition - ISBN 9781466585201, example 4.4, t. 4.30, p. 105
df3         = CSV.File(path*"/csv/df3.csv") |> DataFrame
println("df4...")
#QA 1 Bioequivalence 2x2x4, UnB, NC Dataset
df4         = CSV.File(path*"/csv/df4.csv") |> DataFrame
println("df5...")
#Patterson SD, Jones B. Viewpoint: observations on scaled average bioequivalence. Pharm Stat. 2012; 11(1): 1–7. doi:10.1002/pst.498
df5         = CSV.File(path*"/csv/df5.csv") |> DataFrame
println("df6...")
#Shumaker RC, Metzler CM. The Phenytoin Trial is a Case Study of ‘Individual’ Bioequivalence. Drug Inf J. 1998; 32(4): 1063–72
df6         = CSV.File(path*"/csv/df6.csv") |> DataFrame
