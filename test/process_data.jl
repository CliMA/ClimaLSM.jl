using DelimitedFiles
using Dates
data = readdlm("/Users/katherinedeck/Desktop/ozark_site/AMF_US-MOz_FLUXNET_SUBSET_2004-2019_3-5/AMF_US-MOz_FLUXNET_SUBSET_HH_2004-2019_3-5.csv", ',')

P_indices = [Array(1:1:73)[data[1,:] .== "P_F"][1],Array(1:1:73)[data[1,:] .== "P_F_QC"][1]]
             # mm/interval
LE_indices = [Array(1:1:73)[data[1,:] .== "LE_CORR"][1],Array(1:1:73)[data[1,:] .== "LE_F_MDS_QC"][1]] # W/m^2
dates = DateTime.(string.(data[2:end,1]), dateformat"yyyymmddHHMM") 
t = dates .- dates[1]

chosen_year = [year(d) for d in dates] .== 2019

precip = data[2:end,:][chosen_year, P_indices]# seems fine re: P data
LE = data[2:end,:][chosen_year, LE_indices]
qc_flag = Float64.((LE[:,2] .==2) .+ (LE[:,2] .==3) ) .!= 0.0


time = Dates.value.(t[chosen_year]) ./ 1000 # seconds
seconds = time .- time[1]
LE_time = seconds[qc_flag]
ET = LE[qc_flag,1] ./ 2.257e6 ./ 1000.0 # /λ/ρ -> m/s
P  = (precip[:,1] ./ 1800.00 ./ 1000.0) # m/se


et_data = hcat(LE_time, ET)
p_data = hcat(seconds, P)

writedlm("/Users/katherinedeck/Desktop/ozark_site/p_data_2019.csv", p_data, ',')
writedlm("/Users/katherinedeck/Desktop/ozark_site/et_data_2019.csv", et_data, ',')
