using DelimitedFiles
using Dates
data = readdlm("/Users/katherinedeck/Desktop/ozark_site/AMF_US-MOz_FLUXNET_SUBSET_2004-2019_3-5/AMF_US-MOz_FLUXNET_SUBSET_HH_2004-2019_3-5.csv", ',')

P_indices = [Array(1:1:73)[data[1,:] .== "P_F"][1],Array(1:1:73)[data[1,:] .== "P_F_QC"][1]]
             # mm/interval
LE_indices = [Array(1:1:73)[data[1,:] .== "LE_CORR"][1],Array(1:1:73)[data[1,:] .== "LE_F_MDS_QC"][1]] # W/m^2

soil_water_indices = [Array(1:1:73)[data[1,:] .== "SWC_F_MDS_1"][1],Array(1:1:73)[data[1,:] .== "SWC_F_MDS_1_QC"][1]] # %
soil_temp_indices = [Array(1:1:73)[data[1,:] .== "TS_F_MDS_1"][1],Array(1:1:73)[data[1,:] .== "TS_F_MDS_1_QC"][1]] # deg C

dates = DateTime.(string.(data[2:end,1]), dateformat"yyyymmddHHMM") 
t = dates .- dates[1]

# Pick a year
years = 2004:1:2019
for y in years
    chosen_y = [year(d) for d in dates] .== y
    LE_y = data[2:end,:][chosen_y, LE_indices]
    qc_flag_y = Float64.((LE_y[:,2] .==2) .+ (LE_y[:,2] .==3) ) .== 0.0
    println(y)
    println(sum(qc_flag_y)./sum(chosen_y)*100)
end

    


chosen_year = [year(d) for d in dates] .== 2005

precip = data[2:end,:][chosen_year, P_indices]# seems fine re: P data
LE = data[2:end,:][chosen_year, LE_indices]
swc = data[2:end,:][chosen_year, soil_water_indices]
swc_qc_flag = Float64.((swc[:,2] .==-9999) .+ (swc[:,2] .==3) ) .== 0.0 # all good

ts = data[2:end,:][chosen_year, soil_temp_indices]
ts_qc_flag = Float64.((ts[:,2] .==-9999) .+ (ts[:,2] .==3) ) .== 0.0 # all good

qc_flag = Float64.((LE[:,2] .==2) .+ (LE[:,2] .==3) ) .== 0.0


time = Dates.value.(t[chosen_year]) ./ 1000 # seconds
seconds = time .- time[1]
LE_time = seconds[qc_flag]
ET = LE[qc_flag,1] ./ 2.257e6 ./ 1000.0 # /λ/ρ -> m/s
P  = (precip[:,1] ./ 1800.00 ./ 1000.0) # m/se
θ = swc[:,1]
T = ts[:,1]

et_data = hcat(LE_time, ET)
otherwise_data = hcat(seconds, P, θ, T)

writedlm("/Users/katherinedeck/Desktop/ozark_site/p_data_2005.csv", otherwise_data, ',')
writedlm("/Users/katherinedeck/Desktop/ozark_site/et_data_2005.csv", et_data, ',')
