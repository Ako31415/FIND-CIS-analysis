using CSV, DataFrames

#miss = "NO"
#a619 = "/Data/michael/MOA/Met_final/Normal_res/Raw_data/B73.A619.WW.q0_EG80.lowres.norm_read_count.p01.2AFP.csv"
#peaks = "/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/MOA_QTL/BindingFrequency/NPNRtoNA/p01_2AFP-WW-a619_peak_ratio.csv"

function addA619(a619, peaks)

                                    # bind the a619 ratio files with the A619 infomration 

                                    a619 = CSV.read(a619, DataFrame; delim="\t", header=false)
                                    ratio = CSV.read(peaks, DataFrame; delim=",")


                                    a619.ID .= string.(a619.Column1, "_", a619.Column2)


                                    rename!(a619, :Column6 => :Readdepth_A619)

                                    ratio = outerjoin(a619[!,[7,6]], ratio, on=:ID)

                                    CSV.write(peaks, ratio, delim=",", missingstring="NA")


 end
 
#WW
addA619("/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/B73.A619.WW.2FPs.csv", "/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoValue/WW-ReadDepth_ratio.csv")

# addA619("/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/B73.A619.WW.2FPs.csv", "/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoNA/WW-ReadDepth_ratio.csv")


#DS
addA619("/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/B73.A619.DS.2FPs.csv", "/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoValue/DS-ReadDepth_ratio.csv")

# addA619("/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/B73.A619.DS.2FPs.csv", "/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoNA/DS-ReadDepth_ratio.csv")

