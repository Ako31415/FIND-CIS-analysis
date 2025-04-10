using CSV, DataFrames

#miss = "NO"
#a619 = "/Data/michael/MOA/Met_final/Normal_res/Raw_data/B73.A619.WW.q0_EG80.lowres.norm_read_count.p01.2AFP.csv"
#peaks = "/Data/michael/MOA/Met_final/Normal_res/BindingFrequency/NPNRtoNA/p01_2AFP-WW-a619_peak_ratio.csv"

function addA619(a619, peaks, miss)

                                    # bind the a619 ratio files with the A619 infomration 

                                    a619 = CSV.read(a619, DataFrame; delim="\t", header=false)
                                    ratio = CSV.read(peaks, DataFrame; delim=",")


                                    a619.ID .= string.(a619.Column1, "_", a619.Column2)

                                    if miss != "keep"
                                        a619[occursin.("n.p.", a619[!,5]),5] .= "NA"
                                        a619[occursin.("n.r", a619[!,5]),5] .= "NA"
                                    else
                                        # change n.p or n.r to the values 
                                        a = a619[occursin.("n.p.", a619[!,5]),5]
                                        b = string.(SubString.(a, 5, length.(a)))
                                        a619[occursin.("n.p.", a619[!,5]),5] = b
                                        a619[occursin.("n.r", a619[!,5]),5] .= "NA"


                                    end


                                    rename!(a619, :Column5 => :Ratio_A619)
                                    ratio = outerjoin(a619[!,[7,5]], ratio, on=:ID)

                                    CSV.write(peaks, ratio, delim=",", missingstring="NA")


 end
 
#WW
addA619("/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/B73.A619.WW.2FPs.csv", "/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoValue/WW-MOA_peak_ratio.csv", "keep")

# addA619("/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/B73.A619.WW.2FPs.csv", "/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoNA/WW-MOA_peak_ratio.csv", "not_keep")

#DS
addA619("/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/B73.A619.DS.2FPs.csv", "/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoValue/DS-MOA_peak_ratio.csv", "keep")

# addA619("/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/B73.A619.DS.2FPs.csv", "/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoNA/DS-MOA_peak_ratio.csv", "not_keep")
