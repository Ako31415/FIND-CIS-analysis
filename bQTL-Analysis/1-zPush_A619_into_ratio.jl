using CSV, DataFrames

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
addA619("PATH_TO_A619_DATA_WW/A619_WW.csv", "PATH_TO_DATA/BindingFrequency/NPNRtoValue/WW-MOA_peak_ratio.csv", "keep")

# addA619("PATH_TO_A619_DATA_WW/A619_WW.csv", "PATH_TO_DATA/BindingFrequency/NPNRtoValue/WW-MOA_peak_ratio.csv", "not_keep")

#DS
addA619("PATH_TO_A619_DATA_DS/A619_DS.csv", "PATH_TO_DATA/BindingFrequency/NPNRtoValue/DS-MOA_peak_ratio.csv", "keep")

# addA619("PATH_TO_A619_DATA_DS/A619_DS.csv", "PATH_TO_DATA/BindingFrequency/NPNRtoValue/DS-MOA_peak_ratio.csv", "not_keep")
