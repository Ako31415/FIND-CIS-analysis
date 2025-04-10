using CSV, DataFrames

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
addA619("PATH_TO_A619_DATA_WW/A619_WW.csv", "PATH_TO_DATA/BindingFrequency/NPNRtoValue/WW-ReadDepth_ratio.csv")

# addA619("PATH_TO_A619_DATA_WW/A619_WW.csv", "PATH_TO_DATA/BindingFrequency/NPNRtoNA/WW-ReadDepth_ratio.csv")


#DS
addA619("PATH_TO_A619_DATA_DS/A619_DS.csv", "PATH_TO_DATA/BindingFrequency/NPNRtoValue/DS-ReadDepth_ratio.csv")

# addA619("PATH_TO_A619_DATA_DS/A619_DS.csv", "PATH_TO_DATA/BindingFrequency/NPNRtoNA/DS-ReadDepth_ratio.csv")

