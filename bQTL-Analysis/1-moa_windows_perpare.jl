using CSV, DataFrames, ProgressMeter

## create the moa peak data set for mapping purpose 
# first column - 1 position as 'B73-chr1_2025073'

# specify the #1 folder and the #2 sub name of the files to look up, and #3 if WW or DS
@show ARGS
if ARGS == String[]
    println(" No arguments specified \n first and only argument is the entire path to the file including the filename, e.g. ~/Documents/Variants.vvcf")

end

folder = ARGS[1]
#set = ARGS[2]
ENV = ARGS[2]
miss = ARGS[3]

#folder = "/home/michael-uni/Dokumente/qggp/michael/MOA/MOA_methylation/RawData/DiffMet"
#set = "2AFPs" # not required
#ENV = "WW"
# miss = "keep"

s3 = readdir(folder)

# create an empty matrix of the size length(s1) and norw(x)
ratio  = Dict()
CG = Dict() 
CHG = Dict()
CHH = Dict()
ReadDepth = Dict()

@showprogress for i in s3

            # get the name of the hybrid from the name 'i'
            hyNam = split(i, ".")[2]

            moa = CSV.read("$folder/$i", DataFrame; delim="\t", header=false, missingstring="NA")
            #set the names Chr    Pos    NAM-Pos    Genotype    Bindingfrequency Readdepth-B73+NAM      Methylation-DifferenceCG Methylation-DifferenceCHG    Methylation-DifferenceCHH 
            rename!(moa, ["Chr", "Pos", "$hyNam-coordinates", "GT", "Ratio_$(hyNam)", "Readdepth_$(hyNam)", "CG_$(hyNam)", "CHG_$(hyNam)", "CHH_$(hyNam)",])

            # convert missing to 1
            if miss != "keep"
                                moa[occursin.("n.p.", moa[!,5]),5] .= "NA"
                                moa[occursin.("n.r", moa[!,5]),5] .= "NA"
            else
                                # change n.p or n.r to the values 
                                a = moa[occursin.("n.p.", moa[!,5]),5]
                                b = string.(SubString.(a, 5, length.(a)))
                                moa[occursin.("n.p.", moa[!,5]),5] = b
                                moa[occursin.("n.r", moa[!,5]),5] .= "NA"


            end

            # create a position id 
            moa[!,:ID] .= string.(moa.Chr, "_", moa.Pos)
            ratio[Symbol(hyNam)] = moa[!,[10,5]]
            CG[Symbol(hyNam)] = moa[!,[10,7]]
            CHG[Symbol(hyNam)] = moa[!,[10,8]]
            CHH[Symbol(hyNam)] = moa[!,[10,9]] 
            ReadDepth[Symbol(hyNam)] = moa[!,[10,6]] 

    end 

# merge the dicts in a single table 
ratiO = ratio[collect(keys(ratio))[1]][!,[1,2]]
cg = CG[collect(keys(CG))[1]][!,[1,2]]
chg = CHG[collect(keys(CHG))[1]][!,[1,2]]
chh = CHH[collect(keys(CHH))[1]][!,[1,2]]
readdepth = ReadDepth[collect(keys(ReadDepth))[1]][!,[1,2]]


for i in collect(keys(ratio))
                            if i !=  collect(keys(ratio))[1]
                                                            global ratiO = outerjoin(ratiO, ratio[i], on=:ID)
                                                            global cg = outerjoin(cg, CG[i], on=:ID)
                                                            global chg = outerjoin(chg, CHG[i], on=:ID)
                                                            global chh = outerjoin(chh, CHH[i], on=:ID)
                                                            global readdepth = outerjoin(readdepth, ReadDepth[i], on=:ID)


                            end
    end

if miss == "keep"
                    loc = "NPNRtoValue"
else
                    loc = "NPNRtoNA"
end


set = split(folder, "/")[length(split(folder, "/"))]

# write to file 
CSV.write("/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/$ENV-MOA_peak_ratio.csv", ratiO, delim=",", missingstring="NA")
CSV.write("/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/$ENV-CG_ratio.csv", cg, delim=",", missingstring="NA")
CSV.write("/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/$ENV-CHG_ratio.csv", chg, delim=",", missingstring="NA")
CSV.write("/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/$ENV-CHH_ratio.csv", chh, delim=",", missingstring="NA")
CSV.write("/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/$ENV-ReadDepth_ratio.csv", readdepth, delim=",", missingstring="NA")




# check which genotypes are present
#genotypes = ["A619","B97","CML277","CML322","CML333","CML69","HP301","IL14H","Ki11","Ki3","Ky21","M162W","Mo17","Mo18W","Ms71","NC358","Oh43","Tx303"]
