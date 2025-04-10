## read in the chromosome file from external file


aom = ARGS[1] 
me = ARGS[2] # which chromsome to load
con = ARGS[3]
rd = ARGS[4] #read depth to load
oneg = ARGS[5]  #genotype to load


path="PATH_TO_DATA"

meth = ["$path/BindingFrequency/NPNRtoValue/splitted/$(con)-CG_ratio_$(me)_file.csv", 
"$path/BindingFrequency/NPNRtoValue/splitted/$(con)-CHG_ratio_$(me)_file.csv", 
"$path/BindingFrequency/NPNRtoValue/splitted/$(con)-CHH_ratio_$(me)_file.csv"]

println("list files")
println("$(aom), $(oneg), $(meth)")

using CSV, DataFrames, GLM, StatsModels, ProgressMeter, Statistics, StatsBase, RCall, Suppressor


moa = CSV.read(aom, DataFrame; missingstring="NA")
MOA = nrow(moa)
# sort by position 
moa[!,:Chr] .= ""
moa[!,:Pos] .= ""
for s in 1:MOA 
               moa[s,ncol(moa)-1:ncol(moa)] = split(moa[s,1], "_")
end
moa.Pos = parse.(Int64, moa.Pos)
sort!(moa, [:Chr, :Pos])

geno = CSV.read(oneg, DataFrame; missingstring="NA")
# remove A619 for now 
#geno = geno[!,Not(:A619)]

# sort the moa file in the same column order as the geno file
nam = Symbol.(String.(SubString.(names(moa)[Not([1,ncol(moa)-1,ncol(moa)])], 7, length.(names(moa)[Not([1,ncol(moa)-1,ncol(moa)])]))))
nam = vcat(Symbol("ID"), nam)
moa = moa[!,Not(ncol(moa)-1:ncol(moa))]

rename!(moa, nam)

moa = hcat(moa[!,1], moa[!,Symbol.(names(geno)[Not(1:2)])])
rename!(moa, :x1 => :ID)

# initialize an array to store the lm output in
sink6 = DataFrame(MOA_window = String[], SNP = String[], log10ProbGeno = Float64[], numberOFsnps = Int64[], ExpVariance = Float64[])
sink5 = DataFrame(MOA_window = String[], SNP = String[], log10ProbGeno = Float64[], log10ProbMethy = Union{Missing, Float64}[], numberOFsnps = Int64[], ExpVariance = Float64[])
sink4 = DataFrame(MOA_window = String[], SNP = String[], log10ProbGeno = Float64[], log10ProbMethy = Union{Missing, Float64}[], numberOFsnps = Int64[], ExpVariance = Float64[])
sink3 = DataFrame(MOA_window = String[], SNP = String[], log10ProbGeno = Float64[], log10ProbMethy = Union{Missing, Float64}[], numberOFsnps = Int64[], ExpVariance = Float64[])
sink2 = DataFrame(MOA_window = String[], SNP = String[], log10ProbGeno = Float64[], log10pCG = Union{Missing, Float64}[], log10pCHG = Union{Missing, Float64}[], log10pCHH = Union{Missing, Float64}[],numberOFsnps = Int64[], ExpVariance = Float64[])

sink66 = DataFrame(MOA_window = String[], SNP = String[], log10ProbGeno = Float64[], numberOFsnps = Int64[], ExpVariance = Float64[])
sink55 = DataFrame(MOA_window = String[], SNP = String[], log10ProbGeno = Float64[], log10ProbMethy = Union{Missing, Float64}[], numberOFsnps = Int64[], ExpVariance = Float64[])
sink44 = DataFrame(MOA_window = String[], SNP = String[], log10ProbGeno = Float64[], log10ProbMethy = Union{Missing, Float64}[], numberOFsnps = Int64[], ExpVariance = Float64[])
sink33 = DataFrame(MOA_window = String[], SNP = String[], log10ProbGeno = Float64[], log10ProbMethy = Union{Missing, Float64}[], numberOFsnps = Int64[], ExpVariance = Float64[])
sink22 = DataFrame(MOA_window = String[], SNP = String[], log10ProbGeno = Float64[], log10pCG = Union{Missing, Float64}[], log10pCHG = Union{Missing, Float64}[], log10pCHH = Union{Missing, Float64}[],numberOFsnps = Int64[], ExpVariance = Float64[])
sink77 = DataFrame(MOA_window = String[], SNP = String[], log10pCG = Union{Missing, Float64}[], log10pCHG = Union{Missing, Float64}[], log10pCHH = Union{Missing, Float64}[],numberOFsnps = Int64[], ExpVariance = Float64[])


Results = Dict(:ALLmm => sink22, :CGmm => sink33, :CHGmm => sink44, :CHHmm => sink55, :SNPonlymm => sink66, :METHonly => sink77)



chr = geno[!,1]
Position = geno[!,2]

#oFN = string(split(split(aom, "_")[3], "/")[2], "_Chr", split(aom, "_")[6])

## load the methylation info 
METH = Dict()

METH[:CG] = CSV.read(meth[1], DataFrame; delim=",",  missingstring="NA")
METH[:CHH] = CSV.read(meth[3], DataFrame; delim=",",  missingstring="NA")
METH[:CHG] = CSV.read(meth[2], DataFrame; delim=",",  missingstring="NA")

# change column names and order 
for yxc in collect(keys(METH))
                                newnames = ["ID"]
                                for i in names(METH[yxc])[2:length(names(METH[yxc]))]
                                                                            push!(newnames, split(i, "_")[2])
                                    end

                                rename!(METH[yxc], newnames)


                                ## add the imputed information of A619 methylation data here
                                METH[yxc][!,:A619] .= missing

                                ## change order - until here is correct!!
                                METH[yxc] = hcat(METH[yxc][!,1], METH[yxc][!,Symbol.(names(geno)[Not(1:2)])])
                                rename!(METH[yxc], :x1 => :ID)
                        end

# make sure moa and the methylation information do all have the same number of rows
METH[:CG] = leftjoin(moa[!,[:ID]], METH[:CG], on=:ID) # does not do anything - for whatever reason?! should sort meth by moa, but doesnot do the job
METH[:CHH] = leftjoin(moa[!,[:ID]], METH[:CHH], on=:ID)
METH[:CHG] = leftjoin(moa[!,[:ID]], METH[:CHG], on=:ID)

METH[:CG] = METH[:CG][indexin(moa.ID, METH[:CG].ID),:] # this here performs the correct sorting.
METH[:CHH] = METH[:CHH][indexin(moa.ID, METH[:CHH].ID),:]
METH[:CHG] = METH[:CHG][indexin(moa.ID, METH[:CHG].ID),:]

# read the readdepth information and join it with the moa Data
RD =  CSV.read(rd, DataFrame; delim=",",  missingstring="NA")

RD = leftjoin(moa[!,[:ID]], RD, on=:ID)

#for KEYS in collect(keys(METH)
                               ## remove the loci with no methylation informatio
                               #test = Matrix(METH[KEYS][!,2:ncol(METH[KEYS])]
#
                               ##remove = Bool[
                               ##for i in 1:size(test,1
                               ##                        append!(remove, all(ismissing.(test[i,:]))
                               ##    en
#
                               ##test = test[.!(remove),:
#
                               ##moa = moa[.!(remove),:
#
                               ## impute the missing methylation values with the mean of the other genotypes for all missing informatio
                               #for i in 1:size(test,1
                                                       #test[i,ismissing.(test[i,:])] .= mean(skipmissing(test[i,:])
                                       #en
#
                               ## send it back to the METH dict entry
                               ##test = DataFrame(hcat(METH[KEYS][.!(remove),:ID], test), :auto
                               #test = DataFrame(hcat(METH[KEYS][!,:ID], test), :auto
                               #rename!(test, names(METH[KEYS])
                               #METH[KEYS] = tes
#
           #en
#
# findfirst(x -> x == "B73-chr1_2089898",  moa.ID) # the row Julia wanted to know about 
# findfirst(x -> occursin("B73-chr1_2089898", x), moa.ID)
# itt = 392
# 8297

@showprogress for itt in 1:MOA
                                #println(itt)
                                #if isinteger(itt/10000)
                                #    R"""
                                #    colnames(sink) = c("MOA", "Best_SNP", "logPscore")
                                #    write.table(sink, "/home/michael-uni/Dokumente/localserver//BestSNPs_running.csv", row.names=F, append=TRUE, col.names=F)
                                #    """
                                #end
                    #try
                                moaC = split(moa[itt,:ID],"_")[1]
                                moaP = parse(Int64, split(moa[itt,:ID],"_")[2])
                                candidates = geno[.&(Position .== moaP, chr.== moaC),:]

                            if !isempty(candidates)
                                canM = DataFrame(hcat(transpose(Matrix(moa[[itt], Not([1])])), transpose(Matrix(METH[:CG][[itt], Not([1])])), transpose(Matrix(METH[:CHH][[itt], Not([1])])), 
                                             transpose(Matrix(METH[:CHG][[itt], Not([1])])), transpose(Matrix(RD[[itt], Not([1])])), transpose(Matrix(candidates[!,Not(1,2)]))),:auto)
                                rename!(canM, :x1 => :MOA)
                                canM = hcat(names(geno)[Not(1,2)], canM)
                                rename!(canM, :x1 => :Genotype)
                                rename!(canM, :x2 => :CG)
                                rename!(canM, :x3 => :CHH)
                                rename!(canM, :x4 => :CHG)
                                rename!(canM, :x5 => :ReadDepth)

                                canM[ismissing.(canM.ReadDepth),:ReadDepth] .= mean(skipmissing(canM.ReadDepth))
                                canM.ReadDepth = Int64.(round.(canM.ReadDepth))
                                # impute the moa ratio values
                                if !(all(ismissing.(canM.MOA))) # for some, there is no MOa binding frequency infomration delived - and this will relsut in an error
                                                                canM[ismissing.(canM.MOA),:MOA] .= mean(skipmissing(canM.MOA))
                                                                canM[!,:MOA] .= Float64.(canM.MOA)
                                else
                                    continue
                                end
                            
                                # convert any to Float or int values
                                #for i in 3:5; canM[!,i] = Float64.(skipmissing(canM[!,i])); end
                                #for i in 6:ncol(canM); canM[!,i] = Int64.(skipmissing(canM[!,i])); end
                            # at the moment, the p value is only calculated if all three methlation classes have a sufficent number of ratios to calculate a linear model

                            
                                # calculate a simple linear model for each marker
                             #=   PvaluesLM = DataFrame(snpSingle = Float64[], snpAll = Float64[], snpCHG = Float64[], snpCG = Float64[], snpCHH = Float64[], MethAllCG = Union{Missing, Float64}[], MethAllCHH = Union{Missing, Float64}[], 
                                                        MethAllCHG = Union{Missing, Float64}[], MethCHG = Union{Missing, Float64}[], MethCG = Union{Missing, Float64}[], MethCHH = Union{Missing, Float64}[])
                                eVarLM = DataFrame(all4 = Float64[], CHG = Float64[], CG = Float64[], CHH = Float64[], SNP = Float64[])
                            =#
                                PvaluesMM = DataFrame(snpSingle = Float64[], snpAll = Float64[], snpCHG = Float64[], snpCG = Float64[], 
                                                      snpCHH = Float64[], MethAllCG = Union{Missing, Float64}[], MethAllCHH = Union{Missing, Float64}[], 
                                                        MethAllCHG = Union{Missing, Float64}[], MethCHG = Union{Missing, Float64}[], 
                                                        MethCG = Union{Missing, Float64}[], MethCHH = Union{Missing, Float64}[], 
                                                        methONLYCG = Union{Missing, Float64}[], methONLYCHG = Union{Missing, Float64}[], methONLYCHH = Union{Missing, Float64}[])

                                eVarMM = DataFrame(all4 = Float64[], CHG = Float64[], CG = Float64[], CHH = Float64[], SNP = Float64[], METH = Float64[])
                                
                               

                          
                            allequal(x) = all(y->y==x[1],x) # check if all entries are equal

                        try  # Mixed model 

                             candidates = geno[.&(Position .== moaP, chr.== moaC),:]

                            for i in names(canM)[7:length(names(canM))] 

                                #fm = @formula(yield ~ 1 + (1|batch))
                                #fm1 = fit(MixedModel, fm, dyestuff)

                                                                       # fit(MixedModel, @formula(MOA ~ CHG + CG + CHH + (1 | ReadDepth), canM))






                                #### the R veraion - 2.585 seconds
                                                                        @suppress begin
                                                                        @rput i
                                                                                        # remove the missing genotypes for the itteration we are working on 
                                                                                        canM2 = canM[.!(ismissing.(canM[:,Symbol(i)])),:]
                                                                                        
                                                                                        if allequal(canM2[!,i])
                                                                                            
                                                                                            candidates = candidates[Not(parse(Int64, SubString(i, 2, 3)) - 7),:]
                                                                                            continue 

                                                                                        end
                                                                        
                                                                        @rput canM

                                                                                ##### Mixed Model 
                                                                                R"""
                                                                                #library(lme4)
                                                                                #library(lmerTest)
                                                                                # get the p lvaue of the entire model
                                                                                a1 = summary(lm(MOA ~ canM[,i] + CHG + CG + CHH , data=canM))$coefficients
                                                                                a2 = summary(lm(MOA ~ canM[,i] + CHG , data=canM))$coefficients
                                                                                a3 = summary(lm(MOA ~ canM[,i] + CG , data=canM))$coefficients
                                                                                a4 = summary(lm(MOA ~ canM[,i] + CHH , data=canM))$coefficients
                                                                                a5 = summary(lm(MOA ~ canM[,i] , data=canM))$coefficients
                                                                                a6 = summary(lm(MOA ~ CHG + CG + CHH , data=canM))$coefficients 
                                                                                # when you are sitting on a single SNP, this will run as often as the other a1 to a5 - other than that, it might print multiple times the same result in the result file 

                                                                                pvals = c(a5[2,4], a1[2,4], a2[2,4], a3[2,4], a4[2,4])

                                                                                  # verify the value is not empty 
                                                                                  if (length(a1[row.names(a1) == "CG",4])>0){
                                                                                                                                        pvals = c(pvals, a1[row.names(a1) == "CG",4])
                                                                                                }else {pvals = c(pvals, NA)}

                                                                                            if (length(a1[row.names(a1) == "CHH",4])>0){
                                                                                                                                        pvals = c(pvals, a1[row.names(a1) == "CHH",4])
                                                                                                }else {pvals = c(pvals, NA)}
                                                                                                
                                                                                            if (length(a1[row.names(a1) == "CHG",4])>0){
                                                                                                                                        pvals = c(pvals, a1[row.names(a1) == "CHG",4])
                                                                                                }else {pvals = c(pvals, NA)}
                                                                                                
                                                                                            if (length(a2[row.names(a2) == "CHG",4])>0){
                                                                                                                                        pvals = c(pvals, a2[row.names(a2) == "CHG",4])
                                                                                                }else {pvals = c(pvals, NA)}

                                                                                            if (length(a3[row.names(a3) == "CG",4])>0){
                                                                                                                                        pvals = c(pvals, a3[row.names(a3) == "CG",4])
                                                                                                }else {pvals = c(pvals, NA)}

                                                                                            if (length(a4[row.names(a4) == "CHH",4])>0){
                                                                                                                                            pvals = c(pvals, a4[row.names(a4) == "CHH",4])
                                                                                                }else {pvals = c(pvals, NA)}

                                                                                ## adding the pvalues for the model a6 with the methylation only 
                                                                                if (length(a6[row.names(a6) == "CG",4])>0){
                                                                                                                            pvals = c(pvals, a6[row.names(a6) == "CG",4])
                                                                                                                            }else {pvals = c(pvals, NA)}

                                                                                if (length(a6[row.names(a6) == "CHG",4])>0){
                                                                                                                            pvals = c(pvals, a6[row.names(a6) == "CHG",4])
                                                                                                                            }else {pvals = c(pvals, NA)}
                                                                                                                            
                                                                                if (length(a6[row.names(a6) == "CHH",4])>0){
                                                                                                                            pvals = c(pvals, a6[row.names(a6) == "CHH",4])
                                                                                                                            }else {pvals = c(pvals, NA)}  

                                                                                ###
                                                                                pvals = -log10(pvals)

                                                                                        # calculate the explained variance  
                                                                                        a1 = anova(lm(MOA ~ canM[,i] + CHG + CG + CHH, data=canM))
                                                                                        a2 = anova(lm(MOA ~ canM[,i] + CHG, data=canM))
                                                                                        a3 = anova(lm(MOA ~ canM[,i] + CG , data=canM))
                                                                                        a4 = anova(lm(MOA ~ canM[,i] + CHH, data=canM))
                                                                                        a5 = anova(lm(MOA ~ canM[,i] , data=canM))
                                                                                        a6 = anova(lm(MOA ~ CHG + CG + CHH , data=canM))

                                                                                        A1 = (sum(a1[,2]) - a1[nrow(a1),2]) / sum(a1[,2])
                                                                                        A2 = (sum(a2[,2]) - a2[nrow(a2),2]) / sum(a2[,2])
                                                                                        A3 = (sum(a3[,2]) - a3[nrow(a3),2]) / sum(a3[,2])
                                                                                        A4 = (sum(a4[,2]) - a4[nrow(a4),2]) / sum(a4[,2])
                                                                                        A5 = (sum(a5[,2]) - a5[nrow(a5),2]) / sum(a5[,2])
                                                                                        A6 = (sum(a6[,2]) - a6[nrow(a6),2]) / sum(a6[,2])


                                                                                        expVariance = c(A1, A2, A3, A4, A5, A6)
                                                                                     
                                                                                """
                                                                                @rget pvals expVariance

                                                                                push!(PvaluesMM, pvals)
                                                                                push!(eVarMM, expVariance)
                                                                end # suppress

                            end # for loop mixed models 
                            # concet the candidates with the explained variance, to simple extract the correct row 
                            candidates = hcat(candidates, eVarMM)

                            # save the 4 different models in 4 different data frames 
                            genoSingle = PvaluesMM.snpSingle[.!(isnan.(PvaluesMM.snpSingle))]
                            genolm4 = PvaluesMM.snpAll[.!(isnan.(PvaluesMM.snpAll))]
                            genoCG = PvaluesMM.snpCG[.!(isnan.(PvaluesMM.snpCG))]
                            genoCHG = PvaluesMM.snpCHG[.!(isnan.(PvaluesMM.snpCHG))]
                            genoCHH = PvaluesMM.snpCHH[.!(isnan.(PvaluesMM.snpCHH))]

                            # model with all methylation cofactors
                            pvalmax = maximum(genolm4) # get the best p value
                            genolm4 = hcat(candidates.POS, genolm4)
                            genolm4 = genolm4[genolm4[:,2] .== pvalmax,:]
                            genolm4 = hcat(genolm4, abs.(genolm4[:,1] .- moaP))
                            PLM = genolm4[genolm4[:,3] .== minimum(genolm4[:,3]),1] # get the position of the most SNP with the highest logP value 
                            SNP = string(candidates[1, 1], "_", Int64(PLM[1])) # get the marker name 
                            # get the corresponding methylation values 
                            genolm4 = PvaluesMM.snpAll[.!(isnan.(PvaluesMM.snpAll))]
                            dpmet = PvaluesMM[findfirst(==(pvalmax), genolm4),[6,8,7]]
                            ev = eVarMM[findfirst(==(pvalmax), genolm4), 1]
                            push!(Results[:ALLmm], [moa.ID[itt], SNP, pvalmax, dpmet[1], dpmet[2], dpmet[3], ncol(canM)-6, ev])

                            # model with cg only
                            pvalmax = maximum(genoCG) # get the best p value
                            genoCG = hcat(candidates.POS, genoCG)
                            genoCG = genoCG[genoCG[:,2] .== pvalmax,:]
                            genoCG = hcat(genoCG, abs.(genoCG[:,1] .- moaP))
                            PLM = genoCG[genoCG[:,3] .== minimum(genoCG[:,3]),1] # get the position of the most SNP with the highest logP value 
                            SNP = string(candidates[1, 1], "_", Int64(PLM[1])) # get the marker name 
                            # methylation values 
                            genoCG = PvaluesMM.snpCG[.!(isnan.(PvaluesMM.snpCG))]
                            dpmet = PvaluesMM[findfirst(==(pvalmax), genoCG),10]
                            ev = eVarMM[findfirst(==(pvalmax), genoCG), 3]
                            push!(Results[:CGmm], [moa.ID[itt], SNP, pvalmax, dpmet, ncol(canM)-6, ev])

                            # model with chg only
                            pvalmax = maximum(genoCHG) # get the best p value
                            genoCHG = hcat(candidates.POS, genoCHG)
                            genoCHG = genoCHG[genoCHG[:,2] .== pvalmax,:]
                            genoCHG = hcat(genoCHG, abs.(genoCHG[:,1] .- moaP))
                            PLM = genoCHG[genoCHG[:,3] .== minimum(genoCHG[:,3]),1] # get the position of the most SNP with the highest logP value 
                            SNP = string(candidates[1, 1], "_", Int64(PLM[1])) # get the marker name 
                            # get the corresponding methylation values 
                            genoCHG = PvaluesMM.snpCHG[.!(isnan.(PvaluesMM.snpCHG))]
                            dpmet = PvaluesMM[findfirst(==(pvalmax), genoCHG),9]
                            ev = eVarMM[findfirst(==(pvalmax), genoCHG), 2]
                            push!(Results[:CHGmm], [moa.ID[itt], SNP, pvalmax, dpmet, ncol(canM)-6, ev])

                            # model with chh only
                            pvalmax = maximum(genoCHH) # get the best p value
                            genoCHH = hcat(candidates.POS, genoCHH)
                            genoCHH = genoCHH[genoCHH[:,2] .== pvalmax,:]
                            genoCHH = hcat(genoCHH, abs.(genoCHH[:,1] .- moaP))
                            PLM = genoCHH[genoCHH[:,3] .== minimum(genoCHH[:,3]),1] # get the position of the most SNP with the highest logP value 
                            SNP = string(candidates[1, 1], "_", Int64(PLM[1])) # get the marker name 
                            # get the corresponding methylation values 
                            genoCHH = PvaluesMM.snpCHH[.!(isnan.(PvaluesMM.snpCHH))]
                            dpmet = PvaluesMM[findfirst(==(pvalmax), genoCHH),11]
                            ev = eVarMM[findfirst(==(pvalmax), genoCHH), 4]
                            push!(Results[:CHHmm], [moa.ID[itt], SNP, pvalmax, dpmet, ncol(canM)-6, ev])

                            # model with marker only
                            pvalmax = maximum(genoSingle) # get the best p value
                            genoSingle = hcat(candidates.POS, genoSingle)
                            genoSingle = genoSingle[genoSingle[:,2] .== pvalmax,:]
                            genoSingle = hcat(genoSingle, abs.(genoSingle[:,1] .- moaP))
                            PLM = genoSingle[genoSingle[:,3] .== minimum(genoSingle[:,3]),1] # get the position of the most SNP with the highest logP value 
                            SNP = string(candidates[1, 1], "_", Int64(PLM[1])) # get the marker name 
                            genoSingle = PvaluesMM.snpSingle[.!(isnan.(PvaluesMM.snpSingle))]
                            ev = eVarMM[findfirst(==(pvalmax), genoSingle), 5]
                            push!(Results[:SNPonlymm], [moa.ID[itt], SNP, pvalmax, ncol(canM)-6, ev])

                            # model with methylation only
                            pvalmax = maximum(genolm4) # get the best p value 
                            dpmet = PvaluesMM[findfirst(==(pvalmax), genolm4),[12,13,14]] # using genolm4 as helper to select the correct line - only necessary (will change something in case a window with multiple SNPs is tested)
                            ev = eVarMM[findfirst(==(pvalmax), genolm4), 6]
                            push!(Results[:METHonly], [moa.ID[itt], SNP, dpmet[1], dpmet[2], dpmet[3], ncol(canM)-6, ev])


                        catch
                            println("MOA column $itt not working on Chr $me for MM")
                           continue
                     end
                        

                          end # candidates

                          
           # catch
           #     continue

        #end # most outer catch for Int64(missing)
            
 end # outer for loop 



# specify folder
l1 = split(split(aom, "/")[length(split(aom, "/"))], "MOA")[1] 
l2 = string(split(aom, "_")[9])# chr
disk = "PATH_TO_DATA/Results/NPNRtoValue/snponly/$con/"

for m in string.(collect(keys(Results)))
                                LOC = string(disk, con,"-" ,"MOA_", m, "_", me, "_", ".csv")
                                CSV.write(LOC, Results[Symbol(m)])
                            end
