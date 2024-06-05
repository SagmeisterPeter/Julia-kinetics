struct ExperimentalData
    time_hitec
    flowrate
    flowrate_total 
    Tset
    C_inlet
    time_ftir
    c_ftir
    tmax
end
#using XLSX
#filename = "C:/Users/peter.sagmeister/OneDrive - Research Center Pharmaceutical Engineering/ACS GCIPR_01/100_Data/Amidation/Julia/Julia_Kinetics/experiments/PS-GCIPR-27-20230220_V01.xlsx"

function importexperiment(filename)

    xf = XLSX.readxlsx(filename)
    HITEC = xf["HITEC"]
    FTIR = xf["FTIR"]

    ## determine max. measurement durations
    maxrow_hitec = HITEC.dimension.stop.row_number
    maxrow_ftir = FTIR.dimension.stop.row_number

    t_maxidx_hitec = findfirst(isequal(missing), HITEC["B5:B"*string(maxrow_hitec)])

    if t_maxidx_hitec != nothing
        maxidx_hitec = t_maxidx_hitec.I[1]
    else
        maxidx_hitec = maxrow_hitec
    end

    t_maxidx_ftir = findfirst(isequal(missing), FTIR["B2:B"*string(maxrow_ftir)])

    if t_maxidx_ftir != nothing
        maxidx_ftir = t_maxidx_ftir.I[1]
    else
        maxidx_ftir = maxrow_ftir
    end

    #read data
    time_hitec = (HITEC["C5:C"*string(maxidx_hitec)])       # time in seconds for hitec values
    time_ftir =  (FTIR["D2:D"*string(maxidx_ftir)])         # time in seconds for ftir values
    flowrate = [
        float.(HITEC["D5:D"*string(maxidx_hitec)] * 1E-3 / (60))'; # [l/s]    flowrate ester
        float.(HITEC["E5:E"*string(maxidx_hitec)] * 1E-3 / (60))'; # [l/s]    flowrate amine 
        float.(HITEC["F5:F"*string(maxidx_hitec)] * 1E-3 / (60))'; # [l/s]      flowrate tbd
        float.(HITEC["G5:G"*string(maxidx_hitec)] * 1E-3 / (60))'; # [l/s]  flowrate product
        float.(HITEC["H5:H"*string(maxidx_hitec)] * 1E-3 / (60))'; # [l/s]  flowrate solvent
    ]

    flowrate_total = sum(flowrate, dims = 1)
    
    Tset = float.(HITEC["N5:N"*string(maxidx_hitec)]) #temperature actual of hotcoil
 # calculate inlet concentration of each feed
    C_ester_in_    = flowrate[1,:]
    C_amine_in_    = flowrate[2,:]
    C_TBD_in_      = flowrate[3,:]
    C_product_in_  = flowrate[4,:]
    C_intermediate_in_  = flowrate[4,:]
    C_alcohol_in_  = flowrate[4,:]

    # MT: you could use broadcasting here to write nicer code and maybe improve performace
    for i in range(1,(length(flowrate[1,:])))
        C_ester_in_[i] = (flowrate[1,i] * 2.5 / (flowrate_total[i])) # [l/s]    flowrate ester  
        C_amine_in_[i] = (flowrate[2,i] * 2.5 / (flowrate_total[i]))
    #    C_TBD_in_[i] = ((flowrate[3,i] * 0.5 / (flowrate_total[i])) +  (flowrate[1,i] * 2.2 / (flowrate_total[i])))
        C_TBD_in_[i] = ((flowrate[3,i] * 0.5 / (flowrate_total[i])))
        C_product_in_[i] = (flowrate[4,i] * 1.0 / (flowrate_total[i]))
        C_intermediate_in_[i] = (1e-12)
        C_alcohol_in_[i] = (1e-12)
    end

    C_inlet = [
        (C_ester_in_)';
        (C_amine_in_)'; 
        (C_TBD_in_)';
        (C_product_in_)';
        (C_intermediate_in_)';
        (C_alcohol_in_)';
    ]
# concentration at FTIR
    c_ftir = [
        float.(FTIR["E2:E"*string(maxidx_ftir)])'; # Predicted [Ester]
        float.(FTIR["F2:F"*string(maxidx_ftir)])'; # Predicted [Amine]
        float.(FTIR["G2:G"*string(maxidx_ftir)])'; # Predicted [TBD]
        float.(FTIR["H2:H"*string(maxidx_ftir)])'; # Predicted [Product]
        0 .*float.(FTIR["H2:H"*string(maxidx_ftir)])'; # Predicted [intermediate] ==> zero, because it is not measured
        0 .*float.(FTIR["H2:H"*string(maxidx_ftir)])'; # Predicted [intermediate] ==> zero, because it is not measured
    ]

    c_ftir_sum = [
        float.(FTIR["E2:E"*string(maxidx_ftir)])'; # Predicted [Ester]
        float.(FTIR["H2:H"*string(maxidx_ftir)])'; # Predicted [Product]
    ]

    sum_ftir = sum(c_ftir_sum, dims=1) #sum of predicted [ester] and [product] -> mass balance?

    #c_ftir = [c_ftir;sum_ftir]

print("import complete")
    return ExperimentalData(time_hitec, flowrate, flowrate_total, Tset, C_inlet, time_ftir, c_ftir, max(time_ftir[end],time_hitec[end]))

end

#function toseconds(a::XLSX.Dates.Time)
#    return a.instant.value * 1E-9
#end



function saveparams(params::T, filename::String) where T <: Vector
    serialize(filename,params)
end

