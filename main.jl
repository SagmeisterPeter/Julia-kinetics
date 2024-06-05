using OptSynthesis
using Metaheuristics
using Plots#; plotly()
using DataFrames
using XLSX


run = 1

filenames = [
"./experiments/PS-GCIPR-27-20230220_V01.xlsx"; 
"./experiments/LM-GCIPR-9.xlsx"; 
"./experiments/LM-GCIPR-10_V02.xlsx"; 
"./experiments/LM-GCIPR-11_V02.xlsx"; 
"./experiments/LM-GCIPR-12.xlsx"; 
"./experiments/LM-GCIPR-31.xlsx"; 
"./experiments/LM-GCIPR-32.xlsx"; 
"./experiments/LM-GCIPR-34.xlsx"; 
"./experiments/LM-GCIPR-39.xlsx"; 
"./experiments/LM-GCIPR-40.xlsx"; 
"./experiments/LM-GCIPR-41.xlsx"; 
"./experiments/LM-GCIPR-50.xlsx";
"./experiments/LM-GCIPR-55.xlsx";
]

excel_file_path = "PS-27-test-results.xlsx"

c0 = 
[
0.4;
0.396;
]
 
# returns a SynthesisModel structd
#@time model = OptSynthesis.getreactionmodel("tworeactions_fixed",c0[1],filenames[11]);
@time model = OptSynthesis.getreactionmodel("tworeactions_fixed",c0[1],filenames[2]);
#@time model = OptSynthesis.getreactionmodel("tworeaction_fixed",c0[run],filenames[1]);
# executing the identifyparameters function from ReactionsModels.jl (different param order)
#p_est = OptSynthesis.identifyparameters(model,w = [1000,1000,1000,1000,1000,1000], refine = true)
p_est = OptSynthesis.identifyparameters(model,w = [1,1,0,1,1], refine = true)
#p_est_new = OptSynthesis.identifyparametersnew(model,w = [1,1,0,1,0])

#p_est_ECA = OptSynthesis.identifyparametersnew("tworeactions_fixed", "ECA", model; w = [1,1,0,1,0])
#p_est_DE = OptSynthesis.identifyparametersnew("tworeactions_fixed", "DE", model; w = [1,1,0,1,0])
#p_est_PSO = OptSynthesis.identifyparametersnew("tworeactions_fixed", "PSO", model; w = [1,1,0,1,0])
#p_est_ABC = OptSynthesis.identifyparametersnew("tworeactions_fixed", "ABC", model; w = [1,1,0,1,0])
#p_est_CGSA = OptSynthesis.identifyparametersnew("tworeactions_fixed", "CGSA", model; w = [1,1,0,1,0])
#p_est_WOA = OptSynthesis.identifyparametersnew("tworeactions_fixed", "WOA", model; w = [1,1,0,1,0])
#p_est_MCCGA = OptSynthesis.identifyparametersnew("tworeactions_fixed", "MCCGA", model; w = [1,1,0,1,0])

# gave errors
#p_est_GA = OptSynthesis.identifyparametersnew("tworeactions_fixed", "GA", model; w = [1,1,0,1,0])
#p_est_CCMO = OptSynthesis.identifyparametersnew("tworeactions_fixed", "CCMO", model; w = [1,1,0,1,0])
#p_est_εDE = OptSynthesis.identifyparametersnew("tworeactions_fixed", "εDE", model; w = [1,1,0,1,0])
#p_est_BRKGA = OptSynthesis.identifyparametersnew("tworeactions_fixed", "BRKGA", model; w = [1,1,0,1,0])
#p_est_NSGA2 = OptSynthesis.identifyparametersnew("tworeactions_fixed", "NSGA2", model; w = [1,1,0,1,0])
#p_est_NSGA3 = OptSynthesis.identifyparametersnew("tworeactions_fixed", "NSGA3", model; w = [1,1,0,1,0])
#p_est_SPEA2 = OptSynthesis.identifyparametersnew("tworeactions_fixed", "SPEA2", model; w = [1,1,0,1,0])

#p_est_ABC1 = OptSynthesis.identifyparametersnew("tworeactions_fixed", "ABC", model; w = [1,1,0,1,0])
#p_est_ABC2 = OptSynthesis.identifyparametersnew("tworeactions_fixed", "ABC", model; w = [1,1,0,1,0])
#p_est_ABC3 = OptSynthesis.identifyparametersnew("tworeactions_fixed", "ABC", model; w = [1,1,0,1,0])

print(p_est)
#p_est_new = p_est_MCCGA

#p_est = Metaheuristics.minimizer(p_est_new)
# solve_model is again executing the internal julia function solve! with several params
@time (t,sol_test) = OptSynthesis.solve_model(model,p_est);
# plotting the results via the Plots.jl extension (using at the top)
OptSynthesis.plotidentificationresults(model,t,sol_test)
sol_test
#@time modelSS = OptSynthesis.getreactionmodelSS("tworeactions_fixed", c0[1],p_est);
@time modelSS = OptSynthesis.getreactionmodelSS("tworeaction_fixed", c0[1],p_est);
################################################
df_1 = DataFrames.DataFrame(time_ftir = vec(model.exdata.time_ftir), 
                          Ester_ftir = model.exdata.c_ftir[1,:], 
                          Amine_ftir = model.exdata.c_ftir[2,:],
                          TBD_ftir = model.exdata.c_ftir[3,:],
                          Product_ftir = model.exdata.c_ftir[4,:])

df_2 = DataFrames.DataFrame(time_sim = t,
                          Ester_sim = sol_test[1], 
                          Amine_sim = sol_test[2],
                          TBD_sim = sol_test[3],
                          Product_sim = sol_test[4])
                      #    inter_sim = sol_test[5])                       
XLSX.writetable(("FTIR" * excel_file_path), df_1; sheetname="FTIR")
XLSX.writetable(("sim" * excel_file_path), df_2; sheetname="Sim")

#####################################################
M = 199.17 #--> molar mass product in g/mol
V = 2.79  #--> volume reactor in mL 

# optimal operating points 

q0_init = 5.13E-5;
T0_init = 80+273.15;
C_ester_in_0_init = 0.3
C_amine_in_0_init = 0.3
C_TBD_in_0_init = 0.1

#q0_min = 9.45E-6  
#q0_max = 9.45E-5

q0_max = ((V / 1) * 1E-3) / 60  # Volume reactor / residence time --> convert to [l/s]
q0_min = ((V / 10) * 1E-3) / 60  # Volume reactor / residence time --> convert to [l/s]            
T0_min = 115+273.15;
T0_max = 120+273.15;

#C_ester_in_0_min = 0.499
C_ester_in_0_min = 0.3
C_amine_in_0_min = 0.99
C_TBD_in_0_min = 0.405
#C_TBD_in_0_min = 0.2499
C_ester_in_0_max = 0.5
C_amine_in_0_max = 1.0
C_TBD_in_0_max = 1.0

Eq_amine_min = 1.99
Eq_amine_max = 2.0
#Eq_TBD_min = 0.499
Eq_TBD_min = 1.35
Eq_TBD_max = 2.0

res = OptSynthesis.steadystate(modelSS,[T0_init, q0_max, 0.35, C_amine_in_0_init, C_TBD_in_0_init])
res = OptSynthesis.steadystate(modelSS,[T0_max, q0_min, C_ester_in_0_max, C_amine_in_0_max, C_TBD_in_0_max])
#res = OptSynthesis.steadystate(modelSS,[T0_min,q0_min, C_ester_in_0_min, C_amine_in_0_min, C_TBD_in_0_min])
#prod = OptSynthesis.productivity(res[4], q0_max, M) 
#sel = OptSynthesis.selectivity(res[4],res[1], C_ester_in_0_max)
#con = OptSynthesis.conversion(res[1], C_ester_in_0_max)

# cost function for finding pareto front
function cost(modelSS, p)
    M = 199.17 #--> molar mass product in g/mol
    V = 2.79   #--> volume reactor in mL 

    ss = OptSynthesis.steadystate(modelSS, p)

    prod = OptSynthesis.productivity(ss[4],p[2],M) 
    sel = OptSynthesis.selectivity(ss[4],ss[1],p[3])
    con = OptSynthesis.conversion(ss[1],p[3])
  #  space = Opt.Synthesis.spacetimeyield(ss[4],p[2],M,V)
    yield = ss[4]/p[3]
    if con < 0
       con = 0.000000000001
    end
    if prod < 0
       prod = 0.000000000001
    end


  #  sty = OptSynthesis.spacetimeyield(ss[4],p[2];M,V)
    fx = [(-prod); (-con)]

    gx = [0.0]
    hx = [0.0]

    # order is important
    return fx, gx, hx
end

f(p) = cost(modelSS,p)
res = f([T0_max, q0_max, C_ester_in_0_init, C_amine_in_0_max, C_TBD_in_0_max])

bounds = [T0_min T0_max;
          q0_min q0_max;
          C_ester_in_0_min C_ester_in_0_max;
          C_amine_in_0_min C_amine_in_0_max;
          C_TBD_in_0_min C_TBD_in_0_max;
          ]';

# Common options
options = Options(seed=1, store_convergence = true)

# Optimize
@time result = Metaheuristics.optimize(f, bounds, NSGA2(options=options))

# plot results
X = positions(result)
display(
scatter(X[:,3],(X[:,5] ./ X[:,3]),
xlabel="concentration ester [mol/]",
ylabel="equivalents base [ml/min]",
title="optimal parameters",
label="Pareto opt. parameters",
xlims=[0;1],
ylims=[0;1]
))

# choose random operating points to show pareto optimality
nrand = 10000
Trand = (T0_max-T0_min).*rand(Float64,nrand) .+ T0_min
#Trand = 0 .*rand(Float64,nrand) .+ T0_max
qrand = (q0_max-q0_min).*rand(Float64,nrand) .+ q0_min
C_ester_in_rand = (C_ester_in_0_max - C_ester_in_0_min).*rand(Float64,nrand) .+ C_ester_in_0_min

Eq_amine_rand = (Eq_amine_max - Eq_amine_min).*rand(Float64,nrand) .+ Eq_amine_min
Eq_TBD_rand =   (Eq_TBD_max - Eq_TBD_min).*rand(Float64,nrand) .+ Eq_TBD_min
#Eq_TBD_rand =   0 .*rand(Float64,nrand) .+ Eq_TBD_max

C_amine_in_rand_    = zeros(nrand,1)
C_TBD_in_rand_    = zeros(nrand,1)

for i in range(1,(length(C_ester_in_rand))) 
  C_amine_in_rand_[i] = C_ester_in_rand[i] * Eq_amine_rand[i] 
end

for i in range(1,(length(C_ester_in_rand))) 
  C_TBD_in_rand_[i] = C_ester_in_rand[i] * Eq_TBD_rand[i] 
end


sel_rand = zeros(nrand,1)
prod_rand = zeros(nrand,1)
Objective1_rand = zeros(nrand,1)
Objective2_rand = zeros(nrand,1)


for i = 1:nrand
    fx,gx,hx = f([Trand[i];qrand[i];C_ester_in_rand[i];C_amine_in_rand_[i];C_TBD_in_rand_[i];])

    Objective1_rand[i] = -fx[1]
    Objective2_rand[i] = -fx[2]
end

A = pareto_front(result)
fpareto = scatter(-A[:, 1], -A[:,2]*100, 
label="pareto optimal", 
xlabel = "productivity [g/h]", 
ylabel = "Conversion [%]",
title = "optimal operating points",
#xlims=[0;1],
#ylims=[0;100])
)

# positions in matrix NxD 
scatter!(fpareto,Objective1_rand,Objective2_rand*100,label="randomly sampled")

#vec(sel_rand)
# export data
#using DataFrames
df_random = DataFrames.DataFrame(Temperature = Trand, 
                          Total_Flow = qrand, 
                          Ester = vec(C_ester_in_rand),
                          Amine = vec(C_amine_in_rand_),
                          TBD = vec(C_TBD_in_rand_),
                          tres = (V ./ (qrand .*60 ./1E-3)),
                          productivity = vec(Objective1_rand),
                          conversion = vec(Objective2_rand))
                          

df_pareto= DataFrames.DataFrame(Temperature = X[:,1], 
                                Total_Flow = X[:,2], 
                                Ester = X[:,3],
                                Amine = X[:,4],
                                TBD = X[:,5],
                                tres =(V ./ (X[:,2] .*60 ./1E-3)),
                                productivity = -A[:,1],
                                conversion = -A[:,2])


XLSX.writetable("random_sampled_3_" * excel_file_path, df_random)
XLSX.writetable("pareto_" * excel_file_path, df_pareto)

#Plot1 = OptSynthesis.plot_paracoords(df_random,"paracoord_random")
#Plot2 = OptSynthesis.plot_paracoords(df_,"paracoord_pareto")

# X = positions(result)
#display(scatter!(prod_rand*60*1E3,sel_rand*100,
#label="pareto optimal", 
#xlabel = "productivity [mmol/min]", 
#ylabel = "Conversion [%]",
#title = "optimal operating points for 4-EP (Pd catalyst)"))

