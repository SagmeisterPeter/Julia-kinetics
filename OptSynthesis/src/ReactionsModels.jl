struct SynthesisModel
    odeprob
    materials
    auxmaterial
    exdata
end

struct SynthesisModelSS
    odeprob
    materials
end


"""
    arrhenius(C,ki,Ea,kr,T)

TBW
"""


#function arrhenius(C1, kr1, ki, Ea, T)
#   R = 8.3144598
#  #     kr = 1;
#  return max(0, C1)^kr1 * ki * 10^4 * exp(-(Ea * 10^3) / (R * T))           #fitting negative number thats why concentration is max 0
#end

function arrhenius2(C1, C2, kr1, kr2, ki, Ea, T)
    R = 8.3144598
    #return max(0, C1)^kr1 * max(0, C2)^kr2 * ki * 10^4 * exp(-(Ea * 10^4) / (R * T))         #fitting negative number thats why concentration is max 0    
    return max(0, C1)^kr1 * max(0, C2)^kr2 * ki * exp(-(Ea * 10^3) / (R * T))           #took out the 10^4 for the ki
 end
 
function arrhenius3(C1, C2, C3, kr1, kr2, kr3, ki, Ea, T)
    R = 8.3144598
        
  #  return max(0, C1)^kr1 * max(0, C2)^kr2 * max(0, C3)^kr3* ki * 10^4 * exp(-(Ea * 10^4) / (R * T))     
    return max(0, C1)^kr1 * max(0, C2)^kr2 * max(0, C3)^kr3 * ki * exp(-(Ea * 10^3) / (R * T))    # took out the 10^4       #fitting negative number thats why concentration is max 0
end 


"""
    getReactionModel(catalyst)

returns the reactor model depending on the specific catalyst.
"""
#function getreactionmodel(catalyst, Cin, T, q, tsim=1000)
function getreactionmodel(catalyst, C_ester_in, C_amine_in, C_TBD_in, C_product_in, C_intermediate_in, C_alcohol_in, T, q, tsim)
    # reactor parameters
    V = 5.67E-3 # volume [dm^3]     # reactor parameters
 #   V = 2.79E-3 # volume [dm^3]     # reactor parameters
    numb_dz1 = 10 # number of discretization step
    L = V / 1 # L in [dm] ( --> A = 1 dm^2)
    Li = 0.2E-3 # inert length, second length after reactor   # you need atleast 5 compartments!
    numb_dz2 = 5 # number of discretization step
    dz = 1E-4 # discretization step size

    dz_1 = V / numb_dz1 # calculating discretization steps
    dz_2 = Li / numb_dz2 # calculating discretization steps

    #add different discretization steps to each reactor compartment ---> use number
        # discretization must be divided volume by the discretization step
    #simtime = 1000 # dummy simulation time
    auxmaterial = false

    @parameters t z1 z2

    Dt = Differential(t)
    Dz1 = Differential(z1)
    Dz2 = Differential(z2)

    if catalyst == "tworeactions"

        materials = ["ester"; "amine"; "TBD"; "product"; "intermediate"];
       # materials = ["ester"; "amine"; "TBD"; "product"];
       # materials = ["ester"; "amine"; "product"];
        

        @parameters kr[1:4] ki[1:2] Ea[1:2] # ki = reaction consts, Ea: act energy, kr reaction orders
        @variables  C[1:5](..) Ci[1:5](..) # concentration & conc at end

       # ┌ Warning:# The variable syntax (C[1:5])(..) is deprecated. Use (C(..))[1:5] instead.
       # │                   The former creates an array of functions, while the latter creates an array valued function.
      #  │                   The deprecated syntax will cause an error in the next major release of Symbolics.
       # │                   This change will facilitate better implementation of various features of Symbolics.

        rall = [
            arrhenius2(C[1](z1, t), C[3](z1, t), kr[1], kr[3], ki[1], Ea[1], T(t)), # r1
            arrhenius2(C[2](z1, t), C[5](z1, t), kr[2], kr[4], ki[2], Ea[2], T(t)), # r1

        ]


        # r1 = Ester + Base -> Intermediate  
        # r2 = Intermediate + Anmine -> Base + Product  

        # B are stoichometric factors -->  maybe change it to a different reaction network?? 
             #R1 R2 
        B = [-1   0          # C1  Ester
              0  -1          # C2  amine    
             -1   1          # C3  TBD
              0   1          # C4  Product formation
              1  -1          # intermdeiate
        ]


        eq = [Dt(C[1](z1, t)) ~ -q(t) *  Dz1(C[1](z1, t)) + B[1, :]' * rall, # dC1/dt
              Dt(C[2](z1, t)) ~ -q(t) *  Dz1(C[2](z1, t)) + B[2, :]' * rall, # dC2/dt 
              Dt(C[3](z1, t)) ~ -q(t) *  Dz1(C[3](z1, t)) + B[3, :]' * rall, # dC3/dt 
              Dt(C[4](z1, t)) ~ -q(t) *  Dz1(C[4](z1, t)) + B[4, :]' * rall, # dC4/dt 
              Dt(C[5](z1, t)) ~ -q(t) *  Dz1(C[5](z1, t)) + B[5, :]' * rall, # dC5/dt 
              Dt(Ci[1](z2, t)) ~ -q(t) * Dz2(Ci[1](z2, t)),                  # dC1/dt
              Dt(Ci[2](z2, t)) ~ -q(t) * Dz2(Ci[2](z2, t)),                  # dC2/dt 
              Dt(Ci[3](z2, t)) ~ -q(t) * Dz2(Ci[3](z2, t)),                  # dC3/dt 
              Dt(Ci[4](z2, t)) ~ -q(t) * Dz2(Ci[4](z2, t)),                  # dC4/dt 
              Dt(Ci[5](z2, t)) ~ -q(t) * Dz2(Ci[5](z2, t)),                  # dC4/dt 
        ]


        domains = [t ∈ Interval(0.0, tsim),
            z1 ∈ Interval(0.0, L),
            z2 ∈ Interval(L, L + Li)
        ]
       
       
        # boundary conditions
        bcs = [
            C[1](0, t) ~ C_ester_in(t),
            C[2](0, t) ~ C_amine_in(t), 
            C[3](0, t) ~ C_TBD_in(t), 
            C[4](0, t) ~ C_product_in(t),
            C[5](0, t) ~ C_intermediate_in(t),    #intermediate  set to small nonzero value (dirty bugfix), see https://github.com/SciML/MethodOfLines.jl/issues/194       
            C[1](L, t) ~ Ci[1](L, t),   # interface BCs
            C[2](L, t) ~ Ci[2](L, t), # 
            C[3](L, t) ~ Ci[3](L, t), # 
            C[4](L, t) ~ Ci[4](L, t), # 
            C[5](L, t) ~ Ci[5](L, t), # intermediate
        ]

        ics = [
            C[1](z1, 0) ~ 0,         #initial conditions
            C[2](z1, 0) ~ 0,
            C[3](z1, 0) ~ 0,
            C[4](z1, 0) ~ 0,
            C[5](z1, 0) ~ 0,
            Ci[1](z2, 0) ~ 0,         #initial conditions
            Ci[2](z2, 0) ~ 0,
            Ci[3](z2, 0) ~ 0,
            Ci[4](z2, 0) ~ 0,
            Ci[5](z2, 0) ~ 0,
        ]

        bcs_ics = [bcs; ics]

        @named pdesys = PDESystem(eq, bcs_ics, domains, [z1, z2, t], [
                C[1](z1, t), 
                C[2](z1, t), 
                C[3](z1, t), 
                C[4](z1, t), 
                C[5](z1, t),
                Ci[1](z2, t), 
                Ci[2](z2, t), 
                Ci[3](z2, t), 
                Ci[4](z2, t), 
                Ci[5](z2, t),
            ], 
            [
                kr[1] => 1,
                kr[2] => 1,
                kr[3] => 1,
                kr[4] => 1,
                ki[1] => 7,
                ki[2] => 7,
                Ea[1] => 50,
                Ea[2] => 50,

            ])



    

    elseif catalyst == "tworeactions_fixed"

        materials = ["ester"; "amine"; "TBD"; "product"; "intermediate"];
       # materials = ["ester"; "amine"; "TBD"; "product"];
       # materials = ["ester"; "amine"; "product"];
        

        @parameters ki[1:2] Ea[1:2] # ki = reaction consts, Ea: act energy, kr reaction orders
        @variables  C[1:5](..) Ci[1:5](..) # concentration & conc at end


        rall = [
            arrhenius2(C[1](z1, t), C[3](z1, t), 1, 1, ki[1], Ea[1], T(t)), # r1
            arrhenius2(C[2](z1, t), C[5](z1, t), 1, 1, ki[2], Ea[2], T(t)), # r1

        ]


        # r1 = Ester + Base -> Intermediate  
        # r2 = Intermediate + Anmine -> Base + Product  

        # B are stoichometric factors -->  maybe change it to a different reaction network?? 
             #R1 R2 
        B = [-1   0          # C1  Ester
              0  -1          # C2  amine    
             -1   1          # C3  TBD
              0   1          # C4  Product formation
              1  -1          # intermdeiate
        ]


        eq = [Dt(C[1](z1, t)) ~ -q(t) *  Dz1(C[1](z1, t)) + B[1, :]' * rall, # dC1/dt
              Dt(C[2](z1, t)) ~ -q(t) *  Dz1(C[2](z1, t)) + B[2, :]' * rall, # dC2/dt 
              Dt(C[3](z1, t)) ~ -q(t) *  Dz1(C[3](z1, t)) + B[3, :]' * rall, # dC3/dt 
              Dt(C[4](z1, t)) ~ -q(t) *  Dz1(C[4](z1, t)) + B[4, :]' * rall, # dC4/dt 
              Dt(C[5](z1, t)) ~ -q(t) *  Dz1(C[5](z1, t)) + B[5, :]' * rall, # dC5/dt 
              Dt(Ci[1](z2, t)) ~ -q(t) * Dz2(Ci[1](z2, t)), # dC1/dt
              Dt(Ci[2](z2, t)) ~ -q(t) * Dz2(Ci[2](z2, t)), # dC2/dt 
              Dt(Ci[3](z2, t)) ~ -q(t) * Dz2(Ci[3](z2, t)), # dC3/dt 
              Dt(Ci[4](z2, t)) ~ -q(t) * Dz2(Ci[4](z2, t)), # dC4/dt 
              Dt(Ci[5](z2, t)) ~ -q(t) * Dz2(Ci[5](z2, t)), # dC4/dt 
        ]


        domains = [t ∈ Interval(0.0, tsim),
            z1 ∈ Interval(0.0, L),
            z2 ∈ Interval(L, L + Li)
        ]
       
       
        # boundary conditions
        bcs = [
            C[1](0, t) ~ C_ester_in(t),
            C[2](0, t) ~ C_amine_in(t), 
            C[3](0, t) ~ C_TBD_in(t), 
            C[4](0, t) ~ C_product_in(t),
            C[5](0, t) ~ C_intermediate_in(t),    #intermediate  set to small nonzero value (dirty bugfix), see https://github.com/SciML/MethodOfLines.jl/issues/194
            C[1](L, t) ~ Ci[1](L, t),   # interface BCs
            C[2](L, t) ~ Ci[2](L, t), # 
            C[3](L, t) ~ Ci[3](L, t), # 
            C[4](L, t) ~ Ci[4](L, t), # 
            C[5](L, t) ~ Ci[5](L, t), # intermediate
        ]

        ics = [
            C[1](z1, 0) ~ 0,         #initial conditions
            C[2](z1, 0) ~ 0,
            C[3](z1, 0) ~ 0,
            C[4](z1, 0) ~ 0,
            C[5](z1, 0) ~ 0,
            Ci[1](z2, 0) ~ 0,         #initial conditions
            Ci[2](z2, 0) ~ 0,
            Ci[3](z2, 0) ~ 0,
            Ci[4](z2, 0) ~ 0,
            Ci[5](z2, 0) ~ 0,
        ]

        bcs_ics = [bcs; ics]

        @named pdesys = PDESystem(eq, bcs_ics, domains, [z1, z2, t], [
                C[1](z1, t), 
                C[2](z1, t), 
                C[3](z1, t), 
                C[4](z1, t), 
                C[5](z1, t),
                Ci[1](z2, t), 
                Ci[2](z2, t), 
                Ci[3](z2, t), 
                Ci[4](z2, t), 
                Ci[5](z2, t),
            ], 
            [
                ki[1] => 7,
                ki[2] => 7,
                Ea[1] => 50,
                Ea[2] => 50,

            ])

        elseif catalyst == "threereactions_fixed"

            materials = ["ester"; "amine"; "TBD"; "product"; "intermediate";"alcohol"];
           # materials = ["ester"; "amine"; "TBD"; "product"];
           # materials = ["ester"; "amine"; "product"];
            
    
            @parameters ki[1:3] Ea[1:3] # ki = reaction consts, Ea: act energy, kr reaction orders
            @variables  C[1:6](..) Ci[1:6](..) # concentration & conc at end
    
    
            rall = [
                arrhenius2(C[1](z1, t), C[3](z1, t), 1, 1, ki[1], Ea[1], T(t)), # r1
                arrhenius2(C[5](z1, t), C[6](z1, t), 1, 1, ki[2], Ea[2], T(t)), # r2
                arrhenius2(C[2](z1, t), C[5](z1, t), 1, 1, ki[3], Ea[3], T(t)), # r3
    
            ]
    
    
            # r1 = Ester + Base -> Intermediate  
            # r2 = Intermediate + Anmine -> Base + Product  
    
            # B are stoichometric factors -->  maybe change it to a different reaction network?? 
                 #R1 R2 
            B = [-1  1  0          # C1  Ester
                  0  0 -1          # C2  amine    
                 -1  1  1          # C3  TBD
                  0  0  1          # C4  Product formation
                  1 -1 -1          # intermdeiate
                  1 -1  0          # alcohol
            ]
    
    
            eq = [Dt(C[1](z1, t)) ~ -q(t) *  Dz1(C[1](z1, t)) + B[1, :]' * rall, # dC1/dt
                  Dt(C[2](z1, t)) ~ -q(t) *  Dz1(C[2](z1, t)) + B[2, :]' * rall, # dC2/dt 
                  Dt(C[3](z1, t)) ~ -q(t) *  Dz1(C[3](z1, t)) + B[3, :]' * rall, # dC3/dt 
                  Dt(C[4](z1, t)) ~ -q(t) *  Dz1(C[4](z1, t)) + B[4, :]' * rall, # dC4/dt 
                  Dt(C[5](z1, t)) ~ -q(t) *  Dz1(C[5](z1, t)) + B[5, :]' * rall, # dC5/dt 
                  Dt(C[6](z1, t)) ~ -q(t) *  Dz1(C[6](z1, t)) + B[6, :]' * rall, # dC5/dt 
                  Dt(Ci[1](z2, t)) ~ -q(t) * Dz2(Ci[1](z2, t)), # dC1/dt
                  Dt(Ci[2](z2, t)) ~ -q(t) * Dz2(Ci[2](z2, t)), # dC2/dt 
                  Dt(Ci[3](z2, t)) ~ -q(t) * Dz2(Ci[3](z2, t)), # dC3/dt 
                  Dt(Ci[4](z2, t)) ~ -q(t) * Dz2(Ci[4](z2, t)), # dC4/dt 
                  Dt(Ci[5](z2, t)) ~ -q(t) * Dz2(Ci[5](z2, t)), # dC4/dt 
                  Dt(Ci[6](z2, t)) ~ -q(t) * Dz2(Ci[6](z2, t)), # dC4/dt 
            ]
    
    
            domains = [t ∈ Interval(0.0, tsim),
                z1 ∈ Interval(0.0, L),
                z2 ∈ Interval(L, L + Li)
            ]
           
           
            # boundary conditions
            bcs = [
                C[1](0, t) ~ C_ester_in(t),
                C[2](0, t) ~ C_amine_in(t), 
                C[3](0, t) ~ C_TBD_in(t), 
                C[4](0, t) ~ C_product_in(t),
                C[5](0, t) ~ C_intermediate_in(t),    #intermediate  set to small nonzero value (dirty bugfix), see https://github.com/SciML/MethodOfLines.jl/issues/194
                C[6](0, t) ~ C_alcohol_in(t),
                C[1](L, t) ~ Ci[1](L, t),   # interface BCs
                C[2](L, t) ~ Ci[2](L, t), # 
                C[3](L, t) ~ Ci[3](L, t), # 
                C[4](L, t) ~ Ci[4](L, t), # 
                C[5](L, t) ~ Ci[5](L, t), # intermediate
                C[6](L, t) ~ Ci[6](L, t),
            ]
    
            ics = [
                C[1](z1, 0) ~ 0,         #initial conditions
                C[2](z1, 0) ~ 0,
                C[3](z1, 0) ~ 0,
                C[4](z1, 0) ~ 0,
                C[5](z1, 0) ~ 0,
                C[6](z1, 0) ~ 0,
                Ci[1](z2, 0) ~ 0,         #initial conditions
                Ci[2](z2, 0) ~ 0,
                Ci[3](z2, 0) ~ 0,
                Ci[4](z2, 0) ~ 0,
                Ci[5](z2, 0) ~ 0,
                Ci[6](z2, 0) ~ 0,
            ]
    
            bcs_ics = [bcs; ics]
    
            @named pdesys = PDESystem(eq, bcs_ics, domains, [z1, z2, t], [
                    C[1](z1, t), 
                    C[2](z1, t), 
                    C[3](z1, t), 
                    C[4](z1, t), 
                    C[5](z1, t),
                    C[6](z1, t),
                    Ci[1](z2, t), 
                    Ci[2](z2, t), 
                    Ci[3](z2, t), 
                    Ci[4](z2, t), 
                    Ci[5](z2, t),
                    Ci[6](z2, t),
                ], 
                [
                    ki[1] => 10000,
                    ki[2] => 10000,
                    ki[3] => 10000,
                    Ea[1] => 50,
                    Ea[2] => 50,
                    Ea[3] => 50,
                ])
    



    

    elseif catalyst == "onereaction"

        # materials = ["ester"; "amine"; "TBD"; "product"; "intermediate"];
         materials = ["ester"; "amine"; "TBD"; "product"];
 
         @parameters kr[1:3] ki[1:1] Ea[1:1] # ki = reaction consts, Ea: act energy, kr reaction orders
         @variables  C[1:4](..) Ci[1:4](..) # concentration & conc at end
 
 
         rall = [
             arrhenius3(C[1](z1, t), C[2](z1, t), C[3](z1, t), kr[1], kr[2], kr[3], ki[1], Ea[1], T(t)), # r1
         ]
 
 
         # r1 = Ester + Base + Anmine -> Base + Product  
 
         # B are stoichometric factors -->  maybe change it to a different reaction network?? 
              #R1 R2 
         B = [-1  # 0          # C1  Ester
              -1  #1          # C2  amine    
               0   #1          # C3  TBD
               1   #1          # C4  Product formation
            #   1  -1          # intermdeiate
         ]
 
 
         eq = [Dt(C[1](z1, t)) ~ -q(t) *  Dz1(C[1](z1, t)) + B[1, :]' * rall, # dC1/dt
               Dt(C[2](z1, t)) ~ -q(t) *  Dz1(C[2](z1, t)) + B[2, :]' * rall, # dC2/dt 
               Dt(C[3](z1, t)) ~ -q(t) *  Dz1(C[3](z1, t)) + B[3, :]' * rall, # dC3/dt 
               Dt(C[4](z1, t)) ~ -q(t) *  Dz1(C[4](z1, t)) + B[4, :]' * rall, # dC4/dt 
               Dt(Ci[1](z2, t)) ~ -q(t) * Dz2(Ci[1](z2, t)), # dC1/dt
               Dt(Ci[2](z2, t)) ~ -q(t) * Dz2(Ci[2](z2, t)), # dC2/dt 
               Dt(Ci[3](z2, t)) ~ -q(t) * Dz2(Ci[3](z2, t)), # dC3/dt 
               Dt(Ci[4](z2, t)) ~ -q(t) * Dz2(Ci[4](z2, t)), # dC4/dt  
         ]
 
 
         domains = [t ∈ Interval(0.0, tsim),
             z1 ∈ Interval(0.0, L),
             z2 ∈ Interval(L, L + Li)
         ]
                
         # boundary conditions
         bcs = [
             C[1](0, t) ~ C_ester_in(t),
             C[2](0, t) ~ C_amine_in(t), 
             C[3](0, t) ~ C_TBD_in(t), 
             C[4](0, t) ~ C_product_in(t),
             C[1](L, t) ~ Ci[1](L, t),       # interface BCs
             C[2](L, t) ~ Ci[2](L, t), # 
             C[3](L, t) ~ Ci[3](L, t), # 
             C[4](L, t) ~ Ci[4](L, t), # 
         ]
 
         ics = [
             C[1](z1, 0) ~ 0,         #initial conditions
             C[2](z1, 0) ~ 0,
             C[3](z1, 0) ~ 0,
             C[4](z1, 0) ~ 0,
             Ci[1](z2, 0) ~ 0,         #initial conditions
             Ci[2](z2, 0) ~ 0,
             Ci[3](z2, 0) ~ 0,
             Ci[4](z2, 0) ~ 0,
         ]
 
         bcs_ics = [bcs; ics]
 
         @named pdesys = PDESystem(eq, bcs_ics, domains, [z1, z2, t], [
                 C[1](z1, t), 
                 C[2](z1, t), 
                 C[3](z1, t), 
                 C[4](z1, t), 
                 Ci[1](z2, t), 
                 Ci[2](z2, t), 
                 Ci[3](z2, t), 
                 Ci[4](z2, t), 
             ], 
             [
                 kr[1] => 1,
                 kr[2] => 1,
                 kr[3] => 1,
                 ki[1] => 5,
                 Ea[1] => 10,
             ])
    

    elseif catalyst == "onereaction_fixed"

        # materials = ["ester"; "amine"; "TBD"; "product"; "intermediate"];
        materials = ["ester"; "amine"; "TBD"; "product"];

        @parameters ki[1:1] Ea[1:1] # ki = reaction consts, Ea: act energy, kr reaction orders
        @variables  C[1:4](..) Ci[1:4](..) # concentration & conc at end


        rall = [
            arrhenius3(C[1](z1, t), C[2](z1, t), C[3](z1, t), 1, 1, 1, ki[1], Ea[1], T(t)), # r1
        ]


        # r1 = Ester + Base + Anmine -> Base + Product  

        # B are stoichometric factors -->  maybe change it to a different reaction network?? 
            #R1 R2 
        B = [-1  # 0          # C1  Ester
            -1  #1          # C2  amine    
            0   #1          # C3  TBD
            1   #1          # C4  Product formation
            #   1  -1          # intermdeiate
        ]


        eq = [Dt(C[1](z1, t)) ~ -q(t) *  Dz1(C[1](z1, t)) + B[1, :]' * rall, # dC1/dt
            Dt(C[2](z1, t)) ~ -q(t) *  Dz1(C[2](z1, t)) + B[2, :]' * rall, # dC2/dt 
            Dt(C[3](z1, t)) ~ -q(t) *  Dz1(C[3](z1, t)) + B[3, :]' * rall, # dC3/dt 
            Dt(C[4](z1, t)) ~ -q(t) *  Dz1(C[4](z1, t)) + B[4, :]' * rall, # dC4/dt 
            Dt(Ci[1](z2, t)) ~ -q(t) * Dz2(Ci[1](z2, t)), # dC1/dt
            Dt(Ci[2](z2, t)) ~ -q(t) * Dz2(Ci[2](z2, t)), # dC2/dt 
            Dt(Ci[3](z2, t)) ~ -q(t) * Dz2(Ci[3](z2, t)), # dC3/dt 
            Dt(Ci[4](z2, t)) ~ -q(t) * Dz2(Ci[4](z2, t)), # dC4/dt  
        ]


        domains = [t ∈ Interval(0.0, tsim),
            z1 ∈ Interval(0.0, L),
            z2 ∈ Interval(L, L + Li)
        ]
                
        # boundary conditions
        bcs = [
            C[1](0, t) ~ C_ester_in(t),
            C[2](0, t) ~ C_amine_in(t), 
            C[3](0, t) ~ C_TBD_in(t), 
            C[4](0, t) ~ C_product_in(t),
            C[1](L, t) ~ Ci[1](L, t),       # interface BCs
            C[2](L, t) ~ Ci[2](L, t), # 
            C[3](L, t) ~ Ci[3](L, t), # 
            C[4](L, t) ~ Ci[4](L, t), # 
        ]

        ics = [
            C[1](z1, 0) ~ 0,         #initial conditions
            C[2](z1, 0) ~ 0,
            C[3](z1, 0) ~ 0,
            C[4](z1, 0) ~ 0,
            Ci[1](z2, 0) ~ 0,         #initial conditions
            Ci[2](z2, 0) ~ 0,
            Ci[3](z2, 0) ~ 0,
            Ci[4](z2, 0) ~ 0,
        ]

        bcs_ics = [bcs; ics]

        @named pdesys = PDESystem(eq, bcs_ics, domains, [z1, z2, t], [
                C[1](z1, t), 
                C[2](z1, t), 
                C[3](z1, t), 
                C[4](z1, t), 
                Ci[1](z2, t), 
                Ci[2](z2, t), 
                Ci[3](z2, t), 
                Ci[4](z2, t), 
            ], 
            [
                ki[1] => 500000,
                Ea[1] => 50,
            ])
    end

    order = 2
    #discretization = MOLFiniteDifference([z1 => dz, z2 => dz], t, approx_order=order, grid_align=center_align)
    discretization = MOLFiniteDifference([z1 => dz_1, z2 => dz_2], t, approx_order=order, grid_align=center_align)

    # Convert the PDE problem into an ODE problem
    @time prob = discretize(pdesys, discretization);
    
    return prob, materials, auxmaterial

end

"""
    getReactionModel(catalyst)

returns the reactor model depending on the specific catalyst for a specific experimental run.
"""

#= returns a struct of SynthesisModel with odeprob, materials, auxmaterial, exdata
    the first three params are obtained by the other getreactionmodel function =#
        
function getreactionmodel(catalyst, Cin, filename::String)
    #filename = "C:/Users/peter.sagmeister/OneDrive - Research Center Pharmaceutical Engineering/ACS GCIPR_01/100_Data/Amidation/Julia/Julia_Kinetics/experiments/PS-GCIPR-27-20230220_V01.xlsx"

    data = importexperiment(filename);

    q = DataInterpolations.LinearInterpolation(vcat(data.flowrate_total...),  vcat(data.time_hitec...));
    Ti = DataInterpolations.LinearInterpolation(vcat(data.Tset...) .+ 273.15, vcat(data.time_hitec...)); 
    C_ester_in =   DataInterpolations.LinearInterpolation(vcat(data.C_inlet[1,:]...), vcat(data.time_hitec...));
    C_amine_in =   DataInterpolations.LinearInterpolation(vcat(data.C_inlet[2,:]...), vcat(data.time_hitec...));
    C_TBD_in =     DataInterpolations.LinearInterpolation(vcat(data.C_inlet[3,:]...), vcat(data.time_hitec...));
    C_product_in = DataInterpolations.LinearInterpolation(vcat(data.C_inlet[4,:]...), vcat(data.time_hitec...));
    C_intermediate_in = DataInterpolations.LinearInterpolation(vcat(data.C_inlet[5,:]...), vcat(data.time_hitec...));
    C_alcohol_in = DataInterpolations.LinearInterpolation(vcat(data.C_inlet[6,:]...), vcat(data.time_hitec...));

    tsim = data.tmax

    # C(t) = Cin    #mabye change to C = Cin ---> without C(t)
    #C = Cin    #mabye change to C = Cin ---> without C(t)
    # Cin(t) = [C_ester_in(t); C_amide_in(t); C_TBD_in(t); C_product_in(t)]    #maybe this works

    # prob, materials, auxmaterial = getreactionmodel(catalyst, C, Ti, q, tsim)
    prob, materials, auxmaterial = getreactionmodel(catalyst, C_ester_in, C_amine_in, C_TBD_in, C_product_in, C_intermediate_in, C_alcohol_in, Ti, q, tsim);
    return SynthesisModel(prob, materials, auxmaterial, data)

end


function getreactionmodelSS(catalyst, Cin::T, params) where {T<:Number}

    # reactor parameters
    V = 5.67E-3 # volume [dm^3]     # reactor parameters (for amidations)
   # V = 2.79E-3 # volume [dm^3]     # reactor parameters (for aAPI)
    numb_dz1 = 10 # number of discretization step
    L = V / 1 # L in [dm] ( --> A = 1 dm^2)
    Li = 0.2E-3 # inert length, second length after reactor   # you need atleast 5 compartments!
    numb_dz2 = 5 # number of discretization step
    dz = 1E-4 # discretization step size

    dz_1 = V / numb_dz1 # calculating discretization steps
    dz_2 = Li / numb_dz2 # calculating discretization steps

    tsim = 1000 # dummy simulation time

    @parameters t z1 z2

    Dt = Differential(t)
    Dz1 = Differential(z1)
    Dz2 = Differential(z2)

    if catalyst == "tworeactions"
        materials = ["ester"; "amine"; "TBD"; "product"; "intermediate"];
       # materials = ["ester"; "amine"; "TBD"; "product"];
       # materials = ["ester"; "amine"; "product"];
        if length(params) != 8
            @error "wrong length of parameter vector"
        end

        kr = params[1:4]
        ki = params[5:6]
        Ea = params[7:8]
       

        @parameters T0, q0, C_ester_in_0, C_amine_in_0, C_TBD_in_0
        @variables C[1:5](..) Ci[1:5](..)


        rall = [
            arrhenius2(C[1](z1, t), C[3](z1, t), kr[1], kr[3], ki[1], Ea[1], T0), # r1
            arrhenius2(C[2](z1, t), C[5](z1, t), kr[2], kr[4], ki[2], Ea[2], T0), # r2

        ]

              # r1  r2  
        B = [-1   0          # C1  Ester
              0  -1          # C2  amine    
             -1   1          # C3  TBD
              0   1          # C4  Product formation
              1  -1          # intermdeiate
        ]

        eq = [
            Dt(C[1](z1, t)) ~ -q0 * Dz1(C[1](z1, t)) + B[1, :]' * rall, # dC1/dt
            Dt(C[2](z1, t)) ~ -q0 * Dz1(C[2](z1, t)) + B[2, :]' * rall, # dC2/dt 
            Dt(C[3](z1, t)) ~ -q0 * Dz1(C[3](z1, t)) + B[3, :]' * rall, # dC3/dt 
            Dt(C[4](z1, t)) ~ -q0 * Dz1(C[4](z1, t)) + B[4, :]' * rall, # dC4/dt 
            Dt(C[5](z1, t)) ~ -q0 * Dz1(C[5](z1, t)) + B[5, :]' * rall, # dC5/dt 
            Dt(Ci[1](z2, t)) ~ -q0 * Dz2(Ci[1](z2, t)), # dC1/dt
            Dt(Ci[2](z2, t)) ~ -q0 * Dz2(Ci[2](z2, t)), # dC2/dt 
            Dt(Ci[3](z2, t)) ~ -q0 * Dz2(Ci[3](z2, t)), # dC3/dt 
            Dt(Ci[4](z2, t)) ~ -q0 * Dz2(Ci[4](z2, t)), # dC4/dt 
            Dt(Ci[5](z2, t)) ~ -q0 * Dz2(Ci[5](z2, t)), # dC5/dt 
        ]


        domains = [t ∈ Interval(0.0, tsim),
            z1 ∈ Interval(0.0, L),
            z2 ∈ Interval(L, L + Li)]

        # boundary conditions
        bcs = [
            C[1](0, t) ~ C_ester_in_0,
            C[2](0, t) ~ C_amine_in_0, 
            C[3](0, t) ~ C_TBD_in_0, 
            C[4](0, t) ~ 0,
            C[5](0, t) ~ 0,    #intermediate  set to small nonzero value (dirty bugfix), see https://github.com/SciML/MethodOfLines.jl/issues/194
            C[1](L, t) ~ Ci[1](L, t),   # interface BCs
            C[2](L, t) ~ Ci[2](L, t), # 
            C[3](L, t) ~ Ci[3](L, t), # 
            C[4](L, t) ~ Ci[4](L, t), # 
            C[5](L, t) ~ Ci[5](L, t), # intermediate
        ]

        ics = [
            C[1](z1, 0) ~ 0,         #initial conditions
            C[2](z1, 0) ~ 0,
            C[3](z1, 0) ~ 0,
            C[4](z1, 0) ~ 0,
            C[5](z1, 0) ~ 0,
            Ci[1](z2, 0) ~ 0,         #initial conditions
            Ci[2](z2, 0) ~ 0,
            Ci[3](z2, 0) ~ 0,
            Ci[4](z2, 0) ~ 0,
            Ci[5](z2, 0) ~ 0,
        ]

        bcs_ics = [bcs; ics]

        @named pdesys = PDESystem(eq, bcs_ics, domains, [z1, z2, t], [
                C[1](z1, t), 
                C[2](z1, t), 
                C[3](z1, t),
                C[4](z1, t), 
                C[5](z1, t),
                Ci[1](z2, t), 
                Ci[2](z2, t), 
                Ci[3](z2, t),
                Ci[4](z2, t), 
                Ci[5](z2, t)
            ],
            [
                T0 => 1.0
                q0 => 2.0
                C_ester_in_0 => 3.0 # 0.3  #--> maybe need to remove this
                C_amine_in_0 => 4.0 # 0.3  #--> maybe need to remove this
                C_TBD_in_0 => 5.0   # 0.1    #--> maybe need to remove this
            ])
    
    elseif catalyst == "tworeactions_fixed"
                materials = ["ester"; "amine"; "TBD"; "product"; "intermediate"];
               # materials = ["ester"; "amine"; "TBD"; "product"];
               # materials = ["ester"; "amine"; "product"];
                if length(params) != 4
                    @error "wrong length of parameter vector"
                end
        
               
                ki = params[1:2]
                Ea = params[3:4]
               
        
                @parameters T0, q0, C_ester_in_0, C_amine_in_0, C_TBD_in_0
                @variables C[1:5](..) Ci[1:5](..)
        
        
                rall = [
                    arrhenius2(C[1](z1, t), C[3](z1, t), 1, 1, ki[1], Ea[1], T0), # r1
                    arrhenius2(C[2](z1, t), C[5](z1, t), 1, 1, ki[2], Ea[2], T0), # r2
        
                ]
        
                      # r1  r2  
                B = [-1   0          # C1  Ester
                      0  -1          # C2  amine    
                     -1   1          # C3  TBD
                      0   1          # C4  Product formation
                      1  -1          # intermdeiate
                ]
        
                eq = [
                    Dt(C[1](z1, t)) ~ -q0 * Dz1(C[1](z1, t)) + B[1, :]' * rall, # dC1/dt
                    Dt(C[2](z1, t)) ~ -q0 * Dz1(C[2](z1, t)) + B[2, :]' * rall, # dC2/dt 
                    Dt(C[3](z1, t)) ~ -q0 * Dz1(C[3](z1, t)) + B[3, :]' * rall, # dC3/dt 
                    Dt(C[4](z1, t)) ~ -q0 * Dz1(C[4](z1, t)) + B[4, :]' * rall, # dC4/dt 
                    Dt(C[5](z1, t)) ~ -q0 * Dz1(C[5](z1, t)) + B[5, :]' * rall, # dC5/dt 
                    Dt(Ci[1](z2, t)) ~ -q0 * Dz2(Ci[1](z2, t)), # dC1/dt
                    Dt(Ci[2](z2, t)) ~ -q0 * Dz2(Ci[2](z2, t)), # dC2/dt 
                    Dt(Ci[3](z2, t)) ~ -q0 * Dz2(Ci[3](z2, t)), # dC3/dt 
                    Dt(Ci[4](z2, t)) ~ -q0 * Dz2(Ci[4](z2, t)), # dC4/dt 
                    Dt(Ci[5](z2, t)) ~ -q0 * Dz2(Ci[5](z2, t)), # dC5/dt 
                ]
        
        
                domains = [t ∈ Interval(0.0, tsim),
                    z1 ∈ Interval(0.0, L),
                    z2 ∈ Interval(L, L + Li)]
        
                # boundary conditions
                bcs = [
                    C[1](0, t) ~ C_ester_in_0,
                    C[2](0, t) ~ C_amine_in_0, 
                    C[3](0, t) ~ C_TBD_in_0, 
                    C[4](0, t) ~ 0,
                    C[5](0, t) ~ 0,    #intermediate  set to small nonzero value (dirty bugfix), see https://github.com/SciML/MethodOfLines.jl/issues/194
                    C[1](L, t) ~ Ci[1](L, t),   # interface BCs
                    C[2](L, t) ~ Ci[2](L, t), # 
                    C[3](L, t) ~ Ci[3](L, t), # 
                    C[4](L, t) ~ Ci[4](L, t), # 
                    C[5](L, t) ~ Ci[5](L, t), # intermediate
                ]
        
                ics = [
                    C[1](z1, 0) ~ 0,         #initial conditions
                    C[2](z1, 0) ~ 0,
                    C[3](z1, 0) ~ 0,
                    C[4](z1, 0) ~ 0,
                    C[5](z1, 0) ~ 0,
                    Ci[1](z2, 0) ~ 0,         #initial conditions
                    Ci[2](z2, 0) ~ 0,
                    Ci[3](z2, 0) ~ 0,
                    Ci[4](z2, 0) ~ 0,
                    Ci[5](z2, 0) ~ 0,
                ]
        
                bcs_ics = [bcs; ics]
        
                @named pdesys = PDESystem(eq, bcs_ics, domains, [z1, z2, t], [
                        C[1](z1, t), 
                        C[2](z1, t), 
                        C[3](z1, t),
                        C[4](z1, t), 
                        C[5](z1, t),
                        Ci[1](z2, t), 
                        Ci[2](z2, t), 
                        Ci[3](z2, t),
                        Ci[4](z2, t), 
                        Ci[5](z2, t)
                    ],
                    [
                        T0 => 1.0
                        q0 => 2.0
                        C_ester_in_0 => 3.0 # 0.3  #--> maybe need to remove this
                        C_amine_in_0 => 4.0 # 0.3  #--> maybe need to remove this
                        C_TBD_in_0 => 5.0   # 0.1    #--> maybe need to remove this
                    ])        

    elseif catalyst == "onereaction"

         materials = ["ester"; "amine"; "TBD"; "product"];
          if length(params) != 5
              @error "wrong length of parameter vector"
          end
  
          kr = params[1:3]
          ki = params[4]
          Ea = params[5]
            
          @parameters T0, q0, C_ester_in_0, C_amine_in_0, C_TBD_in_0 
          @variables C[1:5](..) Ci[1:5](..)
  
  
          rall = [
               arrhenius3(C[1](z1, t), C[2](z1, t), C[3](z1, t), kr[1], kr[2], kr[3], ki[1], Ea[1], T0), # r1
          #    arrhenius2(C[2](z1, t), C[5](z1, t), kr[2], kr[4], ki[2], Ea[2], T(t)), # r1
            ]
  
                # r1  
           B = [-1          # C1  Ester
                -1          # C2  amine    
                 0          # C3  TBD
                 1          # C4  Product formation       
           ]
           
  
          eq = [
              Dt(C[1](z1, t)) ~ -q0 * Dz1(C[1](z1, t)) + B[1, :]' * rall, # dC1/dt
              Dt(C[2](z1, t)) ~ -q0 * Dz1(C[2](z1, t)) + B[2, :]' * rall, # dC2/dt 
              Dt(C[3](z1, t)) ~ -q0 * Dz1(C[3](z1, t)) + B[3, :]' * rall, # dC3/dt 
              Dt(C[4](z1, t)) ~ -q0 * Dz1(C[4](z1, t)) + B[4, :]' * rall, # dC4/dt 
              Dt(Ci[1](z2, t)) ~ -q0 * Dz2(Ci[1](z2, t)), # dC1/dt
              Dt(Ci[2](z2, t)) ~ -q0 * Dz2(Ci[2](z2, t)), # dC2/dt 
              Dt(Ci[3](z2, t)) ~ -q0 * Dz2(Ci[3](z2, t)), # dC3/dt 
              Dt(Ci[4](z2, t)) ~ -q0 * Dz2(Ci[4](z2, t)), # dC4/dt 
          ]
  
  
          domains = [t ∈ Interval(0.0, tsim),
              z1 ∈ Interval(0.0, L),
              z2 ∈ Interval(L, L + Li)]
  
          # boundary conditions
          bcs = [
              C[1](0, t) ~ C_ester_in_0,   #this should be the initial concentration and be allowed to changed in the optimization --> changed this instead of  C[1](0, t) ~ C_ester_in(t), thi             
              C[2](0, t) ~ C_amine_in_0,   #this should be the initial concentration and be allowed to changed in the optimization --> changed this instead of  C[2](0, t) ~ C_amine_in(t),
              C[3](0, t) ~ C_TBD_in_0,     #this should be the initial concentration and be allowed to changed in the optimization --> changed this instead of     C[3](0, t) ~ C_TBD_in(t), 
              C[4](0, t) ~ 0,              #this should be the initial concentration and be allowed to changed in the optimization --> changed this instead of       C[4](0, t) ~ C_product_in(t),
              C[1](L, t) ~ Ci[1](L, t),    # interface BCs
              C[2](L, t) ~ Ci[2](L, t),    # 
              C[3](L, t) ~ Ci[3](L, t),    # 
              C[4](L, t) ~ Ci[4](L, t),    # 

          ]
  
          ics = [
              C[1](z1, 0) ~ 0,         #initial conditions
              C[2](z1, 0) ~ 0,
              C[3](z1, 0) ~ 0,
              C[4](z1, 0) ~ 0,
              Ci[1](z2, 0) ~ 0,         #initial conditions
              Ci[2](z2, 0) ~ 0,
              Ci[3](z2, 0) ~ 0,
              Ci[4](z2, 0) ~ 0,
          ]
  
          bcs_ics = [bcs; ics]
  
          @named pdesys = PDESystem(eq, bcs_ics, domains, [z1, z2, t], [
                  C[1](z1, t), 
                  C[2](z1, t), 
                  C[3](z1, t),
                  C[4](z1, t), 
                  Ci[1](z2, t), 
                  Ci[2](z2, t), 
                  Ci[3](z2, t),
                  Ci[4](z2, t), 
              ],
              [
                  T0 => 1.0  # ask markus why we have this here.
                  q0 => 2.0
                  #do i need to have here the inital concentrations?
                  C_ester_in_0 => 3.0 # 0.3  #--> maybe need to remove this
                  C_amine_in_0 => 4.0 # 0.3  #--> maybe need to remove this
                  C_TBD_in_0 => 5.0   # 0.1    #--> maybe need to remove this
                  
              ])

    

    elseif catalyst == "onereaction_fixed"

        materials = ["ester"; "amine"; "TBD"; "product"];
        if length(params) != 2
            @error "wrong length of parameter vector"
        end

       
        ki = params[1]
        Ea = params[2]
        
        @parameters T0, q0, C_ester_in_0, C_amine_in_0, C_TBD_in_0 
        @variables C[1:5](..) Ci[1:5](..)


        rall = [
            arrhenius3(C[1](z1, t), C[2](z1, t), C[3](z1, t), 1, 1, 1, ki[1], Ea[1], T0), # r1
        #    arrhenius2(C[2](z1, t), C[5](z1, t), kr[2], kr[4], ki[2], Ea[2], T(t)), # r1
        ]

            # r1  
        B = [-1          # C1  Ester
            -1          # C2  amine    
                0          # C3  TBD
                1          # C4  Product formation       
        ]
        

        eq = [
            Dt(C[1](z1, t)) ~ -q0 * Dz1(C[1](z1, t)) + B[1, :]' * rall, # dC1/dt
            Dt(C[2](z1, t)) ~ -q0 * Dz1(C[2](z1, t)) + B[2, :]' * rall, # dC2/dt 
            Dt(C[3](z1, t)) ~ -q0 * Dz1(C[3](z1, t)) + B[3, :]' * rall, # dC3/dt 
            Dt(C[4](z1, t)) ~ -q0 * Dz1(C[4](z1, t)) + B[4, :]' * rall, # dC4/dt 
            Dt(Ci[1](z2, t)) ~ -q0 * Dz2(Ci[1](z2, t)), # dC1/dt
            Dt(Ci[2](z2, t)) ~ -q0 * Dz2(Ci[2](z2, t)), # dC2/dt 
            Dt(Ci[3](z2, t)) ~ -q0 * Dz2(Ci[3](z2, t)), # dC3/dt 
            Dt(Ci[4](z2, t)) ~ -q0 * Dz2(Ci[4](z2, t)), # dC4/dt 
        ]


        domains = [t ∈ Interval(0.0, tsim),
            z1 ∈ Interval(0.0, L),
            z2 ∈ Interval(L, L + Li)]

        # boundary conditions
        bcs = [
            C[1](0, t) ~ C_ester_in_0,   #this should be the initial concentration and be allowed to changed in the optimization --> changed this instead of  C[1](0, t) ~ C_ester_in(t), thi             
            C[2](0, t) ~ C_amine_in_0,   #this should be the initial concentration and be allowed to changed in the optimization --> changed this instead of  C[2](0, t) ~ C_amine_in(t),
            C[3](0, t) ~ C_TBD_in_0,     #this should be the initial concentration and be allowed to changed in the optimization --> changed this instead of     C[3](0, t) ~ C_TBD_in(t), 
            C[4](0, t) ~ 0,              #this should be the initial concentration and be allowed to changed in the optimization --> changed this instead of       C[4](0, t) ~ C_product_in(t),
            C[1](L, t) ~ Ci[1](L, t),    # interface BCs
            C[2](L, t) ~ Ci[2](L, t),    # 
            C[3](L, t) ~ Ci[3](L, t),    # 
            C[4](L, t) ~ Ci[4](L, t),    # 

        ]

        ics = [
            C[1](z1, 0) ~ 0,         #initial conditions
            C[2](z1, 0) ~ 0,
            C[3](z1, 0) ~ 0,
            C[4](z1, 0) ~ 0,
            Ci[1](z2, 0) ~ 0,         #initial conditions
            Ci[2](z2, 0) ~ 0,
            Ci[3](z2, 0) ~ 0,
            Ci[4](z2, 0) ~ 0,
        ]

        bcs_ics = [bcs; ics]

        @named pdesys = PDESystem(eq, bcs_ics, domains, [z1, z2, t], [
                C[1](z1, t), 
                C[2](z1, t), 
                C[3](z1, t),
                C[4](z1, t), 
                Ci[1](z2, t), 
                Ci[2](z2, t), 
                Ci[3](z2, t),
                Ci[4](z2, t), 
            ],
            [
                T0 => 1.0  # ask markus why we have this here.
                q0 => 2.0
                #do i need to have here the inital concentrations?
                C_ester_in_0 => 3.0 # 0.3  #--> maybe need to remove this
                C_amine_in_0 => 4.0 # 0.3  #--> maybe need to remove this
                C_TBD_in_0 => 5.0   # 0.1    #--> maybe need to remove this
                
            ])

    end


    order = 2
   # discretization = MOLFiniteDifference([z1 => dz, z2 => dz], t, approx_order=order, grid_align=center_align)
    discretization = MOLFiniteDifference([z1 => dz_1, z2 => dz_2], t, approx_order=order, grid_align=center_align)

    # Convert the PDE problem into an ODE problem
    @time prob = discretize(pdesys, discretization)

    #return prob, materials


    return SynthesisModelSS(prob, materials)

end

function identifyparameters(model; starttime=400, endtime=10202, refine=false, w=[])

    prob = model.odeprob
    materials = model.materials
    exdata = model.exdata

    n = length(materials)
    p_init = prob.p
    np = length(p_init)
    # parameter bounds for global optimization
  #  lower = 1E-6 .* ones(np)      #change here the boundaries???
   # upper = 1E9 .* ones(np)       #change here the boundaries???
  
  # np_k = length([1])
  # np_EA = length([1])

    np_k = length([1;2])
    np_EA = length([1;2])

   # np_k = length([1;2;3])
   # np_EA = length([1;2;3])  

    lower_EA = 1 .* ones(np_EA)                   #change here the boundaries???
    upper_EA = 200 .* ones(np_EA)                    #change here the boundaries???
    lower_k = 1E-6 .* ones(np_k)                   #change here the boundaries???
    upper_k = 1E9 .* ones(np_k)                    #change here the boundaries???    
    lower = [lower_k; lower_EA]                     #change here the boundaries???
    upper = [upper_k; upper_EA]                     #change here the boundaries???

    #hard coded boundaries for global optimization   ---> would need to find a way to better set the values ? pass them on in the function?
  #  np_rate = length([1;2;3])
  #  np_k = length([1])
  #  np_EA = length([1])
  #  lower_rate = 0 .* ones(np_rate)      #change here the boundaries???
  #  upper_rate = 3 .* ones(np_rate)       #change here the boundaries???
  #  lower_EA = 1E-1 .* ones(np_k)         #change here the boundaries???
  #  upper_EA = 1E3 .* ones(np_k)          #change here the boundaries???
  #  lower_k = 1E-6 .* ones(np_EA)         #change here the boundaries???
  #  upper_k = 1E9 .* ones(np_EA)          #change here the boundaries???    
  #  lower = [lower_rate; lower_k; lower_EA]      #change here the boundaries???
  #  upper = [upper_rate; upper_k; upper_EA]        #change here the boundaries???
   
  #  if reactionsystem == "tworeactions"
  #      np_rate = length([1;2;3;4])
  #      np_k = length([1;2])
  #      np_EA = length([1;2])
  #      lower_rate = 0 .* ones(np_rate)                 #change here the boundaries???
  #      upper_rate = 3 .* ones(np_rate)                 #change here the boundaries???
  #      lower_EA = 1E-1 .* ones(np_k)                   #change here the boundaries???
  #      upper_EA = 1E3 .* ones(np_k)                    #change here the boundaries???
  #      lower_k = 1E-6 .* ones(np_EA)                   #change here the boundaries???
  #      upper_k = 1E9 .* ones(np_EA)                    #change here the boundaries???    
  #      lower = [lower_rate;lower_k; lower_EA]          #change here the boundaries???
  #      upper = [upper_rate;upper_k; upper_EA]          #change here the boundaries???
  
  #  elseif reactionsystem == "tworeactions_fixed"
  #      np_k = length([1;2])
  #      np_EA = length([1;2])
  #      lower_EA = 1E-1 .* ones(np_k)                   #change here the boundaries???
  #      upper_EA = 1E3 .* ones(np_k)                    #change here the boundaries???
  #      lower_k = 1E-6 .* ones(np_EA)                   #change here the boundaries???
  #      upper_k = 1E9 .* ones(np_EA)                    #change here the boundaries???    
  #      lower = [lower_k; lower_EA]                     #change here the boundaries???
  #      upper = [upper_k; upper_EA]                     #change here the boundaries???
   # end

 #  np_k = length([1])
 #  np_EA = length([1])
 #  lower_EA = 1E-1 .* ones(np_k)                   #change here the boundaries???
 #  upper_EA = 1E3 .* ones(np_k)                    #change here the boundaries???
 #  lower_k = 1E-6 .* ones(np_EA)                   #change here the boundaries???
 #  upper_k = 1E9 .* ones(np_EA)                    #change here the boundaries???    
 #  lower = [lower_k; lower_EA]                     #change here the boundaries???
 #  upper = [upper_k; upper_EA]                     #change here the boundaries???


    bounds = [lower upper]  
 # original   startidx = findfirst(x -> x >= starttime * 3600, exdata.time_ftir).I[1]
 # original   endidx = findfirst(x -> x >= endtime * 3600, exdata.time_ftir).I[1]
    startidx = findfirst(x -> x >= starttime, exdata.time_ftir).I[1]   #removed 3600 conversion to hours?
    endidx = findfirst(x -> x >= endtime, exdata.time_ftir).I[1]       #removed 3600 conversion to hours?
   
    time_ident = exdata.time_ftir[startidx:endidx]

    if model.auxmaterial
        idx = broadcast(i -> matdict[i], materials[1:end-1]) 
        data = exdata.c_ftir[idx, :]
        sum_ftir = exdata.c_ftir[matdict["sum"],:]
        sum_data = sum(data,dims=1)

        aux = sum_ftir-sum_data[:]
        data = [data;aux']
    else
        idx = broadcast(i -> matdict[i], materials)
        data = exdata.c_ftir[idx, :]
    end
    
    data = data[:, startidx:endidx]

    if isempty(w)                       #why do we need w? --> ask markus
        w = ones(n)
    end

    if length(w) != n
        @error "incorrect length of weight vector!"
    end

    weight = repeat(w, 1, length(time_ident))

    z_outlet =  prob.problem_type.pdesys.domain[3].domain.right # --> right boundary

    # two parameters needed in cost function, adtype is just a dummy for now...
    cost(p_eval, adtype = SciMLBase.NoAD()) = cost_fun(p_eval,prob,time_ident,data,w, z_outlet)
    
    @info "Starting parameter estimation with NLopt-BOBYQA"
    opt = NLopt.Opt(:LN_BOBYQA, np)
    #opt = NLopt.Opt(:LN_COBYLA, np)
    optprob = Optimization.OptimizationProblem(cost, p_init, lb=lower, ub=upper)
   
    @time minx = solve(optprob, opt)
    @info "Parameters NLopt-BOBYQA"
    @info minx

    red_lower_k1 = minx[1] .- (minx[1]*0.05)                  #change here the boundaries???
    red_upper_k1 = minx[1] .+ (minx[1]*0.05)                    #change here the boundaries???   
    red_lower_k2 = minx[2] .- (minx[2]*0.05)                   #change here the boundaries???
    red_upper_k2 = minx[2] .+ (minx[2]*0.05)                    #change here the boundaries???   
  #  red_lower_k2 = minx[3] .- (minx[3]*0.05)                   #change here the boundaries???
   # red_upper_k2 = minx[3] .+ (minx[3]*0.05)                    #change here the boundaries???   

   # red_lower_EA1 = minx[2] .- (minx[2]*0.05)                  #change here the boundaries???
   # red_upper_EA1 = minx[2] .+ (minx[2]*0.05)                    #change here the boundaries???
    red_lower_EA1 = minx[3] .- (minx[3]*0.05)                  #change here the boundaries???
    red_upper_EA1 = minx[3] .+ (minx[3]*0.05)                    #change here the boundaries???
    red_lower_EA2 = minx[4] .- (minx[4]*0.05)                      #change here the boundaries???
    red_upper_EA2 = minx[4] .+ (minx[4]*0.05)                        #change here the boundaries??? 
 #   red_lower_EA2 = minx[6] .- (minx[6]*0.05)                      #change here the boundaries???
 #   red_upper_EA2 = minx[6] .+ (minx[6]*0.05)                        #change here the boundaries??? 
  
 #  red_lower_kr1 = minx[1] .- (minx[1]*0.05)                  #change here the boundaries???
  # red_upper_kr1 = minx[1] .+ (minx[1]*0.05)                    #change here the boundaries???   
 #  red_lower_kr2 = minx[2] .- (minx[2]*0.05)                   #change here the boundaries???
 #  red_upper_kr2 = minx[2] .+ (minx[2]*0.05)                    #change here the boundaries???   
 #  red_lower_kr3 = minx[3] .- (minx[3]*0.05)                   #change here the boundaries???
 #  red_upper_kr3 = minx[3] .+ (minx[3]*0.05)                    #change here the boundaries???   
 #  red_lower_k1 = minx[4] .- (minx[4]*0.05)                  #change here the boundaries???
 #  red_upper_k1 = minx[4] .+ (minx[4]*0.05)                    #change here the boundaries???
 #  red_lower_EA1 = minx[5] .- (minx[5]*0.05)                  #change here the boundaries???
 #  red_upper_EA1 = minx[5] .+ (minx[5]*0.05)                    #change here the boundaries???
  # red_lower_EA2 = minx[4] .- (minx[4]*0.05)                      #change here the boundaries???
  # red_upper_EA2 = minx[4] .+ (minx[4]*0.05)                        #change here the boundaries???

 #   reduced_lower = [red_lower_kr1; red_lower_kr2; red_lower_kr3; red_lower_k1; red_lower_EA1]                     #change here the boundaries???
 #   reduced_upper = [red_upper_kr1; red_upper_kr2; red_upper_kr3; red_upper_k1; red_upper_EA1]                     #change here the boundaries???
    
    reduced_lower = [red_lower_k1; red_lower_k2; red_lower_EA1; red_lower_EA2]                     
    reduced_upper = [red_upper_k1; red_upper_k2; red_upper_EA1; red_upper_EA2] 

   # reduced_lower = [lower_k; lower_EA]                     #change here the boundaries???
   # reduced_upper = [upper_k; upper_EA]                     #change here the boundaries???
    
   # reduced_lower = [red_lower_k1; red_lower_EA1]                     #change here the boundaries???
  #  reduced_upper = [red_upper_k1; red_upper_EA1]                     #change here the boundaries???

    if (refine)
        
        #cost(p_eval, adtype = SciMLBase.NoAD()) = cost_fun(p_eval,prob,time_ident,data,weight, z_outlet)
        @info "Starting parameter refinement with Nelder-Mead"
        opt2 = NLopt.Opt(:LN_NELDERMEAD, np)
        optprob2 = Optimization.OptimizationProblem(cost, minx, lb=reduced_lower, ub=reduced_upper)
      #  optprob2 = Optimization.OptimizationProblem(cost, minx, lb=lower, ub=upper)

        @time minx = solve(optprob2, opt2) # TODO: what does solve ? (solve! is empty function)

    end
    @info "Parameters Nelder Mead"
    @info minx
    return minx.u

end

function identifyparametersnew(reactionsystem, algorithm, model; w=[])
    starttime=60
    endtime=15588

    prob = model.odeprob
    materials = model.materials
    exdata = model.exdata

    n = length(materials)
    p_init = prob.p
    np = length(p_init)
    # parameter bounds for global optimization
    #lower = 1E-2 .* ones(np)      #change here the boundaries???
    #upper = 1E2 .* ones(np)       #change here the boundaries???
   
    #hard coded boundaries for global optimization   ---> would need to find a way to better set the values ? pass them on in the function?
    if reactionsystem == "tworeactions"
        np_rate = length([1;2;3;4])
        np_k = length([1;2])
        np_EA = length([1;2])
        lower_rate = 0 .* ones(np_rate)                 #change here the boundaries???
        upper_rate = 3 .* ones(np_rate)                 #change here the boundaries???
        lower_EA = 1E-1 .* ones(np_k)                   #change here the boundaries???
        upper_EA = 1E3 .* ones(np_k)                    #change here the boundaries???
        lower_k = 1E-6 .* ones(np_EA)                   #change here the boundaries???
        upper_k = 1E9 .* ones(np_EA)                    #change here the boundaries???    
        lower = [lower_rate;lower_k; lower_EA]          #change here the boundaries???
        upper = [upper_rate;upper_k; upper_EA]          #change here the boundaries???
    elseif reactionsystem == "tworeactions_fixed"
        np_k = length([1;2])
        np_EA = length([1;2])
        lower_EA = 1E-1 .* ones(np_k)                   #change here the boundaries???
        upper_EA = 1E3 .* ones(np_k)                    #change here the boundaries???
        lower_k = 1E-6 .* ones(np_EA)                   #change here the boundaries???
        upper_k = 1E9 .* ones(np_EA)                    #change here the boundaries???    
        lower = [lower_k; lower_EA]                     #change here the boundaries???
        upper = [upper_k; upper_EA]                     #change here the boundaries???
    end
    bounds = [lower upper]

    startidx = findfirst(x -> x >= starttime, exdata.time_ftir).I[1]   
    endidx = findfirst(x -> x >= endtime, exdata.time_ftir).I[1]       
   
    time_ident = exdata.time_ftir[startidx:endidx]

    if model.auxmaterial
        idx = broadcast(i -> matdict[i], materials[1:end-1]) 
        data = exdata.c_ftir[idx, :]
        sum_ftir = exdata.c_ftir[matdict["sum"],:]
        sum_data = sum(data,dims=1)

        aux = sum_ftir-sum_data[:]
        data = [data;aux']
    else
        idx = broadcast(i -> matdict[i], materials)
        data = exdata.c_ftir[idx, :]
    end
    
    data = data[:, startidx:endidx]

    if isempty(w)                       #why do we need w? --> ask markus
        w = ones(n)
    end

    if length(w) != n
        @error "incorrect length of weight vector!"
    end

    weight = repeat(w, 1, length(time_ident))

    z_outlet =  prob.problem_type.pdesys.domain[3].domain.right # --> right boundary

    # two parameters needed in cost function, adtype is just a dummy for now...
    cost(p_eval, adtype = SciMLBase.NoAD()) = cost_fun(p_eval,prob,time_ident,data,w, z_outlet)
 
    if algorithm == "ECA"
        @info "Starting parameter estimation with ECA algorithm"
        @time result = Metaheuristics.optimize(cost, bounds, ECA())
        elseif algorithm == "DE"
            @info "Starting parameter estimation with DE algorithm"
            @time result = Metaheuristics.optimize(cost, bounds, DE())        
        elseif algorithm == "PSO"
            @info "Starting parameter estimation with PSO algorithm"
            @time result = Metaheuristics.optimize(cost, bounds, PSO(N = 100, C1=2, C2=2, ω = 0.8))
        elseif algorithm == "ABC"
            @info "Starting parameter estimation with ABC algorithm"
            @time result = Metaheuristics.optimize(cost, bounds, ABC(N = 80, No=20, Ne=50, limit = 5))
        elseif algorithm == "CGSA"     
            @info "Starting parameter estimation with CGSA algorithm"
            @time result = Metaheuristics.optimize(cost, bounds, CGSA())        
        elseif algorithm == "WOA"     
            @info "Starting parameter estimation with WOA algorithm"
            @time result = Metaheuristics.optimize(cost, bounds, WOA())   
        elseif algorithm == "NSGA2"     
            @info "Starting parameter estimation with NSGA2 algorithm"
          #  f(p_eval) = (cost(p_eval), [0.0], [0.0])

            function fcost(p_eval,prob,time_ident,data,w, z_outlet)
                fx = [cost_fun(p_eval,prob,time_ident,data,w, z_outlet); 0.0]
                gx = [0.0]
                hx = [0.0]
                # order is important
                return fx, gx, hx
            end

            f(p_eval) = fcost(p_eval,prob,time_ident,data,w, z_outlet)

            @time result = Metaheuristics.optimize(f, bounds, NSGA2())   
        elseif algorithm == "NSGA3"     
            @info "Starting parameter estimation with NSGA3 algorithm"
            f = (cost, [0.0], [0.0])
            @time result = Metaheuristics.optimize(f, bounds, NSGA3())  
        elseif algorithm == "SMS-EMOA"     
            @info "Starting parameter estimation with SMS-EMOA algorithm"
            @time result = Metaheuristics.optimize(cost, bounds, SMS_EMOA())   
        elseif algorithm == "SPEA2"     
            @info "Starting parameter estimation with SPEA2 algorithm"
            @time result = Metaheuristics.optimize(cost, bounds, SPEA2())   
        elseif algorithm == "MCCGA"     
            @info "Starting parameter estimation with MCCGA algorithm"
            @time result = Metaheuristics.optimize(cost, bounds, MCCGA())   
        elseif algorithm == "GA"     
            @info "Starting parameter estimation with GA algorithm"
            @time result = Metaheuristics.optimize(cost, bounds, GA())   
        elseif algorithm == "CCMO"     
            @info "Starting parameter estimation with CCMO algorithm"
            @time result = Metaheuristics.optimize(cost, bounds, CCMO())       
        elseif algorithm == "εDE"     
            @info "Starting parameter estimation with εDE algorithm"
            @time result = Metaheuristics.optimize(cost, bounds, εDE())     
        elseif algorithm == "BRKGA"     
            @info "Starting parameter estimation with BRKGA algorithm"
            @time result = Metaheuristics.optimize(cost, bounds, BRKGA())    
    end 
    
    # Common options
  #  options = Metaheuristics.Options(seed=1, store_convergence = true)
    # Optimize
   # @time result = Metaheuristics.optimize(cost, bounds, NSGA2(options=options))
   # @time result = Metaheuristics.optimize(cost, bounds, PSO(N = 10, C1=2, C2=2, ω = 0.8))
 #   @time result = Metaheuristics.optimize(cost, bounds, ABC(N = 80, No=20, Ne=50, limit = 5))

    return result

end

function cost_fun(p_eval,prob,time_ident,data,weight,z_outlet)       
    tmp_prob = remake(prob,p = p_eval)
    sol = solve(tmp_prob,TRBDF2(),saveat = time_ident)

    cost = 0.0

    if sol.retcode != ReturnCode.Success
        cost = Inf
    else

    #retrieve solution at the reactor outlet
    sol_val = sol(time_ident,0.0,z_outlet)

    n = length(weight)# only the concentrations of the inert part are needed

    cost = 0.0

    for i = 1:n
        
        simres = sol_val[n+i]
        datares = data[i,:]

        for j = 1:length(simres)
            
            cost = cost + (simres[j]-datares[j])^2*weight[i]

        end

    end

    end
    return cost

end

model_ode(model, p_) = remake(model.odeprob, p=p_);
#solve_model(model, mp_) = solve(model_ode(model, mp_), TRBDF2(), wrap=Val(false), save_idxs=OptSynthesis.outletidx[length(model.materials)])

function solve_model(model, params)

    newprob = remake(model.odeprob, p=params)
    sol = solve(newprob,TRBDF2(),saveat = model.exdata.time_hitec)

    if sol.retcode != ReturnCode.Success
        @warn "Problem solving the model..."
    end

    z_outlet = newprob.problem_type.pdesys.domain[3].domain.right # --> right boundary
    t = vec(model.exdata.time_hitec)
    #retrieve solution at the reactor outlet
    sol_val = sol(t,0.0,z_outlet)

    n = length(model.materials)
    return t,sol_val[n+1:2*n]

end


# solve the model

function steadystate(model::SynthesisModelSS, params; tend=1500)

    newode = remake(model.odeprob, p=params, tspan=(0.0, tend))
    ss = solve(newode, TRBDF2(), saveat=tend)

    if ss.retcode != ReturnCode.Success
        @warn "Problem solving the model..."
    end

    z_outlet = newode.problem_type.pdesys.domain[3].domain.right # --> right boundary
    sol_val = ss(tend,0.0,z_outlet)

    n = length(model.materials) # number of materials
    return sol_val[n+1:2*n]

end


function plotidentificationresults(model, t, sol_model)

    time_exp = model.exdata.time_ftir

    exdata = model.exdata
    materials = model.materials

    if model.auxmaterial
        idx = broadcast(i -> matdict[i], materials[1:end-1]) 
        data = exdata.c_ftir[idx, :]
        sum_ftir = exdata.c_ftir[matdict["sum"],:]
        sum_data = sum(data,dims=1)

        aux = sum_ftir-sum_data[:]
        data = [data;aux']
    else
        idx = broadcast(i -> matdict[i], materials)
        data = exdata.c_ftir[idx, :]
    end

    # #data_exp = model.exdata.c_ftir

    # idx = broadcast(i -> matdict[i], model.materials)  
    
    n = length(model.materials)

    emptystring = repeat([""]; outer=[1, n - 1])

    #plot is a Plots.jl function
    p = plot(time_exp / 3600, data', layout=(n, 1), ylabel=permutedims(model.materials),
        label=hcat("meas.", emptystring),
        linewidth=2,
        xlims=(0, 6.7)
    )
    plot!(p, t / 3600, hcat(sol_model...), layout=(n, 1),
        size=(600, 200 * n),
        left_margin=15mm,
        yguidefontsize=8,
        yformatter=:scientific,
        label=hcat("sim.", emptystring),
        xlabel=hcat(emptystring, "Time t [h]"),
        linewidth=2,
        xlims=(0, 6.7),
        legend=:outertopright,
    )

    return p
end