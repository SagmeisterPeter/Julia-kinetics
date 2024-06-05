# Performance indicator functions
conversion(Ca, Ca0) = (Ca0-Ca)/Ca0
selectivity(C, Ca, Ca0) = C/(Ca0-Ca)    # selectivity, where Ca is the starting material concentration w. input conc. Ca0
productivity(C, q, M) = q*M*C .*60 .*1E3 .*60 ./1E3         # C,q,M are concentration, flowrate & molar Measuresv   l/s
spacetimeyield(C, q, M, V0) = productivity(C, q, M)/V0


