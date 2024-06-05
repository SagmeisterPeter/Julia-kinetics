function loss_abs(sol; ident_data = data, weight = w)
    tot_loss = 0.0

    if any((s.retcode != :Success for s in sol))
        tot_loss = Inf
    else
        tmp = weight.*(hcat(sol.u...) - ident_data)
        tot_loss = sum(abs.(tmp[:]))
    end
    tot_loss
end
