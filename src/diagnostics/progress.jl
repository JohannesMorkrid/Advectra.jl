# ------------------------------------------------------------------------------------------
#                                    Progress Diagnostic                                    
# ------------------------------------------------------------------------------------------

"""
    progress_bar(state, prob, time, progress; t0, dt)

    Update `progress_bar` (`Progress` type from `ProgressMeter`.jl) based on evolution `time`.

    ### Keyword arguments:
    - `t0`: start time (Number)
    - `dt`: timestep (Number)
"""
function progress(state, prob, time, progress_bar; t0, dt)
    update!(progress_bar, round(Int, (time - t0) / dt))
end

function build_diagnostic(::Val{:progress}; schedule=Every(1), tspan, dt, kwargs...)
    args = (Progress(floor(Int, (last(tspan) - first(tspan)) / dt); showspeed=true),)
    kwargs = (; t0=first(tspan), dt=dt)
    Diagnostic(; name="Progress",
               method=progress,
               schedule=schedule,
               metadata="Progressbar.",
               assumes_spectral_state=true,
               stores_data=false,
               args=args,
               kwargs=kwargs)
end
