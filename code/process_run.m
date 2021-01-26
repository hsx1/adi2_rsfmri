function process_run(crun)
% calls model specific function to perform one estimation at a time.
% called by run_or_display.m

if crun.TYPE == "GT"
    RunModelGroupTime(crun)
elseif crun.TYPE == "BMI"
    RunModelBMI(crun)
elseif crun.TYPE == "FD"
    RunModelFD(crun)
elseif crun.TYPE == "ALL"
    RunModelAllTp(crun)
elseif crun.TYPE == "SINGLE"
    RunModelSingleTp(crun)
else
    error("Invalid specification!")

end