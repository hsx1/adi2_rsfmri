function process_run(crun)
% one estimation at a time

if crun.TYPE == "GT"
    RunModelGroupTime(crun)
elseif crun.TYPE == "BMI"
    RunModelBMI(crun)
elseif crun.TYPE == "FD"
    RunModelFD(crun)
else
    error("Invalid specification!")

end