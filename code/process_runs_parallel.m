function process_runs_parallel(runs)
% processes all runs, defined by a struct array, in parallel

parfor crun = 1:length(runs)
    process_run(runs(crun));
end

end