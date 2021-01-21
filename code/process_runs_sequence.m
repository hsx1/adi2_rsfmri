function process_runs_sequence(runs)
% processes all runs, defined by a struct array, in sequence (not in parallel)

for crun = runs
    process_run(crun);
end

end