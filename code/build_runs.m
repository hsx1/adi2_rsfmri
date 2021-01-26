function [runs] = build_runs(param)
% builds structs for each estimation

runs = [];
% loop over all models, only all covariates, all inference types (for
% wild_boot = true) and all ROIs
for l = 1:length(param.MODEL)
    for j = 1:length(param.COVARIATES)
        % skip invalid Model/Covariate combinations
        ccov = param.COVARIATES(j);
        cmodel = param.MODEL(l);

        if startsWith(cmodel,"grouptime") 
            ctype = "GT";
        elseif startsWith(cmodel,"bmi") 
            ctype = "BMI";
        elseif startsWith(cmodel,"fd") 
            ctype = "FD";
        elseif startsWith(cmodel,"all") 
            ctype = "ALL";
        elseif startsWith(cmodel,"single") 
            ctype = "SINGLE";
        else
            error("invalide model name")
        end

        if ((strcmp(ctype,"GT") && (ccov < 10 || ccov >= 20) ) || ...
                (strcmp(ctype,"BMI") && (ccov < 20 || ccov >= 30))|| ...
                (strcmp(ctype,"FD") && (ccov < 30 || ccov >= 40)) || ...
                ((strcmp(ctype,"ALL") || strcmp(ctype,"SINGLE")) && (ccov < 40 || ccov >= 50)) ) 
            continue;
        end

        for k = 1:length(param.INFERENCE_TYPE)
            for m = 1:length(param.ROI_PREP)
                
                crun.ONLY_DISPLAY = param.ONLY_DISPLAY;
                crun.VIEW = param.VIEW;
                crun.OVERWRITE = param.OVERWRITE;
                crun.OUT_DIR = param.OUT_DIR;
                crun.INFO_DIR = param.INFO_DIR;
                crun.MASK_DIR = param.MASK_DIR;
                crun.MASK_GM = param.MASK_GM;
                crun.MASK_B = param.MASK_B;
                crun.MASK = param.MASK; 
                crun.EXCLFD = param.EXCLFD;
                crun.WILD_BOOT = param.WILD_BOOT;
                if ~crun.WILD_BOOT 
                    crun.INFERENCE_TYPE = "parametric"; 
                else
                    crun.INFERENCE_TYPE = param.INFERENCE_TYPE;
                end
                crun.TYPE = ctype;
                crun.MODEL = param.MODEL(l);
                crun.COVARIATES = param.COVARIATES(j);
                crun.INFERENCE_TYPE = param.INFERENCE_TYPE(k);
                crun.ROI_PREP = param.ROI_PREP(m);

                runs = [runs, crun];
            end
        end
    end
end

end