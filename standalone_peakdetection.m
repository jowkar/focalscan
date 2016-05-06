function t = standalone_peakdetection(report_file_path,annot_file_path,peak_level,scorefield,out_file,varargin)
    p = inputParser;
    p.addRequired('report_file_path',@isstr)
    p.addRequired('annot_file_path',@isstr)
    p.addRequired('peak_level',@isnumeric)
    p.addRequired('scorefield',@isstr)
    p.addRequired('out_file',@isstr)
    p.parse(report_file_path,annot_file_path,peak_level,scorefield,out_file);

    if ischar(peak_level)
        peak_level = str2double(peak_level);
    end
    if peak_level > 1 || peak_level < 0
        error('peak_level should be a number between 0 and 1')
    end

    if ~ismember(scorefield,{'fs','fs_hp','sum_cna','sum_cna_hp','spearman_corr'})
        error('Invalid value for the "scorefield" parameter')
    end

    report = readtable(report_file_path,'Delimiter','\t','ReadVariableNames',1,'ReadRowNames',0);
    annot = Annot(annot_file_path);
    T = FocalScan.make_table(report,annot,scorefield);

    t = FocalScan.make_peak_table(T,peak_level,varargin{:});

    if ~isempty(out_file)
        pth = strsplit(out_file,{'/',filesep});
        if length(pth) > 1
            pth = strjoin(pth(1:end-1),filesep);
            mkdir(pth)
        end
        writetable(t,out_file,'Delimiter','\t','FileType','text');
        fprintf('Successfully wrote peak table to %s\n',out_file)
    end
end