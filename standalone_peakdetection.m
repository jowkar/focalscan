function t = standalone_peakdetection(report_file_path,annot_file_path,peak_level,scorefield,out_file,varargin)
%     try
        if ischar(peak_level)
            peak_level = str2num(peak_level);
        end
        if peak_level > 1 || peak_level < 0
            error('peak_level should be a number between 0 and 1')
        end

        if ~ismember(scorefield,{'fs','fs_hp','sum_cna','sum_cna_hp','spearman_corr'})
            error('Invalid value for the "scorefield" parameter')
        end

        report = readtable(report_file_path,'Delimiter','\t','ReadVariableNames',1);
        report(1:10,:)
        annot = Annot(annot_file_path);
        T = FocalScan.make_table(report,annot,scorefield);

        t = FocalScan.make_peak_table(T,peak_level,varargin{:});

        if ~isempty(out_file)
            pth = strsplit(out_file,{'/',filesep});
            if length(pth) > 1
                pth = strjoin(pth(1:end-1),filesep);
                mkdir(pth)
            end
            writetable(t,out_file,'Delimiter','\t');
        end
%     catch ME
%         disp(getReport(ME,'extended','hyperlinks','off'))
%         exit;
%     end
end
