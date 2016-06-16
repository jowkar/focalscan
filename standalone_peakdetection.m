function t = standalone_peakdetection(report_file_path,annot_file_path,peak_level,scorefield,out_file,varargin)

    % STANDALONE_PEAKDETECTION: find peak genes/tiles based on the report created by FocalScan
    %
    % USAGE:
    %     standalone_peakdetection(report_file_path,annot_file_path,peak_level,scorefield,out_file)
    % 
    %         report_file_path: path to main report file created by running FocalScan (ie. report.txt)
    %         annot_file_path: path to annotation file (ie. annotation/gencode17.bed, gene/tile ids must match those in the report file)
    %         peak_level: level of granularity at which to find peaks (0.1-1, where 1 is least granular)
    %         scorefield: metric in the report file to be used as basis for peak detection;
    %             (Valid options:
    %                 fs_hp: the standard FocalScan score (with focality filter)
    %                 fs: FocalScan score without focality filter
    %                 sum_cna_hp: summed copy number amplitudes, with focality filter
    %                 sum_cna: summed copy number amplitudes, without focality filter
    %                 pearson_corr: pearson correlation coefficient)
    %         out_file: name of output file (ie. peaks.txt)
    % 
    % EXAMPLE:
    %     standalone_peakdetection.sh report.txt annotation/gencode17.bed 0.7 fs_hp ./peaks.txt

    isnumericstr = @(x)isnumeric(x)|isnumeric(str2double(x));

    p = inputParser;
    p.addRequired('report_file_path',@isstr)
    p.addRequired('annot_file_path',@isstr)
    p.addRequired('peak_level',isnumericstr)
    expected_values = {'fs','fs_hp','sum_cna','sum_cna_hp','pearson_corr'};
    %expected_values = {'fs','fs_hp','sum_cna','sum_cna_hp','spearman_corr'};
    p.addRequired('scorefield',@(x) any(validatestring(x,expected_values)));
    p.addRequired('out_file',@isstr);
    p.addParameter('optional_gene_annot','',@isstr);
    
    p.addParameter('writedir', '', @isstr);
    p.addParameter('do_plot', 0, isnumericstr);
    p.addParameter('min_genes', 100, isnumericstr);
    p.addParameter('plot_dir', '', @isstr);
    
    try
        p.parse(report_file_path,annot_file_path,peak_level,scorefield,out_file,varargin{:});
    catch ME
%         if isdeployed
%             type help_standalone_peakdetection
%         else
%             help standalone_peakdetection.m
%         end
        disp('Invalid input. Run the program without arguments to display a help page (or use the help command within the MATLAB environment).')
        rethrow(ME)
    end

    peak_level = numstr2num(p.Results.peak_level);
    optional_gene_annot = p.Results.optional_gene_annot;
    writedir = p.Results.writedir;
    do_plot = numstr2num(p.Results.do_plot);
    min_genes = numstr2num(p.Results.min_genes);
    plot_dir = p.Results.plot_dir;
    
    if peak_level > 1 || peak_level < 0
        error('peak_level should be a number between 0 and 1')
    end

    report = readtable(report_file_path,'Delimiter','\t','ReadVariableNames',1,'ReadRowNames',0,'FileType','text');
    annot = Annot(annot_file_path);
    T = FocalScan.make_table(report,annot,scorefield);
       
    t = FocalScan.make_peak_table(T,peak_level,'writedir',writedir,'do_plot',do_plot,'min_genes',min_genes,'plot_dir',plot_dir);
    
    if ~isempty(optional_gene_annot);
        t = FocalScan.get_tile_gene_ids(t,optional_gene_annot);
    end

    if ~isempty(out_file)
        pth = strsplit(out_file,filesep);
        if length(pth) > 1
            pth = strjoin(pth(1:end-1),filesep);
            mkdir(pth)
        end
        writetable(t,out_file,'Delimiter','\t','FileType','text','WriteVariableNames',1,'WriteRowNames',0);
        fprintf('Successfully wrote peak table to %s\n',out_file)
    end
    
    function num = numstr2num(numstr)
        if isnumeric(numstr)
            num = numstr;
        elseif isnumeric(str2double(numstr))
            num = str2double(numstr);
        else
            error('At least one numeric input parameter was not numeric.')
        end
    end
end
