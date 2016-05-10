function t = standalone_peakdetection(report_file_path,annot_file_path,peak_level,scorefield,out_file,varargin)
    p = inputParser;
    p.addRequired('report_file_path',@isstr)
    p.addRequired('annot_file_path',@isstr)
    p.addRequired('peak_level',@isnumeric)
    expected_values = {'fs','fs_hp','sum_cna','sum_cna_hp','pearson_corr'};
    p.addRequired('scorefield',@(x) any(validatestring(x,expected_values)));
    p.addRequired('out_file',@isstr);
    p.addParameter('optional_gene_annot','',@isstr);
    
    p.addParameter('writedir', '', @isstr);
    p.addParameter('do_plot', 0, @isnumeric);
    p.addParameter('min_genes', 100, @isnumeric);
    p.addParameter('plot_dir', '', @isstr);
    
    p.parse(report_file_path,annot_file_path,peak_level,scorefield,out_file,varargin{:});

    optional_gene_annot = p.Results.optional_gene_annot;
    writedir = p.Results.writedir;
    do_plot = p.Results.do_plot;
    min_genes = p.Results.min_genes;
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
end