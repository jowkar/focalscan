classdef FocalScan
    % FOCALSCAN Rank genes based on coordinated focal changes in copy
    % number and expression
    %
    %             Parameter:                    Default:    Type:
    %
    %   [Basic options]
    %             'expr_csv'                    ''              string
    %             'seg_file'                    ''              string
    %             'annot_file'                  ''              string
    %             
    %   [Alternative input options for expression data]
    %             'expr_path'                   ''              string
    %             'index_file'                  ''              string
    %             'file_extension'              ''              string
    %             'expr_ratio_csv'              ''              string
    %             
    %   [Additional options]
    %             'window_size'                 10e6            numeric
    %             'neutral_thresh'              0.1             numeric
    %             'min_neutral'                 20              numeric
    %             'pseudo_expr'                 ''              numeric
    %             'pseudo_expr_relative         10              numeric
    %             'max_nan'                     0.1             numeric
    %             'reportdir'                   '.'             string
    %             'normalization'               'percentile'    string
    %                    {'percentile'
    %                    'library_size'
    %                    'none'}
    %             'percentile'                  95              numeric
    %             'optional_gene_annot'         ''              string
    %             'peak_level'                  0.6             numeric
    %                    {0.0-1.0}
    %             'only_focal'                  ''              string
    %             'scorefield'                  'fs_hp'         string

    properties
        % input parameters
        params
        datasource
        output
%        optional
        
        % input data
        annot
        expr
        cna
        seg
        
        % output
        report
        stats
    end
    
    properties (Hidden = true)
        internal
%         n_tumors
%         n_max_nan
%         offset
%         gene_level
    end
    
    methods
        function obj = FocalScan(varargin)
            % Parse input arguments
                [obj.params, obj.datasource, obj.output] = FocalScan.handle_input(varargin{:});

                % Write log file
                mkdir(obj.output.reportdir)
                logfile = [obj.output.reportdir filesep 'log.txt'];
                if exist(logfile,'file') == 2
                    delete(logfile)
                end
                diary(logfile)
            
                % Display parameter values
                disp(obj.datasource)
                disp(obj.params)
                disp(obj.output)

                % Load data
                obj = obj.read_data();

                if any(isnan(obj.expr.data))
                    error('NaN present in expression data (not allowed).')
                end
                
                % Set some additional internal parameters
                sample = unique(obj.cna.sample_id);
                obj.internal.n_tumors = length(sample);
                obj.internal.n_max_nan = obj.internal.n_tumors*obj.params.max_nan;

                % Run main program
                obj = obj.main();
            diary off
        end
        
        function obj = main(obj)
            fprintf('Setting pseudo_expr = %d\n',obj.params.pseudo_expr)
            disp('Calculating scores')
            tic
            obj.stats = obj.focal_stats;
            toc            

            disp('Writing .wig files')
            obj.make_wigs;

            disp('Writing full report')
            obj.report = obj.make_report;
            report_file = [obj.output.reportdir filesep 'report.txt'];
            writetable(obj.report,report_file,'Delimiter','\t','WriteVariableNames',1,'WriteRowNames',0,'FileType','text');

            disp('Detecting peak genes/tiles')
            T = FocalScan.make_table(obj.report,obj.annot,obj.params.scorefield);
            peak_table = FocalScan.make_peak_table(T,obj.params.peak_level,'min_genes',obj.params.min_genes,'do_plot',obj.output.peak_figure,'plot_dir',obj.output.reportdir);
            if ~strcmp(obj.datasource.optional_gene_annot,'') && ~obj.internal.gene_level
                peak_table = FocalScan.get_tile_gene_ids(peak_table,obj.datasource.optional_gene_annot);
            elseif ~strcmp(obj.datasource.optional_gene_annot,'')
                warning('Ingoring optional gene annotation since gene-level analysis was chosen.')
            end

            disp('Saving peak table')
            peak_file = [obj.output.reportdir filesep 'peaks.txt'];
            writetable(peak_table,peak_file,'Delimiter','\t','WriteVariableNames',1,'WriteRowNames',0,'FileType','text');
        end

        function obj = set_offset(obj)
            tile_distance = obj.annot.start(2) - obj.annot.start(1);
            tile_distance = tile_distance*double(~obj.internal.gene_level);
            if tile_distance == 0
                tile_distance = 1;
            end
            obj.internal.offset = ceil(obj.params.window_size/tile_distance/2);
        end
        
        function [stats] = focal_stats(obj)
            len = size(obj.expr.data,1);
            
            if ~obj.params.only_focal
                fs = nan(len,1);
                sum_cna = nan(len,1);
            end
            sum_cna_hp = nan(len,1);
            mean_expr = nan(len,1);
            fs_hp = nan(len,1);
            num_neutral = nan(len,1);
            num_amplified = nan(len,1);
            num_deleted = nan(len,1);
            pearson_corr = nan(len,1);
            pearson_p_val = nan(len,1);

            for i = 1:len %for each gene/tile
                mean_expr(i) = mean(obj.expr.data(i,~isnan(obj.expr.data(i,:))));
                
                if ~obj.params.only_focal
                    [fs(i), sum_cna(i)] = FocalScan.main_score(obj.cna.data(i,:),obj.expr.data(i,:),obj.params,obj.internal.n_tumors,obj.internal.n_max_nan);
                end

                [cna_left, cna_right] = obj.get_cna_data_at_offset(i);
                cna_hp = FocalScan.focal_score(obj.cna.data(i,:), cna_left, cna_right); % transposed
                num_neutral(i) = sum(abs(cna_hp) <= 0.1);
                num_deleted(i) = sum(cna_hp < -0.1);
                num_amplified(i) = sum(cna_hp > 0.1);
                
                if ~strcmp(obj.expr.datasource.expr_ratio_csv,'')
                    [fs_hp(i), sum_cna_hp(i),pearson_corr(i),pearson_p_val(i)] = FocalScan.main_score_ratios(cna_hp,obj.expr.data(i,:),obj.internal.n_tumors,obj.internal.n_max_nan);
                else
                    [fs_hp(i), sum_cna_hp(i),pearson_corr(i),pearson_p_val(i)] = FocalScan.main_score(cna_hp,obj.expr.data(i,:),obj.params,obj.internal.n_tumors,obj.internal.n_max_nan);
                end
            end
            
            if ~obj.params.only_focal
                stats.fs = fs;
                stats.sum_cna = sum_cna;
            end
            stats.fs_hp = fs_hp;
            stats.sum_cna_hp = sum_cna_hp;
            stats.mean_expr = mean_expr;
            stats.num_neutral = num_neutral;
            stats.num_amplified = num_amplified;
            stats.num_deleted = num_deleted;
            stats.pearson_corr = real(pearson_corr);
            stats.pearson_p_val = pearson_p_val;
        end
        
        function [cna_data_left, cna_data_right] = get_cna_data_at_offset(obj,iter)
            current_chr = obj.annot.chr{iter};
            mid_gene = obj.annot.start(iter) + (obj.annot.stop(iter) - obj.annot.start(iter))/2;
            
            if obj.internal.gene_level
                left_pos = mid_gene - obj.internal.offset;
                right_pos = mid_gene + obj.internal.offset;
                same_chr_left = left_pos >= 0;
                same_chr_right = right_pos <= obj.get_end_of_chr(current_chr);
                
                if isempty(obj.seg)
                    if same_chr_left
                        idx_closest = obj.annot.get_closest(current_chr,left_pos);
                        cna_data_left = obj.cna.data(idx_closest,:);
                    else
                        cna_data_left = zeros(size(obj.cna.data(1,:)));
                    end
                    if same_chr_right
                        idx_closest = obj.annot.get_closest(current_chr,right_pos);
                        cna_data_right = obj.cna.data(idx_closest,:);
                    else
                        cna_data_right = zeros(size(obj.cna.data(1,:)));
                    end
                else
                    if same_chr_left
                        cna_data_left = obj.seg.get_seg_means(current_chr,left_pos);
                    else
                        cna_data_left = zeros(size(obj.cna.data(1,:)));
                    end
                    if same_chr_right
                        cna_data_right = obj.seg.get_seg_means(current_chr,right_pos);
                    else
                        cna_data_right = zeros(size(obj.cna.data(1,:)));
                    end
                end
            else
                if iter - obj.internal.offset < 1 %offset == 10000. 10000 tiles = 10 Mbp. If not 10 Mbp exist to left/right then no chromosome to left/right.
                    chr_left = '';
                else
                    chr_left = obj.annot.chr{iter - obj.internal.offset};
                end
                if iter + obj.internal.offset > length(obj.annot.id)
                    chr_right = '';
                else
                    chr_right = obj.annot.chr{iter + obj.internal.offset};
                end
                
                same_chr_left = strcmp(current_chr, chr_left);
                same_chr_right = strcmp(current_chr,chr_right);
                
                if same_chr_left
                    cna_data_left = obj.cna.data(iter - obj.internal.offset,:); % set all = 0 if past chr edge
                else
                    cna_data_left = zeros(size(obj.cna.data(1,:)));
                end
                if same_chr_right
                    cna_data_right = obj.cna.data(iter + obj.internal.offset,:);
                else
                    cna_data_right = zeros(size(obj.cna.data(1,:)));
                end
            end
        end
        
        % Handle output
        function make_wigs(obj)
            if ~exist(obj.output.reportdir,'file') == 7
                mkdir(obj.output.reportdir)
            end

            if ~obj.params.only_focal
                make_wiggle('score',obj.stats.fs,'%d\t%1.2f\n');
            end
            make_wiggle('score_hp',obj.stats.fs_hp,'%d\t%1.2f\n');
            make_wiggle('rna',obj.stats.mean_expr*1000,'%d\t%1.0f\n');
            make_wiggle('sum_cna',obj.stats.sum_cna,'%d\t%1.2f\n');
            make_wiggle('sum_cna_hp',obj.stats.sum_cna_hp,'%d\t%1.2f\n');
            
            function make_wiggle(filename,score,format_string)
                filename = [obj.output.reportdir filesep filename '.wig'];
                
                last_chr = '';
                fid = fopen(filename, 'w');
                
                for i = 1:length(obj.annot.start),
                    if ~strcmp(obj.annot.chr{i}, last_chr)
                        last_chr = obj.annot.chr{i};
                        fprintf(fid, 'variableStep chrom=%s span=500\n', last_chr);
                    end
                    if ~isnan(score(i))
                        fprintf(fid, format_string, obj.annot.start(i), score(i));
                    end
                end
                fclose(fid);
            end
        end
        
        function report = make_report(obj)
            report = struct2table(obj.stats);
            report.gene_id = obj.annot.id;
            report = [report(:,end) report(:,1:end-1)];
        end
        
        function chr_end = get_end_of_chr(obj,chr)
            last_gene_end = max(obj.annot.stop(strcmp(obj.annot.chr,chr))); %end of last gene in chr
            chr_end = double(last_gene_end);
        end
        
        function obj = read_data(obj)
            my_switch = obj.datasource.input_combination;
            if ismember(my_switch, [1 2 3 4 5 6])
                obj.annot = Annot(obj.datasource.annot_file);
                if ~obj.annot.is_gene_level && isempty(obj.datasource.optional_gene_annot)
                   warning('Tile-level analysis chosen, but no gene annotation has been specified. The peak report will not contain information about genes overlapping each tile.');
                end
            end
            if ismember(my_switch, [1 2])
                disp('Reading expression data')
                obj.expr = Expr('annot',obj.annot,'expr_path',obj.datasource.expr_path,...
                    'index_file',obj.datasource.index_file,...
                    'file_extension',obj.datasource.file_extension);
                obj.expr.data = Expr.normalize(obj.expr.data,obj.params.normalization,obj.params.percentile);
            elseif ismember(my_switch, [3 4])
                disp('Reading expression data')
                obj.expr = Expr('expr_csv',obj.datasource.expr_csv);
                obj.expr.data = Expr.normalize(obj.expr.data,obj.params.normalization,obj.params.percentile);
            elseif ismember(my_switch, [5 6])
                disp('Reading expression data')
                obj.expr = Expr('expr_ratio_csv',obj.datasource.expr_ratio_csv);
            end
            if ismember(my_switch, [1 3 5])
                disp('Reading copy number data')
                obj.seg = Seg(obj.datasource.seg_file);
            end
            if ismember(my_switch, [1 3 5])
                disp('Mapping copy number data to gene/tile IDs')
                obj.cna = CNA('seg',obj.seg,'annot',obj.annot);
            elseif ismember(my_switch, [2 4 6])
                disp('Reading copy number data')
                obj.cna = CNA('cna_csv',obj.datasource.cna_csv);
                if length(obj.annot.id) ~= size(obj.cna.data,1)
                    error('fewer columns in the copy number csv file than genes in the annotation')
                end
            end
            
            % remove non-shared samples from expr and cna
            [samples, idx1, idx2] = intersect(obj.cna.sample_id, obj.expr.sample_id);

            if isempty(samples)
               error('No matching samples found in the CNA and expression data') 
            end
            
            if size(obj.expr.data,1) ~= size(obj.cna.data,1)
               error('Expression and copy number data must contain the same number of genes. (If using CSV input for expression data, make sure that the genes in this file matches (and is in the same order as) those in the annotation file.'); 
            end            
            
            obj.cna.sample_id = obj.cna.sample_id(idx1);
            obj.cna.data = obj.cna.data(:,idx1);
            obj.expr.sample_id = obj.expr.sample_id(idx2);
            obj.expr.data = obj.expr.data(:,idx2);
            
            if ~strcmp(obj.datasource.seg_file,'')
                % remove unnecesary samples from seg
                idx3 = ismember(obj.seg.sample_id,obj.cna.sample_id);
                obj.seg.sample_id = obj.seg.sample_id(idx3);
                obj.seg.chr = obj.seg.chr(idx3);
                obj.seg.start = obj.seg.start(idx3);
                obj.seg.stop = obj.seg.stop(idx3);
                obj.seg.seg_mean = obj.seg.seg_mean(idx3);
                obj.seg.num_mark = obj.seg.num_mark(idx3);
            end

            if ~strcmp(obj.datasource.annot_file,'')
                obj.internal.gene_level = obj.annot.is_gene_level;
                obj = obj.set_offset;
            else
                error('Annotation file not specified')
            end
            
            % all samples in seg must be in expr and vice versa
            if isequal(unique(obj.expr.sample_id),unique(obj.cna.sample_id))
                if ~strcmp(obj.datasource.seg_file,'') 
                    if isequal(unique(obj.expr.sample_id),unique(obj.cna.sample_id),unique(obj.seg.sample_id))
                        [~,~,idx] = unique(obj.seg.sample_id);
                        obj.seg.sample_id = idx;
                        if ~obj.internal.gene_level
                            clear obj.seg; % seg file no longer needed if tile-level
                        end
                    else
                        error('Expression data, CNA data and seg data must all have identical and unique sample ids');
                    end
                end
                
                % replace sample ids with unique numeric ids (to avoid
                % computationally expensive string comparisons later on)
                [~,~,idx] = unique(obj.expr.sample_id);
                obj.cna.sample_id = idx;
                obj.expr.sample_id = idx;
            else
                error('Expression data and CNA data must all have identical and unique sample ids');
            end
            
            if obj.params.pseudo_expr == 0
                obj.params.pseudo_expr = median(obj.expr.data(obj.expr.data>0))*obj.params.pseudo_expr_relative;
            end
            
            if strcmp(obj.params.scorefield,'fs') && obj.params.only_focal
               warning('Both scorefield=fs and only_focal=1 were specified. Ignoring the latter.') 
               obj.params.only_focal = 0;
            end
        end
    end
    
    methods (Static)
        function [fs, sum_cna,pearson_corr,pearson_p_val] = main_score(cna_data,expr_data,params,n_tumors,n_max_nan)
            idx_nan = isnan(cna_data) | isnan(expr_data);
            n_nan = sum(idx_nan);
            
            fs = NaN;
            sum_cna = NaN;
            pearson_corr = NaN;
            pearson_p_val = NaN;
            
            if n_nan <= n_max_nan % skip if to many NaN's
                cna_data_real = cna_data(~idx_nan);
                expr_data_real = expr_data(~idx_nan) + params.pseudo_expr;
                
                norm_factor = n_tumors/(n_tumors - n_nan); % normalize for unobserved;
                sum_cna = sum(cna_data_real)*norm_factor;                
                
                idx_neutral = abs(cna_data_real) < params.neutral_thresh;
                if sum(idx_neutral) >= params.min_neutral % estimate of expression level in neutral/diploid samples should not be based on to few samples
                    median_expr_neut = median(expr_data_real(idx_neutral));
                    cna_norm = cna_data_real*norm_factor;
                    rna_norm = log2(expr_data_real/median_expr_neut);
                    fs = rna_norm*cna_norm';
                    if nargout > 3
                        [tmp_c,tmp_p] = corrcoef(rna_norm, cna_norm,'rows','complete');
                        pearson_corr = tmp_c(1,2);
                        pearson_p_val = tmp_p(1,2);
                    end
                end
            end
        end
        
        function [fs, sum_cna,pearson_corr,pearson_p_val] = main_score_ratios(cna_data,rna_norm,n_tumors,n_max_nan)
            idx_nan = isnan(cna_data) | isnan(rna_norm);
            n_nan = sum(idx_nan);
            
            fs = NaN;
            sum_cna = NaN;
            pearson_corr = NaN;
            pearson_p_val = NaN;
            
            if n_nan <= n_max_nan % skip if to many NaN's
                cna_data_real = cna_data(~idx_nan);
                norm_factor = n_tumors/(n_tumors - n_nan);
                cna_norm = cna_data_real*norm_factor;
                rna_norm = rna_norm(~idx_nan);
                fs = rna_norm'*cna_norm;
                
                sum_cna = sum(cna_data_real)*norm_factor;
                if nargout > 3
                     [tmp_c,tmp_p] = corrcoef(rna_norm, cna_norm,'rows','complete');
                     pearson_corr = tmp_c(1,2);
                     pearson_p_val = tmp_p(1,2);
                end
            end
        end
        
        function fs = focal_score(A, L, R)
            % "High pass filter": For each genomic position and sample,
            % the copy number amplitude with minimal absolute value at two
            % positions, a fix distance apart, is subtracted from the copy
            % number amplitude at the current position.
            X = A - max(L, R);
            Y = A - min(L, R);
            
            absX_le_absY = double(abs(X) <= abs(Y));
            
            fs = X.*absX_le_absY + Y.*~absX_le_absY;
        end

        function t = make_peak_table(T,level,varargin)
            p = inputParser;
            p.addParameter('writedir', '', @isstr);
            p.addParameter('do_plot', 0, @isnumeric);
            p.addParameter('min_genes', 100, @isnumeric);
            p.addParameter('plot_dir', '', @isstr);
            p.parse(varargin{:});
            writedir = p.Results.writedir;
            min_genes = p.Results.min_genes;
            do_plot = p.Results.do_plot;
            plot_dir = p.Results.plot_dir;
            if level <= 0
                warning('peak_level less or equal to 0, setting to 0.1')
                level = 0.1;
            elseif level >= 1
                level = 0.99;
            end

            success = 0;
            while ~success
                try
                    % Perform peak detection separately for amplifications and
                    % deletions
                    T_amp = T;
                    T_amp.Score(T.Sum_CNA_HP<0) = 0;
                    ids_amp = FocalScan.findpeaks_wrapper(T_amp,level,'min_genes',min_genes,'do_plot',do_plot,'plot_dir',plot_dir,'peak_file_name','peaks_amp');
                    clear T_amp
                    T_del = T;
                    T_del.Score(T.Sum_CNA_HP>0) = 0;
                    ids_del = FocalScan.findpeaks_wrapper(T_del,level,'min_genes',min_genes,'do_plot',do_plot,'plot_dir',plot_dir,'peak_file_name','peaks_del');
                    clear T_del

                    ids = [ids_amp,ids_del];

                    t = T(ismember(T.Id,ids),:);
                    t = table(t.Id,t.Score,t.Sum_CNA_HP,t.Chr,t.Start,t.Stop);

                    t.Properties.VariableNames = {'Id','Score','Sum_CNA_HP','Chr','Start','Stop'};
                    t = sortrows(t,'Score','descend');
                    t = t(t.Score>0,:); % retain only peaks with positive scores
                    if ~strcmp(writedir,'')
                        if ~exist(writedir,'file')
                            mkdir(writedir)
                        end
                        out_file = [writedir filesep 'peaks.txt'];
                        writetable(t,out_file,'Delimiter','\t','WriteVariableNames',1,'WriteRowNames',0,'FileType','text')
                    end
                    success = 1;
                catch
                    if level < 1
                        warning('Failed to detect peaks at level %.1f. Trying level %.1f.',level,min(level+0.1,1))
                        level = level + 0.1;
                    else
                       error('Peak detection failed') 
                    end
                end
            end
        end

        function t = get_tile_gene_ids(peak_table_path,optional_gene_annot,varargin)
            p = inputParser;
            p.addParameter('writedir', '', @isstr);
            p.parse(varargin{:});
            writedir = p.Results.writedir;

            if istable(peak_table_path)
                t = peak_table_path;
            else
                t = readtable(peak_table_path,'Delimiter','\t','ReadVariableNames',1,'ReadVariableNames',0);
            end
            gene_annot = Annot(optional_gene_annot); % cannot have headers

            % Find gene(s) overlapping each tile
            peak_genes = cell(length(t.Chr),1);
            for i = 1:length(t.Chr) 
               idx = strcmp(gene_annot.chr,t.Chr(i)) & ...
                   gene_annot.start < t.Stop(i) & ...
                   gene_annot.stop > t.Start(i);
               peak_genes{i} = gene_annot.id(idx);
            end
            peak_genes_joined = cellfun(@(x) strjoin(x), peak_genes,'UniformOutput',0);
            t=table(peak_genes_joined,t.Score,t.Sum_CNA_HP,t.Chr,t.Start,t.Stop);
            t.Properties.VariableNames = {'Id','Score','Sum_CNA_HP','Chr','Start','Stop'};
            t = sortrows(t,'Score','descend');
            if ~strcmp(writedir,'')
                if ~exist(writedir,'file')
                    mkdir(writedir)
                end
                out_file = [writedir filesep 'peaks_gene_ids.txt'];
                writetable(t,out_file,'Delimiter','\t','WriteVariableNames',1,'WriteRowNames',0,'FileType','text')
            end
        end

        function T = make_table(report,annot,varargin)
            t1 = table(annot.id,annot.chr,annot.start,annot.stop);
            t1.Properties.VariableNames = {'Id','Chr','Start','Stop'};
            
            if nargin == 3
                scorefield = varargin{1};
            else
                scorefield = 'fs_hp';
            end
                
            t2 = table(report.gene_id,report.(scorefield),report.sum_cna_hp);
            t2.Properties.VariableNames = {'Id','Score','Sum_CNA_HP'};
            if strcmp(scorefield,'sum_cna_hp') || strcmp(scorefield,'sum_cna')
                t2.Score = abs(t2.Score);
            end
            
            % handle rows without IDs in the report file (if annotated with gene Ids)
            if all(isnumeric(t2.Id)) && all(isnumeric(t1.Id))
                if isempty(intersect(t2.Id,t1.Id))
                    error('None of the gene annotation IDs match those in the report file.')
                end
            elseif all(isnumeric(t2.Id)) && ~all(isnumeric(t1.Id))
                error('The IDs in the annotation file do not match those in the report file.')
            elseif all(iscellstr(t2.Id)) && all(iscellstr(t1.Id))
                if isempty(intersect(t2.Id,t1.Id))
                    error('None of the gene annotation IDs match those in the report file.')
                end
                if length(unique(t1.Id)) < length(t1.Id)
                    warning('Duplicate entries detected in annotation. Ignoring these.');
                    [id,~,numeric_id] = unique(t1.Id);
                    counts=FocalScan.count_occurrences(numeric_id);
                    counts.Id = id;
                    idx1 = ismember(t1.Id,counts.Id(counts.Count>1));
                    idx2 = ismember(t2.Id,counts.Id(counts.Count>1));
                    t1 = t1(~idx1,:);
                    t2 = t2(~idx2,:);
                end
                iidx = cellfun(@(x) isempty(x), t2.Id);
                t2 = t2(~iidx,:);
                idx = ismember(t2.Id,t1.Id);
                t2 = t2(idx,:);
                idx = ismember(t1.Id,t2.Id);
                t1 = t1(idx,:);
            end
            
            T = join(t2,t1);
            T = table(T.Id,T.Score,T.Sum_CNA_HP,T.Chr,T.Start,T.Stop);
            T.Properties.VariableNames = {'Id','Score','Sum_CNA_HP','Chr','Start','Stop'};
            T.Score(isnan(T.Score))=0;
            T = sortrows(T,[4 5]);
        end
        
        function counts = count_occurrences(x)
            a = unique(x(:,1));
            y = zeros(size(a));
            for i = 1:length(a)
                if iscell(a(1))
                    y(i) = sum(strcmp(x(:,1),a(i)));
                else
                    y(i) = sum(x(:,1)==a(i));
                end
            end
            counts = table(a,y);
            counts.Properties.VariableNames = {'Variable','Count'};
        end
        
        function peak_ids_merged = findpeaks_wrapper(T,level,varargin)
            % FINDPEAKS_WRAPPER Used to select parameters for findpeaks_fs and to 
            %   control plotting.

            p = inputParser;
            p.addParameter('do_plot', 0, @isnumeric);
            p.addParameter('min_genes', 100, @isnumeric);
            p.addParameter('plot_dir', '', @isstr);
            p.addParameter('peak_file_name','peaks',@isstr);

            p.parse(varargin{:});
            min_genes = p.Results.min_genes;
            do_plot = p.Results.do_plot;
            plot_dir = p.Results.plot_dir;
            peak_file_name = p.Results.peak_file_name;

            if isempty(level)
                % find the first ("largest scale") level that yields at least min_genes peaks
                max_level = get_max_level(T);
                p = [];
                flag = length(p) <= min_genes;
                i = max_level;
                while flag
                    for i = max_level:-1:1
                        p = findpeaks_wrapper(T,i);
                        if length(p) > min_genes
                            flag = 0;
                            break
                        elseif i == 1 && length(p) <= min_genes
                            flag = 0;
                            break
                        else
                            flag = length(p) <= min_genes;
                        end
                    end
                end
                level = i;
                fprintf('Selecting level %d',level)
            elseif level>0 && level<1 % if level is between 0 and 1, use as percentage of max level
                max_level = get_max_level(T);
                level = floor(max_level*level);
            elseif level == 0
                peak_ids_merged = [];
                return
            elseif level<0
                error('level cannot be less than 0')
            end

            unique_chrs = unique(T.Chr);

            pks = cell(length(unique_chrs),1);
            locs = cell(length(unique_chrs),1);
            peak_ids = cell(length(unique_chrs),1);
            xs_all = cell(length(unique_chrs),1);
            sig_all = cell(length(unique_chrs),1);
            for i = 1:length(unique_chrs)
                idx = strcmp(T.Chr,unique_chrs(i));
                [pks{i},locs{i},peak_ids{i},~,xs_all{i},sig_all{i}] = FocalScan.findpeaks_fs(T.Score(idx),T.Id(idx),level);
            end

            % Collect peak ids from all chromosomes
            for i = 1:length(unique_chrs)
                if i == 1
                    peak_ids_merged = peak_ids{i}';
                else
                    peak_ids_merged = [peak_ids_merged peak_ids{i}'];
                end
            end

            % plot peaks from each level and chromosome
            if do_plot
                if ~strcmp(plot_dir,'')
                    if ~exist(plot_dir,'file') == 7
                        mkdir(plot_dir)
                    end
                    figure('Visible','off');
                else
                    figure;
                end
                hold on
                clear x1
                for i = 1:length(unique_chrs)
                    a=find(strcmp(unique_chrs(i),T.Chr));
                    for j = length(xs_all{i})
                        locs_a = a(xs_all{i}{j});
                        plot(locs_a,sig_all{i}{j})
                    end
                    plot(a(locs{i}),pks{i},'v')
                    x1(i) = min(a);
                end
                for i = 1:length(unique_chrs)
                    y1=get(gca,'ylim');
                    plot([x1(i) x1(i)],y1)
                    chr_name = char(unique_chrs(i));
                    chr_name = chr_name(4:end); % assumes that chromosome names start with "chr"
                    text(x1(i),max(y1),chr_name,'Interpreter','none')
                end
                set(gca,'XTick',[])
                axis tight
                if ~strcmp(plot_dir,'')
                    mkdir(plot_dir)
                    fname = [plot_dir filesep peak_file_name '.pdf'];
                    set(gcf, 'PaperUnits', 'inches');
                    set(gcf, 'PaperSize', [30 25]);
                    set(gcf, 'PaperPositionMode','auto')
                    set(gcf, 'renderer', 'painters');
                    print(gcf,fname,'-dpdf')
                end
            end

            function max_level = get_max_level(T)
                jj = 1;
                while 1
                    pp = FocalScan.findpeaks_wrapper(T,jj);
                    if jj == 1
                        p_old = pp;
                    elseif isequal(p_old,pp)
                        break
                    else
                        p_old = pp;
                    end
                    jj = jj + 1;
                end
                max_level = jj;
            end
        end
        
        function [pks_out,locs_out,peak_ids,num_peaks,xs_all,sig_all] = findpeaks_fs(signal,ids,level)
            sig = signal;
            xs = 1:length(signal);

            for i = 1:level
                if i == 1
                    pks_old = sig;
                    locs_old = xs;
                end
                xs_all{i} = xs;
                sig_all{i} = sig;

                [pks,locs] = findpeaks_simple(sig,xs);

                if i == level
                    pks_out = pks;
                    locs_out = locs;
                    peak_ids = ids(locs);
                end
                num_peaks(i) = length(pks);

                if i ~= level
                    if isempty(pks)
                        %warning('Maximum possible level reached: %d',i-1);
                            pks_out = pks_old;
                            locs_out = locs_old;
                            peak_ids = ids(locs_old);
                    elseif length(pks) == 1
                        %warning('Maximum possible level reached: %d',i);
                            pks_out = pks;
                            locs_out = locs;
                            peak_ids = ids(locs);
                        break
                    else
                        vq = pks;
                        xq = locs;

                        sig = vq;
                        xs = xq;

                        pks_old = pks;
                        locs_old = locs;
                    end
                end
            end

            function [pks,locs] = findpeaks_simple(v,l)
                % a peak is defined if it is the largest of its two neighbors
                % v: signal vector
                % l: position vector
                n = 1;
                for ii = 1:length(v)
                    if ii == 1
                        if length(v) == 1
                            pks = v;
                            locs = l;
                            break;
                        end
                        if v(ii) > v(ii+1)
                            pks(n) = v(ii);
                            locs(n) = l(ii);
                            n = n + 1;
                        end
                    elseif ii == length(v) 
                        if v(ii) > v(ii-1)
                            pks(n) = v(ii);
                            locs(n) = l(ii);
                            n = n + 1;
                        end
                    elseif (v(ii) > v(ii-1) && v(ii) > v(ii+1))
                        pks(n) = v(ii);
                        locs(n) = l(ii);
                        n = n + 1;
                    end
                end
                if n == 1
                    [pks,locs_idx] = max(v);
                    locs = l(locs_idx);
                end
            end
        end

        function [params, datasource, output] = handle_input(varargin)
            p = inputParser;

  %          p.addParameter('install_dir','./',@isstr);
  %          p.addParameter('current_dir','./',@isstr);

            % optional csv-formatted data input?
            p.addParameter('expr_csv', '',@isstr);
            p.addParameter('expr_ratio_csv', '',@isstr);
            p.addParameter('cna_csv', '',@isstr);
            
            % optional unformatted data input
            p.addParameter('expr_path', '', @isstr);
            p.addParameter('index_file', '', @isstr);
            p.addParameter('file_extension','',@isstr);
            p.addParameter('seg_file', '', @isstr);
            p.addParameter('annot_file', '', @isstr);
            
            % optional parameters
            p.addParameter('window_size', 10e6, @isnumeric);
            p.addParameter('neutral_thresh', 0.1, @isnumeric);
            p.addParameter('min_neutral', 20, @isnumeric);
            p.addParameter('pseudo_expr', 0, @isnumeric);
            p.addParameter('pseudo_expr_relative', 10, @isnumeric);
            p.addParameter('max_nan', 0.1, @isnumeric);
            p.addParameter('reportdir', '.', @isstr);
            p.addParameter('normalization', 'percentile', @isstr);
            p.addParameter('percentile', 95, @isnumeric);
            p.addParameter('optional_gene_annot', '', @isstr);
            p.addParameter('peak_level', 0.6, @isnumeric);
            p.addParameter('min_genes', 0, @isnumeric);
            p.addParameter('only_focal',0, @isnumeric);
            p.addParameter('peak_figure',0, @isnumeric);
            p.addParameter('scorefield','fs_hp', @isstr);
            
            % parse parameters
            p.parse(varargin{:});
        
 %           install_dir = p.Results.install_dir;
 %           current_dir = p.Results.current_dir;
 %           cd(current_dir)

%            addpath(genpath(install_dir))

            datasource.expr_path = p.Results.expr_path;
            datasource.index_file = p.Results.index_file;
            datasource.file_extension = p.Results.file_extension;
            datasource.seg_file = p.Results.seg_file;
            datasource.annot_file = p.Results.annot_file;
            
            datasource.expr_csv = p.Results.expr_csv;
            datasource.expr_ratio_csv = p.Results.expr_ratio_csv;
            
            datasource.cna_csv = p.Results.cna_csv;
            datasource.optional_gene_annot = p.Results.optional_gene_annot;
            
            params.window_size = p.Results.window_size;
            params.neutral_thresh = p.Results.neutral_thresh;
            params.min_neutral = p.Results.min_neutral;
            params.pseudo_expr = p.Results.pseudo_expr;
            params.pseudo_expr_relative = p.Results.pseudo_expr_relative;
            params.max_nan = p.Results.max_nan;
            params.normalization = p.Results.normalization;
            params.percentile = p.Results.percentile;
            params.peak_level = p.Results.peak_level;
            params.min_genes = p.Results.min_genes;
            params.only_focal = p.Results.only_focal;
            params.scorefield = p.Results.scorefield;

            output.peak_figure = p.Results.peak_figure;
            output.reportdir = p.Results.reportdir;

            datasource.input_combination = select_input_combination();
            
            function input_combination = select_input_combination()
                expression_set_1 = {...
                    datasource.expr_path,...
                    datasource.index_file,...
                    datasource.file_extension,...
                    datasource.annot_file};
                expression_set_2 = {datasource.expr_csv};
                
                expression_set_3 = {datasource.expr_ratio_csv};
                
                cna_set_1 = {...
                    datasource.seg_file,...
                    datasource.annot_file};
                cna_set_2 = {...
                    datasource.cna_csv,...
                    datasource.annot_file};
                
                combination{1} = [expression_set_1 cna_set_1];
                combination{2} = [expression_set_1 cna_set_2];
                combination{3} = [expression_set_2 cna_set_1];
                combination{4} = [expression_set_2 cna_set_2];
                combination{5} = [expression_set_3 cna_set_1];
                combination{6} = [expression_set_3 cna_set_2];
                
                n = 1;
                input_combination(n) = 0;
                for i = 1:length(combination)
                    if all(cellfun(@(x) ~isempty(x), combination{i}))
                        input_combination(n) = i;
                        n = n + 1;
                    end
                end
                input_combination = max(input_combination);
                
                if ~input_combination
                    disp(datasource)
                    error('Incomplete or redundant specification of input data sources')
                end
            end
        end
    end
end
