classdef Expr < handle
    % EXPR Object used to read, store and normalize expression data
    %             Parameter:                    Default:    Type:
    %
    %   [CSV-formatted gene-level data]
    %             'expr_csv'                    ''              string
    %             'expr_ratio_csv'              ''              string
    %             
    %   [Tab delimited read count files]
    %             'expr_path'                   ''              string
    %             'index_file'                  ''              string
    %             'file_extension'              ''              string
    %             'annot'                       ''              Annot
    %
    %   [Additional options]
    %             'fast_read'                   0               numeric
    
    properties
        sample_id
        data
        gene_id
    end
    
    properties (Hidden = true)
        datasource
        params
    end
    
    methods
        function obj = Expr(varargin)
            obj = obj.handle_input(varargin{:});
            
            switch obj.datasource.input_combination
                case 1 % Separate read count files for each sample
                    expr_path = deblank(obj.datasource.expr_path);
                    index_file = obj.datasource.index_file;
                    annot = obj.datasource.annot;
                    file_extension = obj.datasource.file_extension;
                    if isempty(regexp(file_extension,'\.','match'))
                        file_extension = ['.' file_extension];
                    end
                    
                    fid_index_file = fopen(index_file);
                    if fid_index_file < 0
                        error('Could not open index file %s',index_file)
                    end
                    s = textscan(fid_index_file, '%s%s', 'delimiter', '\t');
                    fclose(fid_index_file);
                    s = cell2struct(s, {'name','sample_id'},2);
                    
                    gene_level = annot.is_gene_level;
                    
                    obj.sample_id = unique(s.sample_id);
                    obj.data = zeros(length(annot.id), length(obj.sample_id), 'single');
                    
                    if obj.params.fast_read
                        disp('fast_read selected: assuming that all expresion data files have identical first columns (gene IDs).')
                    end
                    
                    n = 0;
                    for i = 1:length(s.name),
                        filename = [expr_path filesep s.name{i} file_extension];
                        if exist(filename, 'file') == 2
                            n = n + 1;
                            fid = fopen(filename);
                            if fid < 0
                                error('%s could not be read. Does the file have a valid format?',filename)
                            end
                            
                            if n == 1 || ~obj.params.fast_read
                                if gene_level
                                    disp(filename)
                                    expr_this = textscan(fid, '%s%d', 'delimiter', '\t','HeaderLines',0); %%{1} is tile_id and {2} is expression
                                else
                                    disp(filename)
                                    expr_this = textscan(fid, '%d%d', 'delimiter', '\t','HeaderLines',0); %%{1} is tile_id and {2} is expression
                                end
                                
                                gene_ids = expr_this{1};
                                [isect,idx1,idx2] = intersect(annot.id, gene_ids);
                                if isempty(isect)
                                    error('No gene identifiers matched the annotation. Please check that an appropriate genome annotation file has been provided.')
                                end
                                idx3 = strcmp(obj.sample_id, s.sample_id{i});
                                obj.data(idx1, idx3) = obj.data(idx1, idx3) + single(expr_this{2}(idx2)); % If these is more than one sample or gene with the same name, add the values
                            elseif obj.params.fast_read
                                if gene_level
                                    disp(filename)
                                    expr_this = textscan(fid, '%*s%d', 'delimiter', '\t','HeaderLines',0); %%{1} is tile_id and {2} is expression
                                else
                                    disp(filename)
                                    expr_this = textscan(fid, '%*d%d', 'delimiter', '\t','HeaderLines',0); %%{1} is tile_id and {2} is expression
                                end
                                idx3 = strcmp(obj.sample_id, s.sample_id{i});
                                obj.data(idx1, idx3) = obj.data(idx1, idx3) + single(expr_this{1}(idx2));
                            end
                            fclose(fid);
                        else
                            fprintf('Warning: Sample %s not found\n',s.name{i});
                        end
                    end
                    if n==0
                        error('No samples were read. Check that the correct path to expression files was given.')
                    end
                    fprintf('Read %d samples\n',n);
                    obj.sample_id = Expr.make_valid_var_names(obj.sample_id);
                    obj.gene_id = annot.id;
                case 2 % csv formatted input
                    [obj.data,obj.sample_id,obj.gene_id] = Expr.readcsv(obj.datasource.expr_csv);
                case 3 % csv formatted log2 ratios
                    [obj.data,obj.sample_id,obj.gene_id] = Expr.readcsv(obj.datasource.expr_ratio_csv);
                otherwise
                    error('Invalid expression data input')
            end
            
            if size(obj.data,1) < size(obj.data,2)
                warning('More samples than genes provided. Check that samples correspond to columns and genes to rows in expression data input.')
            end
            disp('Finished reading expression data')
        end
        
        function obj = handle_input(obj,varargin)
            p = inputParser;
            
            p.addParameter('annot','', @(x) isa(x,'Annot'));
            
            % raw data input
            p.addParameter('expr_path', '', @isstr);
            p.addParameter('index_file', '', @isstr);
            p.addParameter('file_extension','',@isstr);
%             p.addParameter('normalization','none',@isstr);
            
            % csv-formatted data input
            p.addParameter('expr_csv', '',@isstr);
            p.addParameter('expr_ratio_csv', '',@isstr);
            
            p.addParameter('fast_read',0,@isnumeric);
            
            p.parse(varargin{:});
            
            obj.datasource.expr_path = p.Results.expr_path;
            obj.datasource.index_file = p.Results.index_file;
            obj.datasource.file_extension = p.Results.file_extension;
            obj.datasource.annot = p.Results.annot;
            
            obj.datasource.expr_csv = p.Results.expr_csv;
            obj.datasource.expr_ratio_csv = p.Results.expr_ratio_csv;
            
            obj.datasource.input_combination = select_input_combination;
            
            obj.params.fast_read = p.Results.fast_read;
            %obj.params.normalization = p.Results.normalization;
            
            function input_combination = select_input_combination()
                expression_set_1 = {...
                    obj.datasource.expr_path,...
                    obj.datasource.index_file,...
                    obj.datasource.file_extension,...
                    obj.datasource.annot};
                expression_set_2 = {obj.datasource.expr_csv};
                
                expression_set_3 = {obj.datasource.expr_ratio_csv};
                
                combination{1} = expression_set_1;
                combination{2} = expression_set_2;
                combination{3} = expression_set_3;

                n = 1;
                input_combination(n) = 0;
                for i = 1:length(combination)
                    if all(cellfun(@(x) ~isempty(x), combination{i}))
                        input_combination(n) = i;
                        n = n + 1;
                        %break
                    end
                end
                input_combination = max(input_combination);
                
                if ~input_combination
                    disp(obj.datasource)
                    error('Invalid expression data input')
                end
            end
        end
    end

    methods (Static)
        function data = normalize(data,mode,varargin)
            num_samples = size(data,2);
            num_genes = size(data,1);
            if num_samples > num_genes
                warning('More samples (columns) than rows (genes) in the expression data. Make sure samples correspond to columns in the data input.')
            else
                fprintf('Assuming %d gene/tile IDs and %d samples\n',num_genes,num_samples)
            end
            fprintf('Performing %s normalization\n',mode)
            switch mode
                case 'percentile'
                    if nargin >= 3
                        percentile = varargin{1};
                        div = 100/(100-percentile);
                    else
                        div = 20;
                    end
                    for i = 1:num_samples % for each sample
                        norm_fact = sort(data(:,i), 'descend');
                        norm_fact = median(norm_fact(1:ceil(end/div)));
                        if norm_fact == 0
                            warning('All genes in sample %d have 0 reads',i)
                            data(:, i) = 0;
                        else
                            data(:, i) = data(:, i)/norm_fact;
                        end
                    end
                case 'percentile_new'
                    if nargin >= 3
                        percentile = varargin{1};
                        div = 100/(100-percentile);
                    else
                        div = 20;
                    end
                    for i = 1:num_samples % for each sample
                        d = data(:,i);
                        d = d(d>0);
                        norm_fact = sort(d, 'descend');
                        norm_fact = median(norm_fact(1:ceil(end/div)));
                        if all(d==0)
                            warning('All genes in sample %d have 0 reads',i)
                            data(:,i) = 0;
                        elseif isnan(norm_fact)
                            warning('Could not calculate a normalization factor for sample %d. Setting all genes of this sample to 0.',i)
                            data(:,i) = 0;
                        elseif norm_fact == 0
                            warning('All genes in sample %d have 0 reads',i)
                            data(:, i) = 0;
                        else
                            data(:, i) = data(:, i)/norm_fact;
                        end
                    end
                case 'none'
                    % Do nothing
                case 'library_size'
                    norm_fact = sum(data,1);
                    if length(norm_fact) ~= num_samples
                        error('Size of normalization factor does not correspond to the number of samples.')
                    end
                    norm_fact = repmat(norm_fact,num_genes,1);
                    data = (data./norm_fact)*10^6;
                otherwise
                    error('Unrecognized normalization option: %s',mode)
            end
        end
        
        function valid_names = make_valid_var_names(names)
            c = char(version);
            c = str2double(c(1:3));
            if c < 9.0
                valid_names = genvarname(names);
            else
                valid_names = matlab.lang.makeValidName(names);
            end
        end
        
        function [data,samples,gene_ids] = readcsv(fname)
            % The first line should countain sample names.
            
            % Skip the first line and read the first column
            fid = fopen(fname,'r');
            firstline = fgetl(fid); % textscan will read the remaining lines
            formatstr = '%s%*[^\n]'; % read only the first column
            firstcolumn = textscan(fid,formatstr,'Delimiter',',');
            fclose(fid);
            firstcolumn_double = str2double(firstcolumn{1}(1:end));
            skipfirst = 0;
            if ~isnumeric(firstcolumn_double) | any(isnan(firstcolumn_double))
                disp('The first column of the csv-file contained strings. Assuming that these are gene IDs.');
                gene_ids = [firstcolumn{:}];
                skipfirst = 1;
            else
                gene_ids = [];
                disp('The first column of the csv-file contained only numbers. Assuming that rows (genes) correspond exactly to rows in the annotation file.');
            end
            
            fid = fopen(fname,'r');
            firstline = fgetl(fid); % textscan will read the remaining lines
            samples = textscan(firstline,'%s','Delimiter',',');
            samples = samples{1};
            formatstring = repmat('%f',1,length(samples));
            if skipfirst
               formatstring = ['%*s' formatstring]; 
            end
            data = textscan(fid,formatstring,'Delimiter', ',', 'CollectOutput',true);
            data = single(data{1});
            fclose( fid );
            
            samples = Expr.make_valid_var_names(samples);
        end
    end
end
