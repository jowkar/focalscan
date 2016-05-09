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
    
    properties
        sample_id
        data
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
                    
                    expr.sample_id = unique(s.sample_id);
                    expr.data = zeros(length(annot.id), length(expr.sample_id), 'single');
                    
                    n = 0;
                    for i = 1:length(s.name),
                        filename = [expr_path filesep s.name{i} file_extension];
                        if exist(filename, 'file') == 2
                            n = n + 1;
                            fid = fopen(filename);
                            if fid < 0
                                error('%s could not be read. Does the file have a valid format?',filename)
                            end
                            if gene_level
                                disp(filename)
                                expr_this = textscan(fid, '%s%d', 'delimiter', '\t'); %%{1} is tile_id and {2} is expression
                            else
                                disp(filename)
                                expr_this = textscan(fid, '%d%d', 'delimiter', '\t'); %%{1} is tile_id and {2} is expression
                            end
                            [isect,idx1,idx2] = intersect(annot.id, expr_this{1});
                            if isempty(isect)
                                error('No gene identifiers matched the annotation. Please check that an appropriate genome annotation file has been provided.')
                            end
                            idx3 = strcmp(expr.sample_id, s.sample_id{i});
                            expr.data(idx1, idx3) = expr.data(idx1, idx3) + single(expr_this{2}(idx2));
                            fclose(fid);
                        else
                            fprintf('Warning: Sample %s not found\n',s.name{i});
                        end
                    end
                    if n==0
                        error('No samples were read. Check that the correct path to expression files was given.')
                    end
                    fprintf('Read %d samples\n',n);
                    obj.data = expr.data;
                    
                    c = char(version);
                    c = str2double(c(1:3));
                    if c < 9.0
                        obj.sample_id = genvarname(expr.sample_id);
                    else
                        obj.sample_id = matlab.lang.makeValidName(expr.sample_id);
                    end
                case 2 % csv formatted input
                    obj.data = readtable(obj.datasource.expr_csv,'Delimiter',',','ReadRowNames',false,'ReadVariableNames',true,'FileType','text');
                    obj.sample_id = obj.data.Properties.VariableNames;
                    obj.data = table2array(obj.data);
                case 3 % csv formatted log2 ratios
                    obj.data = readtable(obj.datasource.expr_ratio_csv,'Delimiter',',','ReadRowNames',false,'ReadVariableNames',true,'FileType','text');
                    obj.sample_id = obj.data.Properties.VariableNames;
                    obj.data = table2array(obj.data);
                otherwise
                    error('Invalid input')
            end
            
            if size(obj.sample_id,1) < size(obj.sample_id,2)
                obj.sample_id = obj.sample_id'; 
            end
            if size(obj.data,1) > size(obj.data,2)
                obj.data = obj.data';
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
            p.addParameter('normalization','none',@isstr);
            
            % csv-formatted data input
            p.addParameter('expr_csv', '',@isstr);
            p.addParameter('expr_ratio_csv', '',@isstr);
            
            p.parse(varargin{:});
            
            obj.datasource.expr_path = p.Results.expr_path;
            obj.datasource.index_file = p.Results.index_file;
            obj.datasource.file_extension = p.Results.file_extension;
            obj.datasource.annot = p.Results.annot;
            
            obj.datasource.expr_csv = p.Results.expr_csv;
            obj.datasource.expr_ratio_csv = p.Results.expr_ratio_csv;
            
            obj.datasource.input_combination = select_input_combination;
            
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
            transpose = 0;
            if size(data,2) > size(data,1)
                fprintf('Assuming %d gene/tile IDs and %d samples\n',size(data,2),size(data,1))
                data = data';
                transpose = 1;
                %error('Either the data matrix is transposed the wrong way, or fewer genes than samples are included')
            else
                fprintf('Assuming %d gene/tile IDs and %d samples\n',size(data,1),size(data,2))
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
                    for i = 1:size(data,2) % for each sample
                        norm_fact = sort(data(:, i), 'descend');
                        norm_fact = median(norm_fact(1:ceil(end/div)));
                        if norm_fact == 0
                            warning('All genes in sample %d have 0 reads',i)
                            data(:, i) = 0;
                        else
                            data(:, i) = data(:, i)/norm_fact;
                        end
                    end
                case 'none'
                    % Do nothing
                case 'library_size'
                    norm_fact = sum(data);
                    norm_fact = repmat(norm_fact,size(data,1),1);
                    data = (data./norm_fact).*10^6;
                otherwise
                    error('Unrecognized normalization option: %s',mode)
            end
            if transpose == 1
                data = data';
            end
        end
    end
end
