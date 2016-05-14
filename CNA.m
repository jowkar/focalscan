classdef CNA < handle
    properties
        sample_id
        data
        gene_id
    end
    
    properties (SetAccess = protected, Hidden = true)
        datasource
    end
    
    methods
        function obj = CNA(varargin)
            obj = obj.handle_input(varargin{:});
            
            switch obj.datasource.input_combination
                case 1 % segmented copy number input
                    annot = obj.datasource.annot;
                    seg = obj.datasource.seg;
                    
                    obj.sample_id = Expr.make_valid_var_names(unique(seg.sample_id));
                    
                    obj.data = NaN(length(annot.id), length(obj.sample_id), 'single');
                    
                    % divide tile_annot into per-chromosome pieces for faster matching later
                    chrs = unique(seg.chr);
                    tile_idx = cell(length(chrs),1);
                    tile_center = cell(length(chrs),1);
                    for i = 1:length(chrs)
                        tile_idx{i} = find(strcmp(annot.chr, chrs{i}));
                        tile_center{i} = (annot.start(tile_idx{i}) + annot.stop(tile_idx{i}))/2;
                    end
                    
                    % map segmented copy-number data to genomic tiles/genes, chromosome-wise
                    for i = 1:length(obj.sample_id)
                        idx = find(strcmp(seg.sample_id, obj.sample_id{i}));
                        for j = 1:length(idx)
                            % now check what tiles/genes are covered by segment j in *.seg
                            idx_chr = strcmp(chrs, seg.chr(idx(j)));
                            % find tiles/genes that are covered or partially covered by this segment
                            % (based on tile center position, so that tiles will be assigned to
                            % whatever segment is covering the most)
                            idx_covered = tile_center{idx_chr} >= seg.start(idx(j)) & tile_center{idx_chr} <= seg.stop(idx(j));
                            obj.data(tile_idx{idx_chr}(idx_covered), i) = single(seg.seg_mean(idx(j)));
                        end
                    end
                    
                    obj.gene_id = annot.id;
                case 2 % csv input
                    [obj.data,obj.sample_id,obj.gene_id] = Expr.readcsv(obj.datasource.cna_csv);
                otherwise
                    error('Invalid input')
            end
            
            if size(obj.data,1) < size(obj.data,2)
                warning('More samples than genes provided. Check that samples correspond to columns and genes to rows in CNA data input.')
            end
        end
        
        function obj = handle_input(obj,varargin)
            p = inputParser;
            
            % "raw" data input
            p.addParameter('seg','',@(x) isa(x,'Seg'));
            p.addParameter('annot','',@(x) isa(x,'Annot'));
            
            % csv-formatted data input
            p.addParameter('cna_csv', '',@isstr);
            
            p.parse(varargin{:});
           
            obj.datasource.annot = p.Results.annot;
            obj.datasource.seg = p.Results.seg;
            obj.datasource.cna_csv = p.Results.cna_csv;
            
            obj.datasource.input_combination = select_input_combination;
            
            function input_combination = select_input_combination()
                cna_set_1 = {obj.datasource.seg,obj.datasource.annot};
                cna_set_2 = {obj.datasource.cna_csv};
                
                combination{1} = cna_set_1;
                combination{2} = cna_set_2;

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
                    error('Invalid copy number data input')
                end
            end
        end
    end
end
