classdef Annot < handle
    properties
        chr
        start
        stop
        id
    end

    methods
        function obj = Annot(varargin)
            if nargin >= 1
                in = varargin{1};
                if isstruct(in) || isa(in,'Annot')
                    obj.chr = in.chr;
                    obj.start = in.start;
                    obj.stop = in.stop;
                    obj.id = in.id;
                elseif istable(in)
                    obj.chr = in.Var1';
                    obj.start = in.Var2';
                    obj.stop = in.Var3';
                    obj.id = in.Var4';
                elseif iscell(in)
                    obj.chr = in{1}';
                    obj.start = in{2}';
                    obj.stop = in{3}';
                    obj.id = in{4}';
%                 elseif regexp(in,'\.gtf')
%                     fid = fopen(in,'r');
%                     line = fgetl(fid);
%                     fclose(fid);
%                     isStrCol = isnan( str2double( regexp( line, '[^\t\s]+', 'match' )));
%                     fmt = cell( 1, numel(isStrCol) );
%                     fmt(isStrCol)  = {'%s'};
%                     fmt(~isStrCol) = {'%d'};
% 
%                     split_line=regexp( line, '[^\t\s]+', 'match' );
%                     idx = strcmp(split_line,'gene_id');
% 
%                     idx_gene_id = find(idx) + 1;
% 
%                     fmt_2 = ['%s', '%*s', '%*s', '%d', '%d', repmat({'%*s'},1,idx_gene_id-6),fmt(idx_gene_id),repmat({'%*s'},1,length(fmt)-idx_gene_id)];
%                     %fmt_2 = ['%s',repmat({'%*s'},1,2),fmt([4,5]),repmat({'%*s'},1,idx_gene_id-6),fmt(idx_gene_id),repmat({'%*s'},1,length(fmt)-idx_gene_id)];
%                     fmt_2 = [fmt_2{:}];
% 
%                     fid = fopen(in,'r');
%                     gtf = textscan(fid,fmt_2,'Delimiter',{'\t',' '});
%                     fclose(fid);
% 
%                     gene_id = cellfun(@(x) regexp(x, '(?<=")[^"]+(?=")', 'match'), gtf{4}); % remove double quotes and semicolon
%                     
%                     obj.chr = gtf(1);
%                     obj.chr = obj.chr{1};
%                     obj.start = gtf{2};
%                     obj.stop = gtf{3};
%                     obj.id = gene_id;
%                     clear gtf;
% 
%                     %if length(strsplit(s.chr{1},'chr')) == 1
%                     if isempty(regexp(obj.chr{1},'chr*', 'once')) % add "chr" to chromosome number if not present
%                         for i = 1:length(obj.chr)
%                             obj.chr{i} = ['chr' obj.chr{i}];
%                         end
%                     end
% 
%                     temp_id = unique(obj.id);
%                     [unique_id,~,numeric_id] = unique(obj.id);
%                     [~,temp_numeric_id] = ismember(temp_id,unique_id); % where is unique_id found in temp_id: length of temp_id
% 
%                     temp_start = zeros(1,length(temp_numeric_id));
%                     temp_stop = temp_start;
%                     temp_chr = cell(1,length(temp_numeric_id));
%                     
%                     isint = @(x)(~isempty(x)&&isnumeric(x)&&isreal(x)&&isfinite(x)&&(x==fix(x)));
%                     if nargin == 2 && isint(varargin{2}) && varargin{2} > 1
%                         n_proc = varargin{2};
%                         parpool(n_proc);
%                         
%                         parfor i = 1:length(temp_numeric_id) % for each gene id
%                             idx = numeric_id == temp_numeric_id(i); % where in numeric_ids is this gene id
%                             t_chr = obj.chr(idx); % what is the chromosome of this gene_id
%                             temp_chr(i) = t_chr(1);
%                             temp_start(i) = min(obj.start(idx));
%                             temp_stop(i) = max(obj.stop(idx));
%                         end
%                         
%                         delete(gcp('nocreate'))
%                     else
%                         for i = 1:length(temp_numeric_id) % for each gene id
%                             idx = numeric_id == temp_numeric_id(i); % where in numeric_ids is this gene id
%                             t_chr = obj.chr(idx); % what is the chromosome of this gene_id
%                             temp_chr(i) = t_chr(1);
%                             temp_start(i) = min(obj.start(idx));
%                             temp_stop(i) = max(obj.stop(idx));
%                         end
%                     end
% 
%                     obj.chr = temp_chr;
%                     obj.start = temp_start;
%                     obj.stop = temp_stop;
%                     obj.id = temp_id';
                elseif regexp(in,'\.bed')
                    if exist(in, 'file') ~= 2
                        error('File %s not found',in)
                    end
                    fid = fopen(in,'r');
                    if fid < 0
                        error('Could not open file %s',in)
                    end
                        
                    filehandle = fopen(in);
                    s = textscan(filehandle,'%s %d %d %d%*[^\n]');
                    obj.chr = s{1};
                    obj.start = s{2};
                    obj.stop = s{3};
                    obj.id = s{4};
                    
                    if isempty(regexp(obj.chr{1},'chr*', 'once'))
                        obj.chr = strcat('chr',obj.chr);
                    end
                else
                    disp('Error: Annotation file needs to have the extension .bed or .gtf or be a struct.')
                end
            end
        end

        function gene_level = is_gene_level(obj)
            tile_distances = diff(obj.start);
            is_chr_start = tile_distances ~= mode(tile_distances);
            if (all(tile_distances(~is_chr_start) == mode(tile_distances))) && (sum(is_chr_start) == length(unique(obj.chr)) - 1) % all tile distances are the same except those of the first position on each chromosome
                gene_level = 0;
            else
                gene_level = 1;
            end
        end
        
        function idx = get_closest(obj,chr,pos)
            idx = strcmp(obj.chr,chr);
            [before_start,gene_before] = max(obj.start(idx & obj.start <= pos));
            [after_start,gene_after] = min(obj.start(idx & obj.start >= pos));
            
            if isempty(before_start) && isempty(after_start)
                error('%s:%d: Found no other annotated position on this chromosome',chr,pos)
            elseif isempty(before_start)
                idx = gene_after;
            elseif isempty(after_start)
                idx = gene_before;
            elseif (pos-before_start) <= (after_start - pos)
                idx = gene_before;
            else
                idx = gene_after;
            end
        end

%         function idx = get_index_of_gene(obj,id)
%             idx = find(strcmp(obj.id,id));
%         end

%         function annot = reduce(obj,idx_or_chr)
%             if ischar(idx_or_chr)
%                 idx = strcmp(obj.chr,idx_or_chr);
%             else
%                 idx = idx_or_chr;
%             end
%             a.chr = obj.chr(idx);
%             a.start = obj.start(idx);
%             a.stop = obj.stop(idx);
%             a.id = obj.id(idx);
% 
%             annot = Annot(a);
%         end
    end
end
