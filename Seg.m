classdef Seg < handle
    properties
        sample_id
        chr
        start
        stop
        num_mark
        seg_mean
    end

    methods
        function obj = Seg(varargin)
            if nargin >= 1
                in = varargin{1};
                if isstruct(in) || isa(in,'Seg')
                    obj.sample_id = in.sample_id;
                    obj.chr = in.chr;
                    obj.start = in.start;
                    obj.stop = in.stop; %change?
                    obj.num_mark = in.num_mark;
                    obj.seg_mean = in.seg_mean;
                elseif ischar(varargin{1})
                    fid = fopen(in);
                    s = textscan(fid, '%s%s%d%d%d%f','delimiter','\t','headerlines',1);
                    fclose(fid);
                    s = cell2struct(s, {'sample_id','chr','start','stop','num_mark','seg_mean'},2);

                    if length(strsplit(s.chr{1},'chr')) == 1 % If the chromosome name does not begin with 'chr'
                        for i = 1:length(s.chr)
                            % convert chromosome names if needed
                            s.chr{i} = ['chr' s.chr{i}];
                        end
                    end

                    obj.chr = s.chr;
                    obj.sample_id = matlab.lang.makeValidName(s.sample_id);
                    obj.start = uint32(s.start);
                    obj.stop = uint32(s.stop);
                    obj.num_mark = s.num_mark;
                    obj.seg_mean = s.seg_mean;
                else
                    error('Input needs to be struct, Seg, matrix or string');
                end
                if size(obj.sample_id,1) < size(obj.sample_id,2)
                    obj.sample_id = obj.sample_id';
                end
            end
        end

        function [cna_data] = get_seg_means(obj,chr,pos) % get CNA data from seg file for gene level
            k = (obj.start <= pos) & (obj.stop >= pos) & strcmp(obj.chr,chr);
%             persistent unique_seg_sample_ids % might cause problems
%             if isempty(unique_seg_sample_ids)
                unique_seg_sample_ids = unique(obj.sample_id);
%             end
            cna_data = NaN(length(unique_seg_sample_ids),1);
            
            [~,locb] = ismember(obj.sample_id(k),unique_seg_sample_ids);
            cna_data(locb) = obj.seg_mean(k);
        end
    end
end
