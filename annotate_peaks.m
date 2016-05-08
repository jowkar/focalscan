function peaks = annotate_peaks(peak_file_path,annot_file_path,out_file)
    p = inputParser;
    p.addRequired('peaks_path',@isstr);
    p.addRequired('gene_annot',@isstr);
    p.addRequired('out_file',@isstr);
    p.parse(peak_file_path,annot_file_path,out_file);
    
    peaks = readtable(peak_file_path,'Delimiter','\t','ReadVariableNames',1,'ReadRowNames',0,'FileType','text');
    
    peaks = FocalScan.get_tile_gene_ids(peaks,annot_file_path);
    
    if ~isempty(out_file)
        pth = strsplit(out_file,filesep);
        if length(pth) > 1
            pth = strjoin(pth(1:end-1),filesep);
            mkdir(pth)
        end
        writetable(peaks,out_file,'Delimiter','\t','FileType','text','WriteVariableNames',1,'WriteRowNames',0);
        fprintf('Successfully wrote peak table to %s\n',out_file)
    end
end