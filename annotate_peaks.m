function peaks = annotate_peaks(peak_file_path,annot_file_path,out_file)

    % ANNOTATE_PEAKS: Add gene IDs to peak report file created by FocalScan. Useful if a tile-level analysis was performed without providing an optional gene annotation.
    %
    % USAGE:
    %     annotate_peaks.sh peaks_file_path annot_file_path out_file
    % 
    %         peak_file_path: path to peak report file created by running FocalScan (ie. peaks.txt)
    %         annot_file_path: path to annotation file (ie. annotation/gencode17.bed)
    %         out_file: name of output file (ie. peaks_annotated.txt)
    % 
    % EXAMPLE:
    %     annotate_peaks.sh peaks.txt annotation/gencode17.bed ./peaks_annotated.txt

    p = inputParser;
    p.addRequired('peaks_path',@isstr);
    p.addRequired('gene_annot',@isstr);
    p.addRequired('out_file',@isstr);
    
    try
        p.parse(peak_file_path,annot_file_path,out_file);
    catch ME
%         if ideployed
%             type help_annotate_peaks.txt
%         else
%             help annotate_peaks.m
%         end
        disp('Invalid input. Run the program without arguments to display a help page (or use the help command within the MATLAB environment).')
        rethrow(ME)
    end
    
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