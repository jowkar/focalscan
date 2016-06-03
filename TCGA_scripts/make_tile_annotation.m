function make_tile_annotation(chromInfo)

    %% generate tiled genome files

    tilesize = 1000;
    stepsize = 500;

    [hg19.chr hg19.size hg19.x] = textread(chromInfo, '%s%d%s', 'delimiter', '\t');

    chr = unique(hg19.chr);
    %chr(search(chr,'_')) = [];
    chr(cellfun(@(x) ~isempty(x), regexp(chr,'_'))) = [];
    fid1 = fopen('1kb_tiles.bed', 'w');
    fid2 = fopen('1kb_tiles_nochr.bed', 'w');
    n = 0;
    for i = 1:length(chr),
        chr_this = chr{i}
        chr_this_stripped = strrep(chr_this, 'chr', '');
        maxpos = hg19.size(strcmp(hg19.chr, chr_this));
        for j = 1:stepsize:(maxpos - tilesize),
            n = n + 1;
            fprintf(fid1, '%s\t%d\t%d\t%d\n', ...
                chr_this, j, j + tilesize - 1, n);
            fprintf(fid2, '%s\t%d\t%d\t%d\n', ...
                chr_this_stripped, j, j + tilesize - 1, n);
        end
    end
    fclose(fid1);
    fclose(fid2);
end
