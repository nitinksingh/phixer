function thresholded_nw = phixer_threshold(pruned_file, threshold)
    if nargin < 1
        error('Missing input file\n');
    end
    % Default threshold -- keep edges whoes weight is in top 1 percentile.
    if nargin < 2
        threshold = 0.99;
    end
    
    % Load input file
    thresholded_nw =  dlmread(pruned_file, '\t');
    
    % Convert to matrix
    thresholded_nw = spconvert(thresholded_nw);
    
    % Drop threshold percentile edges and compute total degree
    x = nonzeros(thresholded_nw); 
    t = quantile(x, threshold); 
    thresholded_nw(thresholded_nw<t) = 0;
    
    % Print number of strongly connected components
    fprintf('Graph SCC: %d\n',  graphconncomp(thresholded_nw));

end