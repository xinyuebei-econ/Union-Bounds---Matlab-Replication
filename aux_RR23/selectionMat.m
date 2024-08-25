function m = selectionMat(selection, size, select)
    if nargin < 3
        select = 'columns';
    end
    
    if strcmp(select, 'rows')
        m = zeros(length(selection), size);
        m(1:length(selection), selection) = eye(length(selection));
    elseif strcmp(select, 'columns')
        m = zeros(size, length(selection));
        m(selection, 1:length(selection)) = eye(length(selection));
    else
        error('Invalid selection type. Must be either "rows" or "columns".');
    end
end


