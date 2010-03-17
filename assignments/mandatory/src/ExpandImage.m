function G = ExpandImage(I, n)
% TODO: Integrate the function ExpandSignal.m
% Continues edge pixel to the specified size
G = AddBorder(I, n);

row_copy = I(1,:);
for i = 1:n
    for j = 1:length(row_copy)
        G(i, n + j) = row_copy(j);
    end
end

row_copy = I(size(I,1),:);
for i = size(I,1) + 1 + n:size(I,1) + 2*n
    for j = 1:length(row_copy)
        G(i, n + j) = row_copy(j);
    end
end

col_copy = G(:,n + 1);
for i = 1:n
    for j = 1:length(col_copy)
        G(j, i) = col_copy(j);
    end
end

col_copy = G(:,size(G,2) - n);
for i = size(G,2) - n:size(G,2)
    for j = 1:length(col_copy)
        G(j, i) = col_copy(j);
    end
end

end
