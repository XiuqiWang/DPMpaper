function mean_A=getMeanOfNonZero(A)
mean_A=zeros(size(A,1),1);
for i=1:size(A,1)
nonzero_elements_index = find(A(i,:) ~= 0 & ~isnan(A(i,:)));
nonzero_elements = A(i,nonzero_elements_index);
mean_A(i) = mean(nonzero_elements);
end
end
