function mean_A=getMeanOfNonNAN(A)
mean_A=zeros(size(A,1),1);
for i=1:size(A,1)
nonNAN_elements_index = find(~isnan(A(i,:)));
nonNAN_elements = A(i,nonNAN_elements_index);
mean_A(i) = mean(nonNAN_elements);
end
end