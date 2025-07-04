function mult = parMult(A,B)
mult = zeros(length(A),1);
parfor i=1:length(A)
    mult(i) = A(i)*B(i);
end

end

