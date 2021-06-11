function n = EucledianNorm(V)
%Function to calculate eucledian norm
tmp = V;
tmp = tmp.*tmp;
n=sum(tmp);
n=sqrt(n);
end

