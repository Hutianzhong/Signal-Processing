function tek=TEK(C)
%Es = 100;
N=length(C);
CC(1)=0;
CC(N)=0;

for j=2:N-1
    CC(j)=C(j)^2-C(j-1)*C(j+1);
end
tek=kurtosis(CC);
end