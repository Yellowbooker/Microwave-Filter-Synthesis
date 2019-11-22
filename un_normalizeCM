%-------------------------------------反归一化耦合矩阵
function unnormalizedCM=un_normalizeCM(CM,CF,FBW)
n=length(CM);
for i=1:n
    for j=1:n
    if j==i
        unnormalizedCM(i,j)=CF*(-CM(i,j)*FBW*0.5+sqrt((0.5*CM(i,j)*FBW)^2+1));
    else
        unnormalizedCM(i,j)=CM(i,j)*FBW*CF;
    end
    end
end
 unnormalizedCM(1,1)=0;
 unnormalizedCM(n,n)=0;
end
