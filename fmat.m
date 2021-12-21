function k = fmat(X0,F, J, n)
%based on f function except that it does calcs for matrix from X where each
%column is a different ensemble member 

k=zeros(J,n);

%first the 3 problematic cases: 1,2,J
k(1,:)=(X0(2,:)-X0(J-1,:)).*X0(J,:)-X0(1,:);
k(2,:)=(X0(3,:)-X0(J,:)).*X0(1,:)-X0(2,:);
k(J,:)=(X0(1,:)-X0(J-2,:)).*X0(J-1,:)-X0(J,:);

%then the general case
for j=3:J-1
 k(j,:)=(X0(j+1,:)-X0(j-2,:)).*X0(j-1,:)-X0(j,:);
end
%add the F    
k=k+F;
