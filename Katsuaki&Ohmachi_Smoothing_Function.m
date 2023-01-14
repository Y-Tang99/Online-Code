% Katsuaki K. and T. Ohmachi (1998). Ground-Motion Characteristics Estimated 
% from Spectral Ratio 779 between Horizontal and Vertical Components of Microtremor. 
% Bull. Seismol. Soc. Am. 88, no. 1, 228-241.


function [FAS_sm]=KOfun(FAS,f,b)

% f=f1;
% FAS=FAS(1,:);

N=length(f);
f_shifted = f/(1+1e-4);

f_c(N)=zeros();    % center frequency
z(N,N)=zeros();
w(N,N)=zeros();

for i=1:N
    f_c(i)=f(i);
    if f_c(i)==0
        w(i,:)=1;
    else
        z(i,:)=f_shifted./f_c(i);
        for j=1:N
            if z(i,j)<0.5
                z(i,j)=0.5;
            elseif z(i,j)>2
                z(i,j)=2;
            end
        end
        z(i,:)=b*(log10(f_shifted./f_c(i)));
        w(i,:)=sin(z(i,:))./z(i,:);
        w(i,:)=w(i,:).^4;
        w(i,:)=w(i,:)./sum(w(i,:));
    end
end
w(isnan(w))=0;  % remove NAN in w
FAS_sm=FAS*w;

FAS_sm(N)=FAS_sm(N-1);
end
