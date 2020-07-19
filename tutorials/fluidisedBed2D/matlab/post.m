%evaluate closures 

%Strain rate
Sr = [2.12434 1.96963 3.4349 4.99752 0.231221 5.00511];

%Settling velocity
Uset = 3.85987;

%Phase field
phase = [0.450033 0.450065 0.999963 0.999969 0.99997 1];

%filterFr
filterFr = 19.1346;

%filter size
fSize = 0.0793701;

rho = 0.579486;

%evaluate Sarkar micro drag
0.00307* (9.81/(Uset*Uset))^(-6/7)
SMicnu = 0.00307* (9.81/(Uset*Uset))^(-6.0/7.0) ...
     *(fSize^(8.0/7.0))*(Sr).*(phase.^1.544)./(0.62-phase);

SmeNufluid = fSize*fSize*Sr.*(0.330+0.218*phase-0.0485*phase.^2);

SmePfluid =  rho* (9.81/(Uset*Uset))^(5/7)*fSize^(19/7)*Sr.*Sr ...
            .* (0.0661 +0.0164*phase -0.194*phase.*phase);