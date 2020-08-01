%Generates case corresponding to Panday et al. 2013 http://dx.doi.org/10.1016/j.powtec.2014.02.010

%Select case 
selectedCase = 2;
datastruct = importdata('Panday_data.dat',' ',1);
alphaRecirc = 0.2;
areaRec     = 0.0417;

%--- END USER INPUT -------------------------
%--------------------------------------------

data = datastruct.data;



file = fopen('particleDiameter.H','w')
fprintf(file,'d %f;', data(selectedCase,1) );
fclose(file);

file = fopen('particleDensity.H','w')
fprintf(file,'rho %f;', data(selectedCase,2));
fclose(file);

file = fopen('Uriser.H','w')
fprintf(file,'Uriser %f;', data(selectedCase,3));
fclose(file);

file = fopen('UrecircAir.H','w')
fprintf(file,'UrecircAir %f;', -data(selectedCase,4)/areaRec);
fclose(file);

file = fopen('alphaRecirc.H','w')
fprintf(file,'alphaRecirc %f;', alphaRecirc);
fclose(file);

file = fopen('UrecircP.H','w')
fprintf(file,'UrecircP %f;',- data(selectedCase,5)/(areaRec*alphaRecirc*data(selectedCase,2)));
fclose(file);

file = fopen('pressure.H','w')
fprintf(file,'refPressure %f;', data(selectedCase,6));
fclose(file);


