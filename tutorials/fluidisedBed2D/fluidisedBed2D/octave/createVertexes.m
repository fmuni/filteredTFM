%Matlab script for creating blockMesh vertexes for the injection
%tutorial for eulerianFilteredTFM

%INPUT
coreH =         16.79;
coreW =        0.3084;
pipeL =            1.;

pipeOneW =       0.23;
pipeOneAxisH =   0.43;

pipeTwoW =      0.203;
pipeTwoAxisH =  15.88;

dpt    =       0.1;

deltax =       0.03;
deltay =       0.1;


%EVALUATE

pipeOneH = pipeOneAxisH - pipeOneW/2.0; 
pipeTwoH = pipeTwoAxisH - pipeTwoW/2.0; 

file      =fopen('vertices.H','w')
fileTwo   =fopen('blocks.H','w')
fileThree =fopen('walls.H','w')
fileFour  =fopen('empty.H','w')
fileFive  =fopen('outlet.H','w')
fileSix   =fopen('inlet.particles.H','w')
fileSeven =fopen('inlet.air.H','w')


%write main block (part 1)
for z=0:1
    
    fprintf(file, ' ( %f %f %f) \n', 0. ,0., dpt*z );
    fprintf(file, ' ( %f %f %f) \n', coreW ,0., dpt*z );
    fprintf(file, ' ( %f %f %f) \n', coreW ,pipeOneH, dpt*z );
    fprintf(file, ' ( %f %f %f) \n', 0. ,pipeOneH, dpt*z );

end

%write block entry
fprintf(fileTwo, ' hex ( %i %i %i %i %i %i %i %i) (%i %i 1) simpleGrading (1 1 1) \n', ...
                 0, 1, 2, 3, 4, 5, 6, 7,...
                floor(coreW/deltax), floor(pipeOneH/deltay)  );
%write wall faces
fprintf(fileThree, '(%i %i %i %i) \n', 0, 4, 7, 3)
fprintf(fileThree, '(%i %i %i %i) \n', 1, 2, 6, 5)

%write empty faces
fprintf(fileFour, '(%i %i %i %i) \n',0, 1, 2,3)
fprintf(fileFour, '(%i %i %i %i) \n', 4,7, 6,5)

%write inlet air
fprintf(fileSeven, '(%i %i %i %i) \n', 0, 1, 5,4)


%write main block (part 2)
for z=0:1
    
  %  fprintf(file, ' ( %f %f %f) \n', 0. ,pipeOneH, dpt*z );
  %  fprintf(file, ' ( %f %f %f) \n', coreW ,pipeOneH, dpt*z );
   
    fprintf(file, ' ( %f %f %f) \n', 0. ,pipeOneH+pipeOneW, dpt*z );
     fprintf(file, ' ( %f %f %f) \n', coreW ,pipeOneH+pipeOneW, dpt*z );

end




fprintf(fileTwo, ' hex ( %i %i %i %i %i %i %i %i) (%i %i 1) simpleGrading (1 1 1) \n', ...
               3, 2, 9, 8, 7, 6, 11, 10, ...
                floor(coreW/deltax), floor(pipeOneW/deltay)  );

%write wall faces
fprintf(fileThree, '(%i %i %i %i) \n', 3, 7, 10,8)

%write empty faces
fprintf(fileFour, '(%i %i %i %i) \n', 3, 2, 9,8)
fprintf(fileFour, '(%i %i %i %i) \n', 7,6, 11,10)



%write main block (part 3)
for z=0:1
    
   % fprintf(file, ' ( %f %f %f) \n', 0. ,pipeOneH+pipeOneW, dpt*z );
   % fprintf(file, ' ( %f %f %f) \n', coreW ,pipeOneH+pipeOneW, dpt*z );
   
    fprintf(file, ' ( %f %f %f) \n', 0. ,pipeTwoH, dpt*z );
     fprintf(file, ' ( %f %f %f) \n', coreW ,pipeTwoH, dpt*z );

end



fprintf(fileTwo, ' hex ( %i %i %i %i %i %i %i %i) (%i %i 1) simpleGrading (1 1 1) \n', ...
               8, 9, 13, 12, 10, 11, 15, 14, ...
                floor(coreW/deltax), floor((pipeTwoH-pipeOneW-pipeOneH)/deltay)  );

 %write wall faces
fprintf(fileThree, '(%i %i %i %i) \n', 8,10,14,12)
fprintf(fileThree, '(%i %i %i %i) \n', 9,11,15,13)

%write empty faces
fprintf(fileFour, '(%i %i %i %i) \n', 8,9,13,12)
fprintf(fileFour, '(%i %i %i %i) \n', 10,11,15,14)
            



%write main block (part 4)
for z=0:1
    
  %  fprintf(file, ' ( %f %f %f) \n', 0., pipeTwoH, dpt*z );
  %  fprintf(file, ' ( %f %f %f) \n', coreW ,pipeTwoH, dpt*z );
    
    fprintf(file, ' ( %f %f %f) \n', 0. ,pipeTwoH+pipeTwoW, dpt*z );
    fprintf(file, ' ( %f %f %f) \n', coreW ,pipeTwoH+pipeTwoW, dpt*z );

end

fprintf(fileTwo, ' hex ( %i %i %i %i %i %i %i %i) (%i %i 1) simpleGrading (1 1 1) \n', ...
               12,13,17,16,14,15,19,18, ...
                floor(coreW/deltax), floor(pipeTwoW/deltay)  );
            
%write wall faces
fprintf(fileThree, '(%i %i %i %i) \n', 12,14,18,16)

%write empty faces
fprintf(fileFour, '(%i %i %i %i) \n', 12,13,17,16)
fprintf(fileFour, '(%i %i %i %i) \n', 14,15,19,18)


%write main block (part 5)
for z=0:1
    
 %   fprintf(file, ' ( %f %f %f) \n', 0., pipeTwoH+pipeTwoW, dpt*z );
 %   fprintf(file, ' ( %f %f %f) \n', coreW ,pipeTwoH+pipeTwoW, dpt*z );
    
    fprintf(file, ' ( %f %f %f) \n', 0. ,coreH, dpt*z );
    fprintf(file, ' ( %f %f %f) \n', coreW ,coreH, dpt*z );

end



fprintf(fileTwo, ' hex ( %i %i %i %i %i %i %i %i) (%i %i 1) simpleGrading (1 1 1) \n', ...
               16,17,21,20,18,19,23,22, ...
                floor(coreW/deltax), floor((coreH - pipeTwoH-pipeTwoW)/deltay)  );

%write wall faces
fprintf(fileThree, '(%i %i %i %i) \n', 16,18,22,20)
fprintf(fileThree, '(%i %i %i %i) \n', 17,19,23,21)
fprintf(fileThree, '(%i %i %i %i) \n', 20,21,23,22)

%write empty faces
fprintf(fileFour, '(%i %i %i %i) \n', 16,17,21,20)
fprintf(fileFour, '(%i %i %i %i) \n', 18,19,23,22)




%write block 1
for z=0:1
    
  %  fprintf(file, ' ( %f %f %f) \n', coreW ,pipeOneH, dpt*z );
    
    fprintf(file, ' ( %f %f %f) \n', coreW+pipeL ,pipeOneH+pipeOneW, dpt*z );
    fprintf(file, ' ( %f %f %f) \n', coreW+pipeL ,pipeOneH, dpt*z );
  %  fprintf(file, ' ( %f %f %f) \n', coreW ,pipeOneH+pipeOneW, dpt*z );

end

node=nodeBase + count;

node(1) = matchOne(2);
node(5) = matchOne(6);
node(8) = matchOne(7);
node(4) = matchOne(3);

fprintf(fileTwo, ' hex ( %i %i %i %i %i %i %i %i) (%i %i 1) simpleGrading (1 1 1) \n', ...
               2,25,24,9,6,27,26,11, ...
                floor(pipeL/deltax), floor(pipeOneW/deltay)  );
            
%write wall faces
fprintf(fileThree, '(%i %i %i %i) \n', 9,24,26,11)
fprintf(fileThree, '(%i %i %i %i) \n', 2,25,27,6)

%write empty faces
fprintf(fileFour, '(%i %i %i %i) \n', 2,25,24,9)
fprintf(fileFour, '(%i %i %i %i) \n', 6,27,26,11)

%write inlet particles
fprintf(fileSix, '(%i %i %i %i) \n', 25,27,26,24)


%write block 2
for z=0:1
    
   % fprintf(file, ' ( %f %f %f) \n', coreW ,pipeTwoH, dpt*z );
    
    fprintf(file, ' ( %f %f %f) \n', coreW+pipeL ,pipeTwoH+pipeTwoW, dpt*z );
    fprintf(file, ' ( %f %f %f) \n', coreW+pipeL ,pipeTwoH, dpt*z );
 %   fprintf(file, ' ( %f %f %f) \n', coreW ,pipeTwoH+pipeTwoW, dpt*z );

end


fprintf(fileTwo, ' hex ( %i %i %i %i %i %i %i %i) (%i %i 1) simpleGrading (1 1 1) \n', ...
               13,29,28,17,15,31,30,19, ...
                floor(pipeL/deltax), floor(pipeTwoW/deltay)  );
%write wall faces
fprintf(fileThree, '(%i %i %i %i) \n', 17,28,30,19)
fprintf(fileThree, '(%i %i %i %i) \n', 13,29,31,15)

%write empty faces
fprintf(fileFour, '(%i %i %i %i) \n', 13,29,28,17)
fprintf(fileFour, '(%i %i %i %i) \n', 15,31,30,19)

%write inlet outlet
fprintf(fileFive, '(%i %i %i %i) \n', 28,29,31,30)
            




