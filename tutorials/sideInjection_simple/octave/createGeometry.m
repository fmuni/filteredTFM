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
nodeBase = [0 1 2 3 4 5 6 7];
count = 0;

%write main block (part 1)

fprintf(file, ' ( %f %f %f) \n', 0. ,pipeOneH*y, 0. );
fprintf(file, ' ( %f %f %f) \n', coreW ,pipeOneH*y, 0. );
fprintf(file, ' ( %f %f %f) \n', coreW ,pipeOneH*y, dpt );
fprintf(file, ' ( %f %f %f) \n', 0. ,pipeOneH*y, dpt );

node =nodeBase;
nodeOld =node;
%write block entry
fprintf(fileTwo, ' hex ( %i %i %i %i %i %i %i %i) (%i 1 %i) simpleGrading (1 1 1) \n', ...
                node(1), node(2), node(3), node(4), node(5), node(6), node(7), node(8), ...
                floor(coreW/deltax), floor(pipeOneH/deltay)  );
%write wall faces
fprintf(fileThree, '(%i %i %i %i) \n', node(1), node(5), node(8),node(4))
fprintf(fileThree, '(%i %i %i %i) \n', node(2), node(6), node(7),node(3))

%write empty faces
fprintf(fileFour, '(%i %i %i %i) \n', node(1), node(5), node(6),node(2))
fprintf(fileFour, '(%i %i %i %i) \n', node(4),node(8), node(7),node(3))

%write air inlet
fprintf(fileSeven, '(%i %i %i %i) \n', node(1),node(2), node(3),node(4))

count = count +4;

%write main block (part 2)
fprintf(file, ' ( %f %f %f) \n', 0. ,pipeOneH+pipeOneW, 0 );
fprintf(file, ' ( %f %f %f) \n', coreW ,pipeOneH+pipeOneW, 0 );
fprintf(file, ' ( %f %f %f) \n', coreW ,pipeOneH+pipeOneW, dpt );
fprintf(file, ' ( %f %f %f) \n', 0 ,pipeOneH+pipeOneW, dpt );


node=nodeBase + count;

nodeOld =node;

%remember
matchOne=node;


fprintf(fileTwo, ' hex ( %i %i %i %i %i %i %i %i) (%i 1 %i) simpleGrading (1 1 1) \n', ...
               node(1), node(2), node(3), node(4), node(5), node(6), node(7), node(8), ...
                floor(coreW/deltax), floor(pipeOneW/deltay)  );

%write wall faces
fprintf(fileThree, '(%i %i %i %i) \n', node(1), node(5), node(8),node(4))

%write empty faces
fprintf(fileFour, '(%i %i %i %i) \n', node(1), node(5), node(6),node(2))
fprintf(fileFour, '(%i %i %i %i) \n', node(4),node(8), node(7),node(3))
count = count +4;

 
fprintf(file, ' ( %f %f %f) \n', 0. ,pipeTwoH, 0);
fprintf(file, ' ( %f %f %f) \n', coreW ,pipeTwoH, 0 );
fprintf(file, ' ( %f %f %f) \n', coreW ,pipeTwoH, dpt );
fprintf(file, ' ( %f %f %f) \n', 0 ,pipeTwoH, dpt );


node=nodeBase + count;
nodeOld =node;



fprintf(fileTwo, ' hex ( %i %i %i %i %i %i %i %i) (%i 1 %i) simpleGrading (1 1 1) \n', ...
               node(1), node(2), node(3), node(4), node(5), node(6), node(7), node(8), ...
                floor(coreW/deltax), floor((pipeTwoH-pipeOneW-pipeOneH)/deltay)  );

 %write wall faces
fprintf(fileThree, '(%i %i %i %i) \n', node(1), node(5), node(8),node(4))
fprintf(fileThree, '(%i %i %i %i) \n', node(2), node(6), node(7),node(3))

%write empty faces
fprintf(fileFour, '(%i %i %i %i) \n', node(1), node(5), node(6),node(2))
fprintf(fileFour, '(%i %i %i %i) \n', node(4),node(8), node(7),node(3))



count = count +4;

%write main block (part 4)

 
 fprintf(file, ' ( %f %f %f) \n', 0. ,pipeTwoH+pipeTwoW, 0 );
 fprintf(file, ' ( %f %f %f) \n', coreW ,pipeTwoH+pipeTwoW, 0 );
 fprintf(file, ' ( %f %f %f) \n', coreW  ,pipeTwoH+pipeTwoW, dpt );
 fprintf(file, ' ( %f %f %f) \n', 0,pipeTwoH+pipeTwoW, dpt );

node=nodeBase + count;

nodeOld =node;

%remember
matchTwo=node;


fprintf(fileTwo, ' hex ( %i %i %i %i %i %i %i %i) (%i 1 %i) simpleGrading (1 1 1) \n', ...
               node(1), node(2), node(3), node(4), node(5), node(6), node(7), node(8), ...
                floor(coreW/deltax), floor(pipeTwoW/deltay)  );
            
%write wall faces
fprintf(fileThree, '(%i %i %i %i) \n', node(1), node(5), node(8),node(4))


%write empty faces
fprintf(fileFour, '(%i %i %i %i) \n', node(1), node(5), node(6),node(2))
fprintf(fileFour, '(%i %i %i %i) \n', node(4),node(8), node(7),node(3))

count = count +4;

%write main block (part 5)
  
fprintf(file, ' ( %f %f %f) \n', 0. ,coreH, 0 );
fprintf(file, ' ( %f %f %f) \n', coreW ,coreH, 0 );
fprintf(file, ' ( %f %f %f) \n', coreW  ,coreH, dpt );
fprintf(file, ' ( %f %f %f) \n', 0 ,coreH, dpt );

node=nodeBase + count;

nodeOld =node;


fprintf(fileTwo, ' hex ( %i %i %i %i %i %i %i %i) (%i 1 %i) simpleGrading (1 1 1) \n', ...
               node(1), node(2), node(3), node(4), node(5), node(6), node(7), node(8), ...
                floor(coreW/deltax), floor((coreH - pipeTwoH-pipeTwoW)/deltay)  );

%write wall faces
fprintf(fileThree, '(%i %i %i %i) \n', node(1), node(5), node(8),node(4))
fprintf(fileThree, '(%i %i %i %i) \n', node(2), node(6), node(7),node(3))
fprintf(fileThree, '(%i %i %i %i) \n', node(5), node(6), node(7),node(8))

%write empty faces
fprintf(fileFour, '(%i %i %i %i) \n', node(1), node(5), node(6),node(2))
fprintf(fileFour, '(%i %i %i %i) \n', node(4),node(8), node(7),node(3))


count = count +4;

%write block 1
for z=0:1
    
  %  fprintf(file, ' ( %f %f %f) \n', coreW ,pipeOneH, dpt*z );
    fprintf(file, ' ( %f %f %f) \n', coreW+pipeL ,pipeOneH, dpt*z );
    fprintf(file, ' ( %f %f %f) \n', coreW+pipeL ,pipeOneH+pipeOneW, dpt*z );

  %  fprintf(file, ' ( %f %f %f) \n', coreW ,pipeOneH+pipeOneW, dpt*z );

end

node=nodeBase + count;

node(2) = node(5);
node(3) = node(7);
node(7) = node(8);

node(1) =matchOne(2);
node(4) =matchOne(3);
node(5) = matchOne(6);
node(8) =matchOne(7);

fprintf(fileTwo, ' hex ( %i %i %i %i %i %i %i %i) (%i 1 %i) simpleGrading (1 1 1) \n', ...
               node(1), node(2), node(3), node(4), node(5), node(6), node(7), node(8), ...
                floor(pipeL/deltax), floor(pipeOneW/deltay)  );
            
%write wall faces
fprintf(fileThree, '(%i %i %i %i) \n', node(5), node(6), node(7),node(8))
fprintf(fileThree, '(%i %i %i %i) \n', node(1), node(2), node(3),node(4))

%write empty faces
fprintf(fileFour, '(%i %i %i %i) \n', node(1), node(5), node(6),node(2))
fprintf(fileFour, '(%i %i %i %i) \n', node(4),node(8), node(7),node(3))

%write inlet particles
fprintf(fileSix, '(%i %i %i %i) \n', node(2), node(3), node(7),node(6))

count = count +4;

%write block 2
for z=0:1
    
   % fprintf(file, ' ( %f %f %f) \n', coreW ,pipeTwoH, dpt*z );
      fprintf(file, ' ( %f %f %f) \n', coreW+pipeL ,pipeTwoH, dpt*z );  
    fprintf(file, ' ( %f %f %f) \n', coreW+pipeL ,pipeTwoH+pipeTwoW, dpt*z );

 %   fprintf(file, ' ( %f %f %f) \n', coreW ,pipeTwoH+pipeTwoW, dpt*z );

end

node=nodeBase + count;

node(2) = node(5);
node(3) = node(7);
node(7) = node(8);

node(1) =matchTwo(2);
node(4) =matchTwo(3);
node(5) = matchTwo(6);
node(8) =matchTwo(7);

fprintf(fileTwo, ' hex ( %i %i %i %i %i %i %i %i) (%i 1 %i) simpleGrading (1 1 1) \n', ...
               node(1), node(2), node(3), node(4), node(5), node(6), node(7), node(8), ...
                floor(pipeL/deltax), floor(pipeTwoW/deltay)  );
%write wall faces
fprintf(fileThree, '(%i %i %i %i) \n', node(5), node(6), node(7),node(8))
fprintf(fileThree, '(%i %i %i %i) \n', node(1), node(2), node(3),node(4))


%write empty faces
fprintf(fileFour, '(%i %i %i %i) \n', node(1), node(5), node(6),node(2))
fprintf(fileFour, '(%i %i %i %i) \n', node(4),node(8), node(7),node(3))

%write  outlet
fprintf(fileFive, '(%i %i %i %i) \n', node(2), node(3), node(7),node(6))
            
count = count +8;



