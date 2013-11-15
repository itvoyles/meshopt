function postnodes

faces = load( 'faces.dat' );
imax = length(faces);
cells = imax-1;
fidnode = fopen('nodes.dat','w');

fprintf(fidnode,'%d\n',cells);
for j = 1:imax
    fprintf(fidnode,'%17.15f\n',faces(j,1));
end

end