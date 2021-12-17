1. the following must be installed:
	. c++ compiler (gcc)
	. povray
	. ffmpeg

2. ./run_me.sh


my code is configured for 3D collisions with periodic boundary conditions.
to add a pairwise attractive potential between particles go to the following lines in collisions.cpp (70-71):

        int_dt(r,v,dt,N);
        //int_dt(r,v,m,CG,dt,N);

comment out the first line and uncomment the second:
	
        //int_dt(r,v,dt,N);
        int_dt(r,v,m,CG,dt,N);


to change from periodic boundary conditions to hard boundary conditions, go to the following lines in collisions.cpp (75-78):

            cb3D_pbc(r,s,N);
            cc3D_pbc(r,v,s,rs,m,dt,N,c);
            //cb3D_hbc(r,v,s,rs,dt,N);
            //cc3D_hbc(r,v,s,rs,m,dt,N,c);

comment out the first two lines and uncomment the last two lines:

            //cb3D_pbc(r,s,N);
            //cc3D_pbc(r,v,s,rs,m,dt,N,c);
            cb3D_hbc(r,v,s,rs,dt,N);
            cc3D_hbc(r,v,s,rs,m,dt,N,c);

 
*note: if you use periodic boundary conditions with the pairwise attractive potential integrator, interactions across the boundary will not be considered.

