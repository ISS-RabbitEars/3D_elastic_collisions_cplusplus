#include "colors.inc"

background{color Black}

#declare sx=5;
#declare sy=5;
#declare sz=5;

#fopen ndt "cam_temp.dat" read
#read(ndt,temp)
#declare cam=temp;


camera{
    location cam
    look_at <0,0,0>
    right <-4/3,0,0>
    sky<0,0,1>
}

light_source {2*cam color White}
light_source {-2*cam color White}

#fopen ndt "N.dat" read
#read(ndt,temp)
#declare N=temp;


#declare r=array[N];
#declare rs=array[N];

#fopen Rsd "rs.dat" read
#declare i=0;
#while (defined(Rsd))
    #read(Rsd,rs1)
    #declare rs[i]=rs1;
    #declare i=i+1;
#end

#fopen Pts "temp.dat" read
#declare i=0;
#while (defined(Pts))
    #read(Pts,vec1)
    #declare r[i]=vec1;
    #declare i=i+1;
#end

box{
    <-sx,-sy,-sz>,<sx,sy,sz>
    texture { pigment { color Blue transmit 0.9 } }
}

#declare j=0;
#while(j < i)
sphere {
  r[j] , rs[j]
  bounded_by { box {<-sx,-sy,-sz>,<sx,sy,sz>} }
  clipped_by { bounded_by }
  pigment { color Red }
}
    #if (r[j].x-rs[j]<-sx)
    sphere {
        <r[j].x+2*sx,r[j].y,r[j].z>, rs[j]
        bounded_by { box {<-sx,-sy,-sz>,<sx,sy,sz>} }
        clipped_by { bounded_by }
        pigment { color Red }
    }
    #end
    #if (r[j].x+rs[j]>sx)
    sphere {
        <r[j].x-2*sx,r[j].y,r[j].z>, rs[j]
        bounded_by { box {<-sx,-sy,-sz>,<sx,sy,sz>} }
        clipped_by { bounded_by }
        pigment { color Red }
    }
    #end
    #if (r[j].y-rs[j]<-sy)
    sphere {
        <r[j].x,r[j].y+2*sy,r[j].z>, rs[j]
        bounded_by { box {<-sx,-sy,-sz>,<sx,sy,sz>} }
        clipped_by { bounded_by }
        pigment { color Red }
    }
    #end
    #if (r[j].y+rs[j]>sy)
    sphere {
        <r[j].x,r[j].y-2*sy,r[j].z>, rs[j]
        bounded_by { box {<-sx,-sy,-sz>,<sx,sy,sz>} }
        clipped_by { bounded_by }
        pigment { color Red }
    }
    #end
    #if (r[j].z-rs[j]<-sz)
    sphere {
        <r[j].x,r[j].y,r[j].z+2*sz>, rs[j]
        bounded_by { box {<-sx,-sy,-sz>,<sx,sy,sz>} }
        clipped_by { bounded_by }
        pigment { color Red }
    }
    #end
    #if (r[j].z+rs[j]>sz)
    sphere {
        <r[j].x,r[j].y,r[j].z-2*sz>, rs[j]
        bounded_by { box {<-sx,-sy,-sz>,<sx,sy,sz>} }
        clipped_by { bounded_by }
        pigment { color Red }
    }
    #end
#declare j=j+1;
#end

