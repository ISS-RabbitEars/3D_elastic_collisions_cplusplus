#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "ran3.h"
#include "vec.h"
#include "list.h"
using namespace std;
#define PI 3.141592653589793238462643383

void init_pos3D(vec*,vec&,double*,int,int);
void init_vel3D(vec*,double,double,int,int);
void init_mrs(double*,double*,double,double,int,int);
void int_dt(vec*,vec*,double,int);
void cb3D_pbc(vec*,vec&,int);
void cb3D_pbc(vec&,vec&);
void cc3D_pbc(vec*,vec*,vec&,double*,double*,double,int,bool);
void cb3D_hbc(vec*,vec*,vec&,double*,double,int);
void cb3D_hbc(vec&,vec&,vec&,double&,double);
void cc3D_hbc(vec*,vec*,vec&,double*,double*,double,int,bool);
vec af(vec*,double*,double,int);
void int_dt(vec*,vec*,double*,double,double,int);

int main()
{
    int seed=3447342,N=100,fps,ns,fr;
	vec r[N],s,v[N];
	double rs[N],va,vb,m[N],ma,mb,dt,a,GC;
    bool c;
    ofstream file;
    
    a=5;
    s(a,a,a);
    va=1;
    vb=5;
    ma=3;
    mb=5;
    GC=1;
    dt=0.01;
    fps=60;
    ns=60;
    fr=fps*ns;
    c=true;
    
    init_mrs(m,rs,ma,mb,N,seed);
    init_pos3D(r,s,rs,N,seed);
    init_vel3D(v,va,vb,N,seed);
    
    file.open("pos/N.dat");
    file<<N<<endl;
    file.close();
    
    file.open("pos/nframes.dat");
    file<<fr<<endl;
    file.close();
    
    file.open("pos/rs.dat");
    for(int i=0;i<N;i++)
    {
        file<<rs[i]<<endl;
        if(i<N-1) file<<",";
        file<<endl;
    }
    file.close();
    
    string pts="pos/pos_";
    for(int i=0;i<fr;i++)
    {
        int_dt(r,v,dt,N);
        //int_dt(r,v,m,GC,dt,N);
        do
        {
            c=false;
            cb3D_pbc(r,s,N);
            cc3D_pbc(r,v,s,rs,m,dt,N,c);
            //cb3D_hbc(r,v,s,rs,dt,N);
            //cc3D_hbc(r,v,s,rs,m,dt,N,c);
        } while(c);
        
        string pu,snum;
        stringstream ssnum;
        ssnum<<i+1;
        snum=ssnum.str();
        ssnum.str("");
        pu=pts+snum;
        file.open(pu.c_str());
        for(int j=0;j<N;j++)
        {
            file<<r[j];
            if(j<N-1) file<<",";
            file<<endl;
        }
        file.close();
    }
    
	return 0;
}

void init_pos3D(vec* r, vec& s, double* rs, int N, int seed)
{
    vec rt,dr;
    bool col;
    double ds;
    
    r[0]((s.x-rs[0])*(1-2*ran3(&seed)),(s.y-rs[0])*(1-2*ran3(&seed)),(s.z-rs[0])*(1-2*ran3(&seed)));
    for(int i=1;i<N;i++)
    {
        col=false;
        do
        {
            rt((s.x-rs[i])*(1-2*ran3(&seed)),(s.y-rs[i])*(1-2*ran3(&seed)),(s.z-rs[0])*(1-2*ran3(&seed)));
            int j=0;
            do
            {
                col=false;
                dr=rt-r[j];
                ds=rs[i]+rs[j];
                if(dr.mag()<ds) {col=true;}
                j++;
            }while(j<i&&!col);
        }while(col);
        r[i]=rt;
    }
}

void init_vel3D(vec* v,double a,double b, int N, int seed)
{
    double alpha,theta,phi,vi;
    
    alpha=b-a;
    for(int i=0;i<N;i++)
    {
        theta=2*PI*ran3(&seed);
        phi=PI*ran3(&seed);
        vi=alpha*ran3(&seed)+a;
        v[i](vi*cos(theta)*sin(phi),vi*sin(theta)*sin(phi),vi*cos(phi));
    }
}

void init_mrs(double* m, double* rs, double a, double b, int N, int seed)
{
    double alpha=b-a;
    for(int i=0;i<N;i++)
    {
        m[i]=alpha*ran3(&seed)+a;
        rs[i]=m[i]/(2*b);
    }
}

void int_dt(vec* r, vec* v, double dt, int N)
{
    for(int i=0;i<N;i++)
    {
        r[i]+=v[i]*dt;
    }
}

void cb3D_pbc(vec* r, vec& s, int N)
{
    for(int i=0;i<N;i++)
    {
        cb3D_pbc(r[i],s);
    }
}

void cb3D_pbc(vec& r, vec& s)
{
    if(r.x<-s.x) {r.x+=2*s.x;}
    if(r.x>=s.x) {r.x-=2*s.x;}
    if(r.y<-s.y) {r.y+=2*s.y;}
    if(r.y>=s.y) {r.y-=2*s.y;}
    if(r.z<-s.z) {r.z+=2*s.z;}
    if(r.z>=s.z) {r.z-=2*s.z;}
}

void cc3D_pbc(vec* r, vec* v, vec& s, double* rs, double* m, double t, int n, bool c)
{
    for(int i=0;i<n;i++)
    {
        list drc,drj,drlmi;
        bool col=false;
        for(int j=0;j<n;j++)
        {
            if(j!=i)
            {
                vec dr0,dr1,dr2,temp;
                list drl;
                
                dr0=r[j]-r[i];
                dr1=dr0+2*s;
                dr2=dr0-2*s;
                drl.append(dr0.mag());
                temp(dr0.x,dr0.y,dr1.z);
                drl.append(temp.mag());
                temp(dr0.x,dr0.y,dr2.z);
                drl.append(temp.mag());
                temp(dr0.x,dr1.y,dr0.z);
                drl.append(temp.mag());
                temp(dr0.x,dr1.y,dr1.z);
                drl.append(temp.mag());
                temp(dr0.x,dr1.y,dr2.z);
                drl.append(temp.mag());
                temp(dr0.x,dr2.y,dr0.z);
                drl.append(temp.mag());
                temp(dr0.x,dr2.y,dr1.z);
                drl.append(temp.mag());
                temp(dr0.x,dr2.y,dr2.z);
                drl.append(temp.mag());
                temp(dr1.x,dr0.y,dr0.z);
                drl.append(temp.mag());
                temp(dr1.x,dr0.y,dr1.z);
                drl.append(temp.mag());
                temp(dr1.x,dr0.y,dr2.z);
                drl.append(temp.mag());
                temp(dr1.x,dr1.y,dr0.z);
                drl.append(temp.mag());
                drl.append(dr1.mag());
                temp(dr1.x,dr1.y,dr2.z);
                drl.append(temp.mag());
                temp(dr1.x,dr2.y,dr0.z);
                drl.append(temp.mag());
                temp(dr1.x,dr2.y,dr1.z);
                drl.append(temp.mag());
                temp(dr1.x,dr2.y,dr2.z);
                drl.append(temp.mag());
                temp(dr2.x,dr0.y,dr0.z);
                drl.append(temp.mag());
                temp(dr2.x,dr0.y,dr1.z);
                drl.append(temp.mag());
                temp(dr2.x,dr0.y,dr2.z);
                drl.append(temp.mag());
                temp(dr2.x,dr1.y,dr0.z);
                drl.append(temp.mag());
                temp(dr2.x,dr1.y,dr1.z);
                drl.append(temp.mag());
                temp(dr2.x,dr1.y,dr2.z);
                drl.append(temp.mag());
                temp(dr2.x,dr2.y,dr0.z);
                drl.append(temp.mag());
                temp(dr2.x,dr2.y,dr1.z);
                drl.append(temp.mag());
                drl.append(dr2.mag());
                if(drl.min<=rs[i]+rs[j])
                {
                    drc.append(drl.min);
                    drj.append(j);
                    drlmi.append(drl.mini+1);
                    col=true;
                    c=true;
                }
            }
        }
        if(col)
        {
            int j,id;
            vec dr,dv,drh,vri,vrj;
            double qa,qb,qc,dt1,dt2,tr,pr,sm,dm,vifx,vjfx;
            
            j=drj.a[drc.mini];
            id=drlmi.a[drc.mini];
            r[i]-=v[i]*t;
            r[j]-=v[j]*t;
            switch(id)
            {
                case 1: dr=r[j]-r[i]; break;
                case 2: dr=r[j]-r[i]+vec(0,0,2*s.z); break;
                case 3: dr=r[j]-r[i]+vec(0,0,-2*s.z); break;
                case 4: dr=r[j]-r[i]+vec(0,2*s.y,0); break;
                case 5: dr=r[j]-r[i]+vec(0,2*s.y,2*s.z); break;
                case 6: dr=r[j]-r[i]+vec(0,2*s.y,-2*s.z); break;
                case 7: dr=r[j]-r[i]+vec(0,-2*s.y,0); break;
                case 8: dr=r[j]-r[i]+vec(0,-2*s.y,2*s.z); break;
                case 9: dr=r[j]-r[i]+vec(0,-2*s.y,-2*s.z); break;
                case 10: dr=r[j]-r[i]+vec(2*s.x,0,0); break;
                case 11: dr=r[j]-r[i]+vec(2*s.x,0,2*s.z); break;
                case 12: dr=r[j]-r[i]+vec(2*s.x,0,-2*s.z); break;
                case 13: dr=r[j]-r[i]+vec(2*s.x,2*s.y,0); break;
                case 14: dr=r[j]-r[i]+2*s; break;
                case 15: dr=r[j]-r[i]+vec(2*s.x,2*s.y,-2*s.z); break;
                case 16: dr=r[j]-r[i]+vec(2*s.x,-2*s.y,0); break;
                case 17: dr=r[j]-r[i]+vec(2*s.x,-2*s.y,2*s.z); break;
                case 18: dr=r[j]-r[i]+vec(2*s.x,-2*s.y,-2*s.z); break;
                case 19: dr=r[j]-r[i]+vec(-2*s.x,0,0); break;
                case 20: dr=r[j]-r[i]+vec(-2*s.x,0,2*s.z); break;
                case 21: dr=r[j]-r[i]+vec(-2*s.x,0,-2*s.z); break;
                case 22: dr=r[j]-r[i]+vec(-2*s.x,2*s.y,0); break;
                case 23: dr=r[j]-r[i]+vec(-2*s.x,2*s.y,2*s.z); break;
                case 24: dr=r[j]-r[i]+vec(-2*s.x,2*s.y,-2*s.z); break;
                case 25: dr=r[j]-r[i]+vec(-2*s.x,-2*s.y,0); break;
                case 26: dr=r[j]-r[i]+vec(-2*s.x,-2*s.y,2*s.z); break;
                case 27: dr=r[j]-r[i]-2*s; break;
            }
            dv=v[j]-v[i];
            qa=dv.mag2();
            qb=2*(dr*dv);
            qc=dr.mag2()-pow(rs[i]+rs[j],2);
            dt1=(-qb-sqrt(qb*qb-4*qa*qc))/(2*qa);
            dt2=t-dt1;
            r[i]+=v[i]*dt1;
            r[j]+=v[j]*dt1;
            switch(id)
            {
                case 1: dr=r[j]-r[i]; break;
                case 2: dr=r[j]-r[i]+vec(0,0,2*s.z); break;
                case 3: dr=r[j]-r[i]+vec(0,0,-2*s.z); break;
                case 4: dr=r[j]-r[i]+vec(0,2*s.y,0); break;
                case 5: dr=r[j]-r[i]+vec(0,2*s.y,2*s.z); break;
                case 6: dr=r[j]-r[i]+vec(0,2*s.y,-2*s.z); break;
                case 7: dr=r[j]-r[i]+vec(0,-2*s.y,0); break;
                case 8: dr=r[j]-r[i]+vec(0,-2*s.y,2*s.z); break;
                case 9: dr=r[j]-r[i]+vec(0,-2*s.y,-2*s.z); break;
                case 10: dr=r[j]-r[i]+vec(2*s.x,0,0); break;
                case 11: dr=r[j]-r[i]+vec(2*s.x,0,2*s.z); break;
                case 12: dr=r[j]-r[i]+vec(2*s.x,0,-2*s.z); break;
                case 13: dr=r[j]-r[i]+vec(2*s.x,2*s.y,0); break;
                case 14: dr=r[j]-r[i]+2*s; break;
                case 15: dr=r[j]-r[i]+vec(2*s.x,2*s.y,-2*s.z); break;
                case 16: dr=r[j]-r[i]+vec(2*s.x,-2*s.y,0); break;
                case 17: dr=r[j]-r[i]+vec(2*s.x,-2*s.y,2*s.z); break;
                case 18: dr=r[j]-r[i]+vec(2*s.x,-2*s.y,-2*s.z); break;
                case 19: dr=r[j]-r[i]+vec(-2*s.x,0,0); break;
                case 20: dr=r[j]-r[i]+vec(-2*s.x,0,2*s.z); break;
                case 21: dr=r[j]-r[i]+vec(-2*s.x,0,-2*s.z); break;
                case 22: dr=r[j]-r[i]+vec(-2*s.x,2*s.y,0); break;
                case 23: dr=r[j]-r[i]+vec(-2*s.x,2*s.y,2*s.z); break;
                case 24: dr=r[j]-r[i]+vec(-2*s.x,2*s.y,-2*s.z); break;
                case 25: dr=r[j]-r[i]+vec(-2*s.x,-2*s.y,0); break;
                case 26: dr=r[j]-r[i]+vec(-2*s.x,-2*s.y,2*s.z); break;
                case 27: dr=r[j]-r[i]-2*s; break;
            }
            drh=dr.hat();
            tr=acos(drh.x/sqrt(drh.x*drh.x+drh.y*drh.y));
            if(drh.y<0) {tr=-tr;}
            pr=asin(drh.z);
            vri=v[i].rotate(tr,pr);
            vrj=v[j].rotate(tr,pr);
            sm=m[i]+m[j];
            dm=m[j]-m[i];
            vifx=(-dm/sm)*vri.x+(2*m[j]/sm)*vrj.x;
            vjfx=(2*m[i]/sm)*vri.x+(dm/sm)*vrj.x;
            vri.x=vifx;
            vrj.x=vjfx;
            v[i]=vri.rotate(-tr,-pr);
            v[j]=vrj.rotate(-tr,-pr);
            r[i]+=v[i]*dt2;
            r[j]+=v[j]*dt2;
            cb3D_pbc(r[i],s);
            cb3D_pbc(r[j],s);
        }
    }
}

void cb3D_hbc(vec* r, vec* v, vec& s, double* rs, double t, int n)
{
    for(int i=0;i<n;i++)
    {
        cb3D_hbc(r[i],v[i],s,rs[i],t);
    }
}

void cb3D_hbc(vec& r, vec& v, vec& s, double& rs, double t)
{
    double dt1,dt2;
    if(r.x-rs<-s.x)
    {
        r-=v*t;
        dt1=(-s.x+rs-r.x)/v.x;
        dt2=t-dt1;
        r+=v*dt1;
        v.x=-v.x;
        r+=v*dt2;
    }
    if(r.x+rs>=s.x)
    {
        r-=v*t;
        dt1=(s.x-rs-r.x)/v.x;
        dt2=t-dt1;
        r+=v*dt1;
        v.x=-v.x;
        r+=v*dt2;
    }
    if(r.y-rs<-s.y)
    {
        r-=v*t;
        dt1=(-s.y+rs-r.y)/v.y;
        dt2=t-dt1;
        r+=v*dt1;
        v.y=-v.y;
        r+=v*dt2;
    }
    if(r.y+rs>=s.y)
    {
        r-=v*t;
        dt1=(s.y-rs-r.y)/v.y;
        dt2=t-dt1;
        r+=v*dt1;
        v.y=-v.y;
        r+=v*dt2;
    }
    if(r.z-rs<-s.z)
    {
        r-=v*t;
        dt1=(-s.z+rs-r.z)/v.z;
        dt2=t-dt1;
        r+=v*dt1;
        v.z=-v.z;
        r+=v*dt2;
    }
    if(r.z+rs>=s.z)
    {
        r-=v*t;
        dt1=(s.z-rs-r.z)/v.z;
        dt2=t-dt1;
        r+=v*dt1;
        v.z=-v.z;
        r+=v*dt2;
    }
}

void cc3D_hbc(vec* r, vec* v, vec& s, double* rs, double* m, double t, int n, bool c)
{
    for(int i=0;i<n;i++)
    {
        list drc,drj;
        bool col=false;
        for(int j=0;j<n;j++)
        {
            if(j!=i)
            {
                vec drt;
                double drtm,rss;
                
                drt=r[j]-r[i];
                drtm=drt.mag();
                rss=rs[i]+rs[j];
                if(drtm<=rss)
                {
                    drc.append(drtm/rss);
                    drj.append(j);
                    col=true;
                    c=true;
                }
            }
        }
        if(col)
        {
            int j=drj.a[drc.mini];
            vec dr,dv,drh,vri,vrj;
            double qa,qb,qc,dt1,dt2,tr,pr,sm,dm,vifx,vjfx;
            
            r[i]-=v[i]*t;
            r[j]-=v[j]*t;
            dr=r[j]-r[i];
            dv=v[j]-v[i];
            qa=dv.mag2();
            qb=2*(dr*dv);
            qc=dr.mag2()-pow(rs[i]+rs[j],2);
            dt1=(-qb-sqrt(qb*qb-4*qa*qc))/(2*qa);
            dt2=t-dt1;
            r[i]+=v[i]*dt1;
            r[j]+=v[j]*dt1;
            drh=dr.hat();
            tr=acos(drh.x/sqrt(drh.x*drh.x+drh.y*drh.y));
            if(drh.y<0) {tr=-tr;}
            pr=asin(drh.z);
            vri=v[i].rotate(tr,pr);
            vrj=v[j].rotate(tr,pr);
            sm=m[i]+m[j];
            dm=m[j]-m[i];
            vifx=(-dm/sm)*vri.x+(2*m[j]/sm)*vrj.x;
            vjfx=(2*m[i]/sm)*vri.x+(dm/sm)*vrj.x;
            vri.x=vifx;
            vrj.x=vjfx;
            v[i]=vri.rotate(-tr,-pr);
            v[j]=vrj.rotate(-tr,-pr);
            r[i]+=v[i]*dt2;
            r[j]+=v[j]*dt2;
            cb3D_hbc(r[i],v[i],s,rs[i],t);
            cb3D_hbc(r[j],v[j],s,rs[j],t);
        }
    }
}

vec af(int i, vec* r, double* m, double G, int n)
{
    vec dr,sum;
    for(int j=0;j<n;j++)
    {
        if(j!=i)
        {
            dr=r[j]-r[i];
            sum+=(G*m[j]/dr.mag2())*dr.hat();
        }
    }
    return sum;
}

void int_dt(vec* r, vec* v, double* m, double G, double dt, int n)
{
    for(int i=0;i<n;i++)
    {
        v[i]+=af(i,r,m,G,n)*dt;
        r[i]+=v[i]*dt;
    }
}

#undef PI
