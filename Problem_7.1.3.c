#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define M1 0.08*1.989*pow(10,30)
#define M2 0.85*5.9722*pow(10,24)
#define G 6.67*pow(10,-11)

double f_r12(double x1, double x2, double y1, double y2)
{
    double v;
    v=sqrt(((x1-x2)*(x1-x2))+((y1-y2)*(y1-y2)));
    return v;
}

double f_r1(double x1,double y1,double x2,double y2)
{
    double x;
    x=-1*((G*M2*(x1-x2))/pow(f_r12(x1,x2,y1,y2),3));

    return x;
}

double f_r2(double x1,double y1,double x2,double y2)
{
    double x;
    x=-1*((G*M1*(x2-x1))/pow(f_r12(x1,x2,y1,y2),3));

    return x;
}

double a(double v)
{
    return v;
}

double dtChange(double dt,double x1,double x2,double y1,double y2,double x11,double x21,double y11,double y21,
                double v1_x,double v2_x,double v1_y,double v2_y,double v1_x1,double v2_x1,double v1_y1,double v2_y1)
{
    double u[8],w[8], e = pow(10,-8), deltadesired,d1,d2,d3,d4, phi, phi2;

    u[0] = fabs(x11-x1);        w[0]=x11;
    u[1] = fabs(x21-x2);        w[1]=x21;
    u[2] = fabs(y11-y1);        w[2]=y11;
    u[3] = fabs(y21-y2);        w[3]=y21;

    u[4] = fabs(v1_x1-v1_x);    w[4]=v1_x1;
    u[5] = fabs(v1_y1-v1_y);    w[5]=v1_y1;
    u[6] = fabs(v2_x1-v2_x);    w[6]=v2_x1;
    u[7] = fabs(v2_y1-v2_y);    w[7]=v2_y1;

    double max= -1000000000;
    int place=0;
    for (int t=0;t<3;t++)
    {
        if (u[t]>max)
            max=u[t];
            place=t;
    }

    //deltadesired = e*(fabs(w[place])+dt)/abs(u[place+4]);
    deltadesired=e*(fabs(max)+dt)/abs(u[place+4]);
    //printf("deltadesired: %lf",deltadesired);
    /*
    d1=(e*(fabs(x11)+dt*fabs(v1_x1)))/(fabs(v1_x1-v1_x));
    d2=(e*(fabs(y11)+dt*fabs(v1_y1)))/(fabs(v1_y1-v1_y));
    d3=(e*(fabs(x21)+dt*fabs(v2_x1)))/(fabs(v2_x1-v2_x));
    d4=(e*(fabs(y21)+dt*fabs(v2_y1)))/(fabs(v2_y1-v2_y));

    //printf("%lf %lf %lf %lf\n",d1,d2,d3,d4);

    if(d1>d2&&d1>d3&&d1>d4)
        phi=d1;
    else if(d2>d1&&d2>d3&&d2>d4)
        phi=d2;
    else if(d3>d1&&d3>d2&&d3>d4)
        phi=d3;
    else if(d4>d1&&d4>d2&&d4>d3)
        phi=d4;*/

    phi=dt+deltadesired/max;

    /*if(phi>1)
        dt=dt*pow(phi,.2);

    if(phi<1)
        dt=dt*pow(phi,.25);*/

    if(max<deltadesired)
        dt=dt*pow(phi,.2);

    if(max>deltadesired)
        dt=dt*pow(phi,.25);

    return(dt);
}

double RK(double max,FILE *fptr)
{
    double X1star,X1starstar,X1starstarstar,Y1star,Y1starstar,Y1starstarstar;
    double X2star,X2starstar,X2starstarstar,Y2star,Y2starstar,Y2starstarstar;
    double V1_x_star,V1_x_starstar,V1_x_starstarstar,V1_y_star,V1_y_starstar,V1_y_starstarstar;
    double V2_x_star,V2_x_starstar,V2_x_starstarstar,V2_y_star,V2_y_starstar,V2_y_starstarstar;
    double x1=0.0,y1=0.0,v1_x=0.0,v1_y=0.0, x2=1.9889*pow(10,8),y2=-1.651*pow(10,9),v2_x=79348.8,v2_y=9078.256,DT=60.0,DT1=DT,t=0.0,t1=0.0;

    fprintf(fptr,"%lf, %lf, %lf, %lf, %lf\n",t,x1,y1,x2,y2);
    printf("x1(%lf): %lf\t\ty1(%lf): %lf\t\tx2(%lf): %lf\t\ty2(%lf): %lf\n\n",t,x1,t,y1,t,x2,t,y2);

    for(t=DT;t<max-5*DT;t+=DT1)
    {
        // Step 1 Runge Kutta
        X1star=x1+DT/2*a(v1_x);
        X2star=x2+DT/2*a(v2_x);
        Y1star=y1+DT/2*a(v1_y);
        Y2star=y2+DT/2*a(v2_y);


        V1_x_star=v1_x+DT/2*f_r1(x1,y1,x2,y2);
        V1_y_star=v1_y+DT/2*f_r1(y1,x1,y2,x2);

        V2_x_star=v2_x+DT/2*f_r2(x1,y1,x2,y2);
        V2_y_star=v2_y+DT/2*f_r2(y1,x1,y2,x2);


        // Step 2 Runge Kutta
        X1starstar=x1+DT/2*a(V1_x_star);
        X2starstar=x2+DT/2*a(V2_x_star);
        Y1starstar=y1+DT/2*a(V1_y_star);
        Y2starstar=y2+DT/2*a(V2_y_star);

        V1_x_starstar=v1_x+DT/2*f_r1(X1star,Y1star,X2star,Y2star);
        V1_y_starstar=v1_y+DT/2*f_r1(Y1star,X1star,Y2star,X2star);

        V2_x_starstar=v2_x+DT/2*f_r2(X1star,Y1star,X2star,Y2star);
        V2_y_starstar=v2_y+DT/2*f_r2(Y1star,X1star,Y2star,X2star);


        // Step 3 Runge Kutta
        X1starstarstar=x1+DT*a(V1_x_starstar);
        X2starstarstar=x2+DT*a(V2_x_starstar);
        Y1starstarstar=y1+DT*a(V1_y_starstar);
        Y2starstarstar=y2+DT*a(V2_y_starstar);

        V1_x_starstarstar=v1_x+DT*f_r1(X1starstar,Y1starstar,X2starstar,Y2starstar);
        V1_y_starstarstar=v1_y+DT*f_r1(Y1starstar,X1starstar,Y2starstar,X2starstar);

        V2_x_starstarstar=v2_x+DT*f_r2(X1starstar,Y1starstar,X2starstar,Y2starstar);
        V2_y_starstarstar=v2_y+DT*f_r2(Y1starstar,X1starstar,Y2starstar,X2starstar);


        // Step 4 Runge Kutta
        x1=x1+(DT/6)*(a(v1_x)+2*a(V1_x_star)+2*a(V1_x_starstar)+a(V1_x_starstarstar));
        y1=y1+(DT/6)*(a(v1_y)+2*a(V1_y_star)+2*a(V1_y_starstar)+a(V1_y_starstarstar));
        x2=x2+(DT/6)*(a(v2_x)+2*a(V2_x_star)+2*a(V2_x_starstar)+a(V2_x_starstarstar));
        y2=y2+(DT/6)*(a(v2_y)+2*a(V2_y_star)+2*a(V2_y_starstar)+a(V2_y_starstarstar));

        v1_x=v1_x+(DT/6)*(f_r1(x1,y1,x2,y2)+2*f_r1(X1star,Y1star,X2star,Y2star)+2*f_r1(X1starstar,Y1starstar,X2starstar,Y2starstar)+
                           f_r1(X1starstarstar,Y1starstarstar,X2starstarstar,Y2starstarstar));
        v1_y=v1_y+(DT/6)*(f_r1(y1,x1,y2,x2)+2*f_r1(Y1star,X1star,Y2star,X2star)+2*f_r1(Y1starstar,X1starstar,Y2starstar,X2starstar)+
                           f_r1(Y1starstarstar,X1starstarstar,Y2starstarstar,X2starstarstar));

        v2_x=v2_x+(DT/6)*(f_r2(x1,y1,x2,y2)+2*f_r2(X1star,Y1star,X2star,Y2star)+2*f_r2(X1starstar,Y1starstar,X2starstar,Y2starstar)+
                           f_r2(X1starstarstar,Y1starstarstar,X2starstarstar,Y2starstarstar));
        v2_y=v2_y+(DT/6)*(f_r2(y1,x1,y2,x2)+2*f_r2(Y1star,X1star,Y2star,X2star)+2*f_r2(Y1starstar,X1starstar,Y2starstar,X2starstar)+
                           f_r2(Y1starstarstar,X1starstarstar,Y2starstarstar,X2starstarstar));

        double x11=x1;
            double x21=x2;
            double y11=y1;
            double y21=y2;

            double v1_x1=v1_x;
            double v1_y1=v1_y;
            double v2_x1=v2_x;
            double v2_y1=v2_y;

        for(t1=DT1;t1<(DT1*2);t1+=(DT1/2))
        {

            // Step 1 Runge Kutta
        X1star=x11+DT1/4*a(v1_x1);
        X2star=x21+DT1/4*a(v2_x1);
        Y1star=y11+DT1/4*a(v1_y1);
        Y2star=y21+DT1/4*a(v2_y1);


        V1_x_star=v1_x1+DT1/4*f_r1(x11,y11,x21,y21);
        V1_y_star=v1_y1+DT1/4*f_r1(y11,x11,y21,x21);

        V2_x_star=v2_x1+DT1/4*f_r2(x11,y11,x21,y21);
        V2_y_star=v2_y1+DT1/4*f_r2(y11,x11,y21,x21);


        // Step 2 Runge Kutta
        X1starstar=x11+DT1/4*a(V1_x_star);
        X2starstar=x21+DT1/4*a(V2_x_star);
        Y1starstar=y11+DT1/4*a(V1_y_star);
        Y2starstar=y21+DT1/4*a(V2_y_star);

        V1_x_starstar=v1_x1+DT1/4*f_r1(X1star,Y1star,X2star,Y2star);
        V1_y_starstar=v1_y1+DT1/4*f_r1(Y1star,X1star,Y2star,X2star);

        V2_x_starstar=v2_x1+DT1/4*f_r2(X1star,Y1star,X2star,Y2star);
        V2_y_starstar=v2_y1+DT1/4*f_r2(Y1star,X1star,Y2star,X2star);


        // Step 3 Runge Kutta
        X1starstarstar=x11+DT1/2*a(V1_x_starstar);
        X2starstarstar=x21+DT1/2*a(V2_x_starstar);
        Y1starstarstar=y11+DT1/2*a(V1_y_starstar);
        Y2starstarstar=y21+DT1/2*a(V2_y_starstar);

        V1_x_starstarstar=v1_x1+DT1/2*f_r1(X1starstar,Y1starstar,X2starstar,Y2starstar);
        V1_y_starstarstar=v1_y1+DT1/2*f_r1(Y1starstar,X1starstar,Y2starstar,X2starstar);

        V2_x_starstarstar=v2_x1+DT1/2*f_r2(X1starstar,Y1starstar,X2starstar,Y2starstar);
        V2_y_starstarstar=v2_y1+DT1/2*f_r2(Y1starstar,X1starstar,Y2starstar,X2starstar);


        // Step 4 Runge Kutta
        x11=x11+(DT1/12)*(a(v1_x1)+2*a(V1_x_star)+2*a(V1_x_starstar)+a(V1_x_starstarstar));
        y11=y11+(DT1/12)*(a(v1_y1)+2*a(V1_y_star)+2*a(V1_y_starstar)+a(V1_y_starstarstar));
        x21=x21+(DT1/12)*(a(v2_x1)+2*a(V2_x_star)+2*a(V2_x_starstar)+a(V2_x_starstarstar));
        y21=y21+(DT1/12)*(a(v2_y1)+2*a(V2_y_star)+2*a(V2_y_starstar)+a(V2_y_starstarstar));

        v1_x1=v1_x1+(DT1/12)*(f_r1(x11,y11,x21,y21)+2*f_r1(X1star,Y1star,X2star,Y2star)+2*f_r1(X1starstar,Y1starstar,X2starstar,Y2starstar)+
                           f_r1(X1starstarstar,Y1starstarstar,X2starstarstar,Y2starstarstar));
        v1_y1=v1_y1+(DT1/12)*(f_r1(y11,x11,y21,x21)+2*f_r1(Y1star,X1star,Y2star,X2star)+2*f_r1(Y1starstar,X1starstar,Y2starstar,X2starstar)+
                           f_r1(Y1starstarstar,X1starstarstar,Y2starstarstar,X2starstarstar));

        v2_x1=v2_x1+(DT1/12)*(f_r2(x11,y11,x21,y21)+2*f_r2(X1star,Y1star,X2star,Y2star)+2*f_r2(X1starstar,Y1starstar,X2starstar,Y2starstar)+
                           f_r2(X1starstarstar,Y1starstarstar,X2starstarstar,Y2starstarstar));
        v2_y1=v2_y1+(DT1/12)*(f_r2(y11,x11,y21,x21)+2*f_r2(Y1star,X1star,Y2star,X2star)+2*f_r2(Y1starstar,X1starstar,Y2starstar,X2starstar)+
                           f_r2(Y1starstarstar,X1starstarstar,Y2starstarstar,X2starstarstar));
        }

        fprintf(fptr,"%lf, %lf, %lf, %lf, %lf\n",t,x1,y1,x2,y2);
        printf("x1(%lf): %lf  y1(%lf): %lf  x2(%lf): %lf  y2(%lf): %lf  dt: %lf\n\n",t,x1,t,y1,t,x2,t,y2,DT1);
        DT1=dtChange(DT,x1,x2,y1,y2,x11,x21,y11,y21,v1_x,v2_x,v1_y,v2_y,v1_x1,v2_x1,v1_y1,v2_y1);
    }

}

int main()
{
    FILE *fptr;
    remove("Planet_Data.csv");
    fptr=fopen("Planet_Data.csv","w");
    fprintf(fptr,"t, x1, y1, x2, y2\n");
    printf("Problem 7.1.1\n\n");

    RK(1.509*(24*60*60),fptr);

    fclose(fptr);

    return 0;
}
