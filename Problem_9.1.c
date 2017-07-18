#include <stdio.h>
#include <math.h>

int main()
{
    double alpha[20],g[20],a[20],f[20],b[20],c[20],u[20],e[20];
    a[0]=2.04;a[1]=2.04;a[2]=2.04;a[3]=2.04;
    f[0]=40.8;f[1]=0.8;f[2]=0.8;f[3]=200.8;
    b[1]=-1.0;b[2]=-1.0;b[3]=-1.0;
    c[0]=-1.0,c[1]=-1.0;c[2]=-1.0;
    alpha[0]=a[0];g[0]=f[0];
    int num=0;

    printf("Enter number of columns in matrix: ");
    scanf("%d",&num);

    for(int j=1;j<num+1;j++)
    {
        a[j]=a[j]-((b[j]*c[j-1])/a[j-1]);
        f[j]=f[j]-(b[j]/a[j-1])*f[j-1];
    }

    u[num]=f[num]/a[num];

    for(int k=1;k<num-1;k++)
    {
        u[num-k]=(f[num-k]-(c[num-k]*u[num-k+1]))/a[num-k];
    }

    for(int h=num-1;h>1;h--)
    {
        u[h]=(g[h]-(c[h]*u[h+1]))/a[h];
    }

    for(int i=0;i<num;i++)
    {
        e[i]=(fabs(b[i]*u[i-1]+a[i]*u[i]+c[i]*u[i+1]-f[i]))/(fabs(b[i]*u[i-1])+fabs(a[i]*u[i])+fabs(c[i]*u[i+1])+fabs(f[i]));
        printf("\n\n%lf\n\n",e[i]);
    }
    return 0;
}
