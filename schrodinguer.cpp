//PROGRAMA QUE RESUELVE LA ECUACION DE SCHRODINGUER PARA UNA BARRERA DE POTENCIAL CUYA ALTURA
//PODREMOS IR VARIANDO CON EL PARÁMETRO LAMBDA

#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

using namespace std;

static const double PI= 3.14159265359; //Definicion de pi
static const double h=6.62607015e-34; //Cte de Planck
static const complex<double> i(0.0,1.0);//Definimos la unidad imaginaria


void generador (double s, double V[], double k0,int N,double lambda,complex<double> phi[],complex<double>a[]);
void Beta (double s, double V[],int N,complex <double> b[],complex<double>phi[],complex<double>a[]);
void Chi (int N,complex<double>a[],complex<double>b[],complex<double>chi[]);
void Phi (int N,complex<double>chi[],complex<double>phi[],double modulo[]);


int main()
{
    double lambda,s,k0,norma;
    int ciclos,N;
    double *V,*modulo; //Defino el potencial y el módulo con memoria dinámica
    int j,t=0;//Indices: el indice i queda reservado ya que se llama asi a la unidad imaginaria
    ofstream fich2,fich3;//Definimos los ficheros que se necesitarán
    stringstream pot;
    string pot2;

    //Introducción de N
    cout<<"Introduzca N"<<endl;
    cin>>N;

    //Vectores necesarios para el cálculo
    complex<double> phi[N];
    complex<double> a[N-1];
    complex<double> b[N-1];
    complex<double> chi[N];

    /** Recogida de parametros, inicialización de memoria dínamica y de algunas variables iniciales**/

    //Inicialización de la memoria dínamica de los vectores

    V=new double [N];
    modulo=new double [N];

    //Introducción de la longitud de onda
    cout<<"Introduzca el valor de lambda"<<endl;
    cin>>lambda;

    fich2.open("Norma.txt");
    fich3.open("Probabilidad.txt");

    //Introducción del número de ciclos
    do
    {
        cout<<"Introduzca el numero de ciclos a calcular entre 1 y "<<N/4 <<endl;
        cin>>ciclos;
    }
    while ((ciclos<1) || (ciclos >N/4));



    //Calculos necesarios para generar variables y condiciones
    k0=2.0*PI*ciclos/N;
    s=1.0/(4.0*k0*k0);

    //Generación de las funciones de onda iniciales y las alfas
    generador(s,V,k0,N,lambda,phi,a);

    //Calculo la norma inicial
    norma=0.0;
    for(j=0;j<N;j++)
        modulo[j]=norm(phi[j]);
    for(j=0;j<N;j++)
        norma+=modulo[j];



    //Este paso que hacemos a continuación no es necesario estrictamente pero es para tener el fichero
    //con norma 1 en todos los pasos
    for(j=0;j<N;j++)
    {
        modulo[j]=norm(phi[j]);
        fich3<<j<<" "<<modulo[j]<<"\t";
    }
    fich3<<endl;
    norma=0.0;
    for(j=0;j<N;j++)
        norma+=modulo[j];
    fich2<<norma<<endl;




    //Ahora empezamos el algoritmo iterativo
    while (t<10*ciclos)
    {

        //Step 2: Calculamos beta
        Beta(s,V,N,b,phi,a);

        //Step 3: Calculamos chi
        Chi(N,a,b,chi);

        //Step 4: Actualización de las funciones de onda
        Phi(N,chi,phi,modulo);

        //Se guardan los modulos de las funciones de onda al cuadrado para su futura representación
        if(t%10==0)
        {
            for(j=0;j<N;j++)
                fich3<<j<<" "<<modulo[j]<<"\t";
            fich3<<endl;
        }

        //A modo de comprobación vemos si la norma se mantiene constante
        norma=0.0;

        for(j=0;j<N;j++)
            norma+=modulo[j];
        fich2<<norma<<endl;

        t++; //Aumentamos el tiempo en una unidad
    }

    /** Borrado de la memoria dínamica y fin del programa **/

    pot<<V[(2*N/5)+1];
    pot>>pot2;


    //Borrado de memoria dínamica
    delete[] V;
    delete[] modulo;

    //Cierre de ficheros
    fich2.close();
    fich3.close();



    return 0;
}

void generador (double s, double V[], double k0,int N,double lambda,complex<double> phi[],complex<double>a[])
{
    int j;
    complex<double>Phi;

    //Condiciones de contorno
    phi[0]=phi[N-1]=complex<double>(0.0,0.0);
    V[0]=V[N-1]=0.0;
    a[N-1]=complex<double>(0.0);

    //Calculo de la función de onda
    for(j=1;j<(N-1);j++)
    {
        phi[j]=exp(j*k0*i)*exp(-8.0*(4.0*j-N)*(4.0*j-N)/(N*N*1.0));
    }

    //Inicialización del potencial
    for (j=0;j<N;j++)
    {
        if ((j>=(2.0*N/5.0)) && (j<=(3.0*N/5.0)))
            V[j]=(lambda*k0*k0);
        else
            V[j]=0.0;
    }

    //Inicialización de alfa
    for(j=(N-2);j>0;j--)
    {
        a[j]=-1.0/(-2.0-V[j+1]+a[j+1]+(2.0*i)/s);
    }

    return;
}

void Beta (double s, double V[],int N,complex <double> b[],complex<double>phi[],complex<double>a[])
{
    int j;

    //Condiciones de contorno
    b[N-1]=complex<double>(0.0,0.0);

    for(j=(N-2);j>0;j--)
    {
        b[j]=(((4.0*i*phi[j+1])/s)-b[j+1])/(-2.0-V[j+1]+a[j+1]+(2.0*i)/s);
    }

}

void Chi (int N,complex<double>a[],complex<double>b[],complex<double>chi[])
{
    int j;

    //Condiciones de contorno
    chi[0]=complex<double>(0.0,0.0);
    chi[N-1]=complex<double>(0.0,0.0);

    for(j=0;j<(N-2);j++)
    {
        chi[j+1]=a[j]*chi[j]+b[j];
    }
    return;
}

void Phi (int N,complex<double>chi[],complex<double>phi[],double modulo[])
{
    int j;
    for(j=0;j<N;j++)
    {
        phi[j]=chi[j]-phi[j];
        modulo[j]=norm(phi[j]);
    }
    return;
}

