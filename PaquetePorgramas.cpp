#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iomanip> 

int const Tam=101;


int dim;
float norma(float vector1[],float vector2[]);
float suma_jacobi(float Matriz[], float vector[], int componente);



void PideDatos(int *Dim, float Mat[][Tam]);
void EscribeDatos(int Dim, float Mat[][Tam]);
void CalcDet(int Dim, float Mat[][Tam]);
void ResuelveGauss(int Dim, float Sist[][102]);
void PideDatos2(int *Dim, float Sist[][102]);


using namespace std;
double funcion (double x)
{
    return cos(x)*pow(x,2)-2*x;
}

double derivada (double x)
{
    return 2*x*(cos(2*x) - x*sin(2*x));
}

double funcion2 (double x)
{
    return 4/(pow(x,2)*5);
}

double derivada2 (double x)
{
    return -8/(5*(pow(x,3)));
}


double funcion3 (double x)
{
    return pow(x,3)-3*sin(pow(x,2))+1;
}

double derivada3 (double x)
{
    return 3*x*(x - 2*cos(pow(x,2)));
}

double funcion4 (double x)
{
    return pow(x,3) + 6*pow(x,2) + 9.4*x + 2.5;
}

double derivada4 (double x)
{
    return 3*pow(x,2) + 12*x + 9.4;
}



int main (int argc, char *argv[]) {
    int p;
    double a, b, tolerancia, a2, b2;
	cout<<"Paquete de programas "<<endl;
	cout<<"Profesora: Teresa Carrillo"<<endl;
	cout<<"Alumno: David Diaz Dorantes"<<endl;
	cout<<"Materia: Metodos numericos"<<endl;
	cout<<"Facultad de Estudios Superiores Acatlan"<<endl;	
	cout<<"Que programa desea usar? ------>  \n 1.-Programa 1 Metodo de Newton y Biseccion \n 2.-Programa 2 Calculo del determinate de una matriz \n 3.-Programa 3 Metodo de Jacobi \n 4.-Programa 4 Valores propios. Metodo de potencias  "<<endl; 
	cout<<"Para ejecutar el programa 1 ingrese un 1"<<endl;
	cout<<"Para ejecutar el programa 2 ingrese un 2"<<endl;
	cout<<"Para ejecutar el programa 3 ingrese un 3"<<endl;
	cout<<"Para ejecutar el programa 4 ingrese un 4"<<endl;
	cin >> p ;
	
	if(p == 1){
		
		cout <<"-------------------Primero se ejcutara el metodo de newton---------------"<<endl;
		double x, x_1, err;
    	int i=0;
    	cout << "M\202todo de Newton-Raphson, C\240lculo de la funci\242n: cos(x)*x^2-2*x" << endl << endl;
    	cout << "Ingrese el valor inicial de x0: ";
    	cin >> x;
    	cout << endl;
    	do {
        	x_1 = x;
        	x = x - funcion (x) / derivada (x);
        	err = fabs ((x - x_1) / x);
        	cout << "x" << i << " = " << x << "\t\terror = " << err << endl;
        	i++;
    	} while (x!=x_1 && i<100);
    		if (i==100){
        		cout << endl << "La soluci\242n no es convergente." << endl << endl;
    	}else{
        	cout << endl << "La soluci\242n es " << x << endl << endl;
		}	

   		
   
   cout << "M\202todo de Newton-Raphson, C\240lculo de la funci\242n: ((6-2)/x^2)/(4+1)" << endl << endl;
    	cout << "Ingrese el valor inicial de x0: ";
    	cin >> x;
    	cout << endl;
    	do {
        	x_1 = x;
        	x = x - funcion2 (x) / derivada2 (x);
        	err = fabs ((x - x_1) / x);
        	cout << "x" << i << " = " << x << "\t\terror = " << err << endl;
        	i++;
    	} while (x!=x_1 && i<100);
    		if (i==100){
        		cout << endl << "La soluci\242n no es convergente." << endl << endl;
    	}else{
        	cout << endl << "La soluci\242n es " << x << endl << endl;
		}
		
		
		
	 cout << "M\202todo de Newton-Raphson, C\240lculo de la funci\242n: x^3 - 3*sin(x^2)+1" << endl << endl;
    	cout << "Ingrese el valor inicial de x0: ";
    	cin >> x;
    	cout << endl;
    	do {
        	x_1 = x;
        	x = x - funcion3 (x) / derivada3 (x);
        	err = fabs ((x - x_1) / x);
        	cout << "x" << i << " = " << x << "\t\terror = " << err << endl;
        	i++;
    	} while (x!=x_1 && i<100);
    		if (i==100){
        		cout << endl << "La soluci\242n no es convergente." << endl << endl;
    	}else{
        	cout << endl << "La soluci\242n es " << x << endl << endl;
		}
		
		cout << "M\202todo de Newton-Raphson, C\240lculo de la funci\242n: x**3 + 6*x**2 + 9.4*x + 2.5" << endl << endl;
    	cout << "Ingrese el valor inicial de x0: ";
    	cin >> x;
    	cout << endl;
    	do {
        	x_1 = x;
        	x = x - funcion3 (x) / derivada3 (x);
        	err = fabs ((x - x_1) / x);
        	cout << "x" << i << " = " << x << "\t\terror = " << err << endl;
        	i++;
    	} while (x!=x_1 && i<100);
    		if (i==100){
        		cout << endl << "La soluci\242n no es convergente." << endl << endl;
    	}else{
        	cout << endl << "La soluci\242n es " << x << endl << endl;
		}
		
	cout <<"-------------------Ahora se ejcutara el metodo de biseccion---------------"<<endl;
	float xi, xu, xr, xa, ea;
	int iteraciones;
	cout << "M\202todo de biseccion, C\240lculo de la funci\242n: cos(x)*x^2-2*x" << endl << endl;
	printf("Limite superior de la raiz: ");
	scanf("%f", &xi);
	printf("Limite inferior de la raiz: ");
	scanf("%f", &xu);
	
	iteraciones = 0;
	xa = 0.0;
	ea = 0;
	printf("%12s %10s %10s %10s %10s\n", "Iteraciones", "Xi", "Xu", "Xr", "Ea");
 
	do{
		iteraciones++;
		xr = (xi + xu) / 2;
		if(iteraciones >= 2) {
			ea = ((xr - xa) / xr) * 100; // Calculamos el error aproximado a partir de la segunda iteración
		}
		ea = fabs(ea); // Calculamos el valor absoluto del número
		printf("%12d %10f %10f %10f %10f\n", iteraciones, xi, xu, xr, fabs(ea));
		xa = xr;
		if(funcion(xr) > 0) {
			xi = xr;
		} else{
			xu = xr;
		}
 
	} while(ea > 0.5 || (ea == 0.0 && iteraciones == 1));
 
	printf("\nLa raiz es de f(x) = cos(x)*x^2-2*x; es: %f\n", xr);
	
	cout << "M\202todo de biseccion, C\240lculo de la funci\242n: ((6-2)/x^2)/(4+1)" << endl << endl;
	printf("Limite superior de la raiz: ");
	scanf("%f", &xi);
	printf("Limite inferior de la raiz: ");
	scanf("%f", &xu);
	
	iteraciones = 0;
	xa = 0.0;
	ea = 0;
	printf("%12s %10s %10s %10s %10s\n", "Iteraciones", "Xi", "Xu", "Xr", "Ea");
 
	do{
		iteraciones++;
		xr = (xi + xu) / 2;
		if(iteraciones >= 2) {
			ea = ((xr - xa) / xr) * 100; // Calculamos el error aproximado a partir de la segunda iteración
		}
		ea = fabs(ea); // Calculamos el valor absoluto del número
		printf("%12d %10f %10f %10f %10f\n", iteraciones, xi, xu, xr, fabs(ea));
		xa = xr;
		if(funcion2(xr) > 0) {
			xi = xr;
		} else{
			xu = xr;
		}
 
	} while(ea > 0.5 || (ea == 0.0 && iteraciones == 1));
 
	printf("\nLa raiz es de f(x) = ((6-2)/x^2)/(4+1); es: %f\n", xr);
	
	cout << "M\202todo de biseccion, C\240lculo de la funci\242n: x^3 - 3*sin(x^2)+1" << endl << endl;
	
	printf("Limite superior de la raiz: ");
	scanf("%f", &xi);
	printf("Limite inferior de la raiz: ");
	scanf("%f", &xu);
	
	iteraciones = 0;
	xa = 0.0;
	ea = 0;
	printf("%12s %10s %10s %10s %10s\n", "Iteraciones", "Xi", "Xu", "Xr", "Ea");
 
	do{
		iteraciones++;
		xr = (xi + xu) / 2;
		if(iteraciones >= 2) {
			ea = ((xr - xa) / xr) * 100; // Calculamos el error aproximado a partir de la segunda iteración
		}
		ea = fabs(ea); // Calculamos el valor absoluto del número
		printf("%12d %10f %10f %10f %10f\n", iteraciones, xi, xu, xr, fabs(ea));
		xa = xr;
		if(funcion3(xr) > 0) {
			xi = xr;
		} else{
			xu = xr;
		}
 
	} while(ea > 0.5 || (ea == 0.0 && iteraciones == 1));
 
	printf("\nLa raiz es de f(x) = x^3 - 3*sin(x^2)+1; es: %f\n", xr);
	cout << "M\202todo de biseccion, C\240lculo de la funci\242n: x**3 + 6*x**2 + 9.4*x + 2.5" << endl << endl;
	printf("Limite superior de la raiz: ");
	scanf("%f", &xi);
	printf("Limite inferior de la raiz: ");
	scanf("%f", &xu);
	
	iteraciones = 0;
	xa = 0.0;
	ea = 0;
	printf("%12s %10s %10s %10s %10s\n", "Iteraciones", "Xi", "Xu", "Xr", "Ea");
 
	do{
		iteraciones++;
		xr = (xi + xu) / 2;
		if(iteraciones >= 2) {
			ea = ((xr - xa) / xr) * 100; // Calculamos el error aproximado a partir de la segunda iteración
		}
		ea = fabs(ea); // Calculamos el valor absoluto del número
		printf("%12d %10f %10f %10f %10f\n", iteraciones, xi, xu, xr, fabs(ea));
		xa = xr;
		if(funcion4(xr) > 0) {
			xi = xr;
		} else{
			xu = xr;
		}
 
	} while(ea > 0.5 || (ea == 0.0 && iteraciones == 1));
 
	printf("\nLa raiz es de f(x) = x**3 + 6*x**2 + 9.4*x + 2.5; es: %f\n", xr);	
	

	};
	
	if(p == 2){
		cout <<"-------------------Bienvenido al progama 2---------------"<<endl;
		int C,Dimension;
 		float Matriz[Tam][Tam];
 		float Sistema[101][102];
 		PideDatos(&Dimension,Matriz);
		printf("\n\n\nCalcula DETERMINANTE: \n\n");
 		
		EscribeDatos(Dimension,Matriz);
		CalcDet(Dimension,Matriz);
		PideDatos2(&Dimension,Sistema);
		ResuelveGauss(Dimension,Sistema);
    	printf("\n\n\nLas soluciones son:\n");
    	for(C=1;C<=Dimension;C++) printf("\n X%d=%f\n",C,Sistema[C][Dimension+1]);
    	
 		scanf("%d");
	};
	
	
	if( p == 3){
    int i,j,iteraciones=0;
    float error,epsilon;
    printf("\n METODO DE JACOBI DE RESOLUCION DE SISTEMAS Ax=b \n");

    printf("Dimension de la matriz A: ");
    scanf("%d",&dim);
    float A[dim][dim],b[dim],x[dim],x_prev[dim],aux[dim];

    printf("\n Elementos de la matriz A: \n");
    for(i=0;i<dim;i++) for(j=0;j<dim;j++){
        printf("A(%d,%d)=",i,j); scanf("%f",&A[i][j]);
    }

    printf("\n Elementos del vector b: \n");
    for(i=0;i<dim;i++){
        printf("b(%d)=",i); scanf("%f",&b[i]);
    }

    printf("\n Error de parada: \n");
    printf("E=",i); scanf("%f",&epsilon);
    error=epsilon+1;

    
    printf("\n Valor inicail de la iteracion: \n");
    for(i=0;i<dim;i++){
        printf("x0(%d)=",i); scanf("%f",&x_prev[i]);
    }
    while (error>epsilon){
        for(i=0;i<dim;i++){
            for(j=0;j<dim;j++) aux[j]=A[i][j];
            x[i]=(1/A[i][i])*(b[i]-suma_jacobi(aux,x_prev,i));
        }
        error=norma(x,x_prev);

        printf("\n\n Iteracion %d: \n",iteraciones);
        for(i=0;i<dim;i++){
            x_prev[i]=x[i];
            printf("X(%d)=%f \n",i,x[i]);
        }

        iteraciones++;
        if (iteraciones==10) error=epsilon-1;
    }

    printf("Solucion del sistema\n");
    printf("Numero de iteraciones: %d \n", iteraciones);
    for(i=0;i<dim;i++){
        printf("x(%d)=%f\n",i,x[i]);
    }
    return 1;	
	}
	
	if(p ==4){
	printf("Este programa aplica el metodo de las potencias .\n"); 
	printf("\n");


	printf("Introduzca la dimension de la matriz. La matriz tiene que ser cuadrada por definicion: \n");
	int n=0,i=0,j=0;

	scanf("%d", &n); 
	float coeficientes[n][n], inicial[n], vector_iteraciones[n], resultado;

	printf("\n");

	printf("Introduzca la matriz que desea: \n");
	for(i=0; i<n;i++){
		for(j=0; j<n; j++){
			printf("A[%d][%d]=",i,j);
			scanf("%f", &coeficientes[i][j]);
		}
	}


	printf("Su matriz queda de la siguiente manera: \n");	
	for(i=0; i<n; i++){
	 	printf("|");
	 	for(j=0; j<n; j++){
	 		printf(" %.2f ", coeficientes[i][j]);
	 	}
	 	printf(" | \n");
	}

	printf("\n");	

	
	printf("Introduzca la condicion inicial en orden (x,y,...): \n");
	for(i=0; i<n;i++){
		printf("X[%d]=",i);
		scanf("%f", &inicial[i]);
	}

	printf("Su matriz queda de la siguiente manera: \n");	
	for(i=0; i<n; i++){
	 	printf("|");
	 	printf(" %.2f ", inicial[i]);
	 	printf(" | \n");
	}

	printf("Procedemos a aplicar el metodo....\n");

	back1 :


	for(i=0; i<n;i++){
		vector_iteraciones[i]=0; 
		for(j=0; j<n; j++){
			vector_iteraciones[i]=vector_iteraciones[i]+coeficientes[i][j]*inicial[j]; 
       	}                                                                             																	
 	}                                                                                

	if(inicial[0]!=0){ 
		resultado=vector_iteraciones[0]/inicial[0];
	}
	printf("El resultado de la iteracion: %f \n", resultado);

	for(i=0; i<n; i++){
		inicial[i]=vector_iteraciones[i]; 
	}

	
	}
	return 0;
}





float norma(float vector1[],float vector2[]){
    float aux=0;
    int i;
    for(i=0;i<dim;i++){
        aux=aux+(vector1[i]-vector2[i])*(vector1[i]-vector2[i]);
    }
    return aux;
}

float suma_jacobi(float Matriz[], float vector[], int componente){
    float aux=0;
    int i;
    for(i=0;i<dim;i++){
        if (componente!=i){
            aux=aux+Matriz[i]*vector[i];
        }
    }
    return aux;
}



void PideDatos(int *Dim,float Mat[][Tam]){
 int A,B;
 printf("\n\n ||CALCULA DETERMINANTE ESCALONANDO MATRIZ||");
 printf("\n\n\n Introduce la dimension de la matriz:");
 scanf("%d",&*Dim);
 printf("\n\n INTRODUCIR CADA COMPONENTE DE LA MATRIZ:");
 for(A=1;A<=*Dim;A++) for(B=1;B<=*Dim;B++){
  printf("\n Termino A(%d,%d):",A,B); scanf("%f",&Mat[A][B]);}
}

void PideDatos2(int *Dim,float Sist[][102]){
    int A,B;
    printf("\n\n ||RESUELVE SISTEMAS LINEALES DETERMINADOS POR GAUSS||");
    printf("\n\n\n Introduce el numero de incognitas:(menor que 100)");
    scanf("%d",&*Dim);
    printf("\n\n PASE A INTRODUCIR CADA COMPONENTE DEL SISTEMA (A|B):");
    printf("\n\n MATRIZ A:\n");
    for(A=1;A<=*Dim;A++) for(B=1;B<=*Dim;B++){
        printf("\n Termino A(%d,%d):",A,B); scanf("%f",&Sist[A][B]);}
    printf("\n\n\n VECTOR B:\n");
    for(A=1;A<=*Dim;A++){
        printf("\n Termino B(%d):",A);scanf("%f",&Sist[A][*Dim+1]);
    }}

void EscribeDatos(int Dim, float Mat[][Tam]){
 int A,B;
 for(A=1;A<=Dim;A++){
  for(B=1;B<=(Dim);B++)
   printf("%7.2f",Mat[A][B]);
  printf("\n");
 }}

void CalcDet(int Dim, float Mat[][Tam]){
 int NoCero,Col,C1,C2,A,NoReg,Perm=0/*permutaciones*/;
 float Pivote,V1,Det=1;

 for(Col=1;Col<=Dim;Col++){
  NoCero=0;A=Col;
  while((NoCero==0)&&(A<=Dim)){
   if((Mat[A][Col]>0.0000001)||((Mat[A][Col]<-0.0000001))){
    NoCero=1;}
   else A++;}
  if (A>Dim) NoReg=1;
  if (A!=Col) Perm++;
  Pivote=Mat[A][Col];
  for(C1=1;C1<=(Dim);C1++){
   V1=Mat[A][C1];
   Mat[A][C1]=Mat[Col][C1];
   Mat[Col][C1]=V1;}
  for(C2=Col+1;C2<=Dim;C2++){
   V1=Mat[C2][Col];
   for(C1=Col;C1<=(Dim);C1++){
    Mat[C2][C1]=Mat[C2][C1]-((V1/Pivote)*Mat[Col][C1]);}}
 }
    for(C2=1;C2<=Dim;C2++) Det=Det*Mat[C2][C2];
 A=Perm;
 if ((A%2)==1) Det=-Det; /*Caso de permutaciones impares*/
 if (NoReg==1) Det=0;
    printf("El determinante de la matriz es:   %f", Det);

}

void ResuelveGauss(int Dim, float Sist[][102]){
    int NoCero,Col,C1,C2,A;
    float Pivote,V1;
   
    for(Col=1;Col<=Dim;Col++){
        NoCero=0;A=Col;
        while(NoCero==0){
            if(Sist[A][Col]!=0){
                NoCero=1;}
            else A++;}
        Pivote=Sist[A][Col];
        for(C1=1;C1<=(Dim+1);C1++){
            V1=Sist[A][C1];
            Sist[A][C1]=Sist[Col][C1];
            Sist[Col][C1]=V1/Pivote;}
        for(C2=Col+1;C2<=Dim;C2++){
            V1=Sist[C2][Col];
            for(C1=Col;C1<=(Dim+1);C1++){
                Sist[C2][C1]=Sist[C2][C1]-V1*Sist[Col][C1];}
    }}
   
    for(Col=Dim;Col>=1;Col--) for(C1=(Col-1);C1>=1;C1--){
        Sist[C1][Dim+1]=Sist[C1][Dim+1]-Sist[C1][Col]*Sist[Col][Dim+1];
        Sist[C1][Col]=0;
    }
}

