#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//si assume massa unitaria

typedef struct{

  double x;
  double v;
  double a;

} Punto;


typedef struct{

  double dt;
  double w;

} Sistema;


Punto oscillatore(Punto* p, Sistema* sist) {

  Punto d;
  d.x = p->v;
  d.v = -sist->w * sist->w * p->x;
  return d;
  
}

double E_oscillatore(Punto* p, Sistema* sist) {
  
  double E = 0;
  E = 0.5 * p->v * p->v + 0.5 * sist->w * sist->w * p->x * p->x;
  return E;
  
}


typedef Punto (*Derivatore)(Punto*, Sistema*); //ti da le derivate


typedef void (*Integratore)(Punto*, Sistema*, Derivatore); //fa un passo


typedef double (*Energia)(Punto*, Sistema*); //comodo, opzionale


void eulero(Punto* p, Sistema* sist, Derivatore Der) {
  
  Punto d = Der(p, sist);
  p->x = p->x + d.x * sist->dt;
  p->v = p->v + d.v * sist->dt;
  
}

void eulero_cromer(Punto* p, Sistema* sist, Derivatore Der) {
  
  Punto d = Der(p, sist);
  p->v = p->v + d.v * sist->dt;
  p->x = p->x + p->v * sist->dt;
  
}



int main(int argc, char* argv[]) {

  
  if (argc != 8) {
    
    printf("Il programma usa la sintassi ./oscillatore.exe [metodo] [dt] [nsteps] [stride] [x0] [v0] [omega]\n");
    exit(1);
    
  }

  Sistema sist;
  Sistema* p_sist = &sist;

  Punto p;
  Punto* p_p = &p;
  
  int metodo = atoi(argv[1]);
  sist.dt = atof(argv[2]);
  int nsteps = atoi(argv[3]);
  int stride = atoi(argv[4]);
  double x0 = atof(argv[5]);
  double v0 = atof(argv[6]);
  p.x = atof(argv[5]);
  p.v = atof(argv[6]);
  double omega0 = atof(argv[7]);
  sist.w = atof(argv[7]);

  Derivatore der_corrente = oscillatore;
  Integratore passo;
  Energia E_corrente = E_oscillatore; //cambia qui l'energia
  
  if(metodo == 1) {
    
    passo = eulero;
    
  } else {
    
    if(metodo == 2) {
      
      passo = eulero_cromer;
      
    }else {printf("metodo invalido\n");}
    
  }

  int i = 0;
  double t = 0;

  FILE *filer = fopen("dati.dat", "w");
  
  for(i = 0; i < nsteps; i++) {

    passo(p_p, p_sist, der_corrente);
    t = t + sist.dt;

    if(((i + 1) % stride) == 0) {
      printf("%lf %lf %lf %lf\n", t, p.x, p.v, E_corrente(p_p, p_sist));
      fprintf(filer, "%lf %lf %lf %lf\n", t, p.x, p.v, E_corrente(p_p, p_sist));
    }
    
  }

  fclose(filer);
  
  return 0;
}
