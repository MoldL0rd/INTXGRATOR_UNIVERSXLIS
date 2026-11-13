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

Punto pendolo(Punto* p, Sistema* sist) { //x come angolo e w come lunghezza
  
  Punto d;
  d.x = p->v;
  d.v = -( 9.81 / sist->w ) * sin(p->x);
  return d;
  
}

double E_oscillatore(Punto* p, Sistema* sist) {
  
  double E = 0;
  E = 0.5 * p->v * p->v + 0.5 * sist->w * sist->w * p->x * p->x;
  return E;
  
}

double E_pendolo(Punto* p, Sistema* sist) { //w è la lunghezza
  
  double E = 0;
  E = 0.5 * sist->w * sist->w * p->v * p->v + 9.81 * sist->w * (1 - cos(p->x));
  return E;
  
}

double errore_globale_oscillatore(double Tf, double x_f, double x0, double v0, double w0) {

  double E_g = fabs( (x_f) - (x0 * cos(w0 * Tf) + (v0 / w0) * sin(w0 * Tf)) ) / sqrt( (x0 * x0) + (v0/w0) * (v0/w0) );
  return E_g;
  
}

typedef Punto (*Derivatore)(Punto*, Sistema*); //ti da le derivate


typedef void (*Integratore)(Punto*, Sistema*, Derivatore); //fa un passo


typedef double (*Energia)(Punto*, Sistema*); //comodo, opzionale


typedef double (*Errore_globale)(double, double, double, double, double); //comodo, opzionale


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

void velocity_verlet(Punto* p, Sistema* sist, Derivatore Der) {

  double vm = 0;
  Punto d = Der(p, sist);
  vm = d.x + 0.5 * d.v * sist->dt;
  p->x = p->x + vm * sist->dt;
  p->a = Der(p, sist).v;
  p->v = vm + 0.5 * p->a * sist->dt;
  
}

void RK2(Punto* p, Sistema* sist, Derivatore Der) {

  Punto k1 = Der(p, sist);
  Punto m1;
  m1.x = p->x + k1.x * 0.5 * sist->dt;
  m1.v = p->v + k1.v * 0.5 * sist->dt;

  Punto k2 = Der(&m1, sist);
  p->x = p->x + k2.x * sist->dt;
  p->v = p->v + k2.v * sist->dt;
  
}


void RK4(Punto * p, Sistema* sist, Derivatore Der) {

  Punto k1 = Der(p, sist);
  
  Punto m1;
  m1.x = p->x + k1.x * 0.5 * sist->dt;
  m1.v = p->v + k1.v * 0.5 * sist->dt;
  Punto k2 = Der(&m1, sist);
  
  Punto m2;
  m2.x = p->x + 0.5 * k2.x * sist->dt;
  m2.v = p->v + 0.5 * k2.v * sist->dt;
  Punto k3 = Der(&m2, sist);

  Punto m3;
  m3.x = p->x + k3.x * sist->dt;
  m3.v = p->v + k3.v * sist->dt;
  Punto k4 = Der(&m3, sist);

  p->x = p->x + (sist->dt / 6.0) * (k1.x + 2 * k2.x + 2 * k3.x + k4.x);
  p->v = p->v + (sist->dt / 6.0) * (k1.v + 2 * k2.v + 2 * k3.v + k4.v);
  
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
  double w0 = atof(argv[7]);
  sist.w = atof(argv[7]);


  //pannello di controllo per cambiare il sistema
  Derivatore der_corrente = oscillatore;
  Integratore passo;
  Energia E_corrente = E_oscillatore; //cambia qui l'energia
  Errore_globale err_corrente = errore_globale_oscillatore;
  
  if(metodo == 1) {
    
    passo = eulero;
    
  } else {
    
    if(metodo == 2) {
      
      passo = eulero_cromer;
      
    } else {
      
      if(metodo == 3) {
	
	passo = velocity_verlet;
	
      } else {
	
	if(metodo == 4) {
	  
	  passo = RK2;
	  
	} else {
	  
	  if(metodo == 5) {
	    
	    passo = RK4;
	    
	  } else {
	    
	    printf("metodo invalido\n");
	    
	  }
	  
	}
	
      }
      
    }
    
  }
  
  printf("Integrazione principale:\n\n");
  
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
  
  double xf = p.x; 
  double Tf = t;
  
  printf("L'errore globale è di %lf\n\n\n", err_corrente(Tf, xf, x0, v0, w0));



  printf("Ciclo sui diversi delta t:\n\n");
  
  t = 0;
  int j = 0;

  FILE* filer1 = fopen("dt_errore.dat", "w");
  
  for(j = 0; j < 7; j++) {

    p.x = x0;
    p.v = v0;
    
    for(i = 0; i < nsteps; i++) {
      
      passo(p_p, p_sist, der_corrente);
      t = t + sist.dt;
      
      if(((i + 1) % (stride * 2)) == 0) {
	printf("%lf %lf %lf %lf\n", t, p.x, p.v, E_corrente(p_p, p_sist));
      }
      
    }
    
    printf("nsteps dell'integrazione appena fatta = %d\n", nsteps);
    
    Tf = t;
    xf = p.x;
    printf("l'errore per dt = %lf è: %lf\n\n", sist.dt, err_corrente(Tf, xf, x0, v0, w0));
    fprintf(filer1, "%lf %lf\n", sist.dt, err_corrente(Tf, xf, x0, v0, w0));
    sist.dt = sist.dt * 2;
    nsteps = Tf / sist.dt;
    t = 0;
 
  }

  fclose(filer1);
  
  return 0;
}
