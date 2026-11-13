#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define EULERO 1
#define VELOCITY_VERLET 3
#define RK2 4
#define RK4 5
#define M 1.0
#define DT_MAX 0.064 /* il dt piu' grande */
#define N_DT 8      /* numero di valori di dt: 0.064, 0.032, ..., 0.064 / 2^7 = 0.0005 */

/* Lo stato del sistema, con i campi x (posizione), v (velocita') e
 * a (accelerazione). Il campo a e' l'accelerazione calcolata nello stato
 * (x, v): la si conserva insieme a x e v cosi' che velocity Verlet non debba
 * ricalcolarla all'inizio di ogni passo. Eulero, RK2 e RK4 possono ignorarlo.
 * La stessa struttura si usa anche per le derivate restituite da f
 * (dx/dt nel campo x e dv/dt nel campo v; il campo a non si usa).
 */
typedef struct {
  
  double x;
  double v;
  double a;
  
} Punto;

/* I parametri della simulazione e del sistema fisico: almeno il passo di
 * integrazione dt e la pulsazione omega (con m = 1 non serve altro). Per il
 * pendolo della seconda parte potete usare lo stesso campo omega, con
 * omega = sqrt(g / L), oppure aggiungere i campi che vi servono. */
typedef struct {
  
  double dt;
  double w;
  
} Sistema;

/* Accelerazione e' il nome di un tipo: "puntatore a una funzione che riceve
 * (double t, Punto s, const Sistema *sist) e restituisce un Punto".
 * Una variabile di questo tipo puo' contenere l'indirizzo di funzioni come
 * oscillatore (vedi sotto): gli integratori la ricevono come argomento f
 * e la chiamano con f(t, s, sist) senza sapere quale sistema stanno integrando. */
typedef Punto (*Accelerazione)(double t, Punto s, const Sistema *sist);

/* Puntatore a un integratore: una funzione che avanza lo stato *s di un passo
 * sist->dt, a partire dal tempo t, usando solo la funzione f per conoscere le
 * derivate. Tutti gli integratori qui sotto hanno questa forma. */
typedef void (*Integratore)(double t, Punto *s, const Sistema *sist, Accelerazione f);

/* La funzione da integrare per l'oscillatore armonico. Ha la forma
 * richiesta dal tipo Accelerazione e restituisce dx/dt = v nel campo x e
 * dv/dt = a = -omega^2 x nel campo v. Sostituite il valore provvisorio. */

double trueX(double tf, double x0, double v0, double w0) {
  double trueX = x0 * cos(w0 * tf) + (v0 / w0) *  sin(w0 * tf);
  return trueX;
}

Punto oscillatore(double t, Punto s, const Sistema *sist) {
  
  Punto d;

  d.x = s.v;
  d.v = -sist->w * sist->w * s.x;
  d.a = 0;
  
  return d;
}

Punto vanderpol(double t, Punto s, const Sistema *sist) {

  Punto d;
  
  d.x = s.v;
  d.v = sist->w * (1 - s.x * s.x) * s.v - s.x; // uso w come mu
  return d;
}

/* Energia meccanica, E = m v^2 / 2 + m omega^2 x^2 / 2.
 * Potete riusare la funzione dell'esercitazione 3. */
double energia(Punto s, const Sistema *sist) {
  double E = 0.5 * (s.x * s.x + s.v * s.v);
  return E;
}

/* Esempio: un passo del metodo di Eulero scritto nella nuova forma.
 * L'integratore non sa quale sistema sta integrando: chiede le derivate a f. */
void eulero(double t, Punto *s, const Sistema *sist, Accelerazione f) {
    Punto d = f(t, *s, sist);
    s->x = s->x + d.x * sist->dt;
    s->v = s->v + d.v * sist->dt;
}

/* Un passo di velocity Verlet nella nuova forma. Partite dal vostro
 * codice dell'esercitazione 3.
 * All'ingresso s->a contiene GIA' l'accelerazione al tempo t (calcolata al
 * passo precedente, oppure in main prima del ciclo): non va ricalcolata.
 * Serve una sola chiamata a f, nella nuova posizione: l'accelerazione e' nel
 * campo v di cio' che f restituisce. Prima di uscire salvatela in s->a, cosi'
 * il passo successivo la trova pronta. */
void velocity_verlet(double t, Punto *s, const Sistema *sist, Accelerazione f) {
  
  double vm = s->v + 0.5 * s->a * sist->dt;
  s->x = s->x + vm * sist->dt;
  s->a = f(t, *s, sist).v;
  s->v = vm + s->a * 0.5 * sist->dt;
  
}

/* Un passo di RK2 (punto medio). Non deve contenere nulla di specifico dell'oscillatore: solo chiamate a f. */
void rk2(double t, Punto *s, const Sistema *sist, Accelerazione f) {
  Punto k1 = f(t, *s, sist);
  Punto sm;
  sm.x = s->x + 0.5 * sist->dt * k1.x;
  sm.v = s->v + 0.5 * sist->dt * k1.v;
  
  Punto k2 = f(t + sist->dt * 0.5, sm, sist);
  s->x = s->x + sist->dt * k2.x;
  s->v = s->v + sist->dt * k2.v;

  
}

/* Un passo di RK4. */
void rk4(double t, Punto *s, const Sistema *sist, Accelerazione f) {

 Punto k1 = f(t, *s, sist);
  Punto sm1;
  sm1.x = s->x + 0.5 * sist->dt * k1.x;
  sm1.v = s->v + 0.5 * sist->dt * k1.v;
  
  Punto k2 = f(t + sist->dt * 0.5, sm1, sist);
  Punto sm2;
  sm2.x = s->x + 0.5 * sist->dt * k2.x;
  sm2.v = s->v + 0.5 * sist->dt * k2.v;

  Punto k3 = f(t + sist->dt * 0.5, sm2, sist);
  Punto sm3;
  sm3.x = s->x + sist->dt * k3.x;
  sm3.v = s->v + sist->dt * k3.v;

  Punto k4 = f(t + sist->dt, sm3, sist);

  s->x = s->x + sist->dt / 6 * (k1.x + 2 * k2.x + 2 * k3.x + k4.x);
  s->v = s->v + sist->dt / 6 * (k1.v + 2 * k2.v + 2 * k3.v + k4.v);

}



int main(int argc, char *argv[]) {
    /* Come nell'esercitazione 3 il programma fa un ciclo su dt e stampa due
     * colonne, dt ed eps_G. Al posto del numero di passi si passa il tempo
     * finale tf, lo stesso per tutti i dt:
     *     ./oscillatore metodo tf x0 v0 omega0 */
    if(argc != 7) {
        fprintf(stderr, "Uso: %s metodo tf x0 v0 omega0\n", argv[0]);
        return 2;
    }
 
    double dt = atof(argv[1]);
    int nsteps = atoi(argv[2]);
    int stride = atoi(argv[3]);
    double x0 = atof(argv[4]);
    double v0 = atof(argv[5]);
    double w = atof(argv[6]);
    Sistema sist;
    const Sistema *p_sist = &sist;
    sist.dt = dt;
    sist.w = w;
    double t = 0;
    Punto p;
    Punto *q = &p;
    p.x = x0;
    p.v = v0;
    p.a = 0;
    Accelerazione p_vanderpol = vanderpol;
    int i = 0;
    int j = 0;

    double peak_x = 0;
    double lowahh_x = 0;
    double peak_t = 0;
    double lowahh_t = 0;
    double current_v = 0;
    double prev_v = 0;
    
    /* si sceglie l'integratore una volta sola, prima del ciclo, assegnando al puntatore a funzione `passo` la funzione giusta. */
    Integratore passo;
    passo = rk4;
    
    FILE* filer = fopen("xt6.dat", "w");
    
    for(i = 0; i < nsteps; i++) {
      
      passo(t, q, p_sist, p_vanderpol);

      prev_v = current_v;
      current_v = p.v;

      if ((prev_v * current_v) < 0) {
	
	lowahh_x = peak_x;
	peak_x = p.x;

	lowahh_t = peak_t;
	peak_t = t;
	
      }
      
      if (i % (stride) == 0) {
	printf("%lf %lf %lf %lf\n", t, p.x, p.v, energia(p, p_sist));
	fprintf(filer, "%lf %lf %lf %lf\n", t, p.x, p.v, energia(p, p_sist));
      }
      
    t = t + sist.dt;
    
    }

    fclose(filer);

    double Tf = fabs(peak_t - lowahh_t);
    double Af = fabs(peak_x - lowahh_x);
    
    printf("il periodo finale è %lf, l'ampiezza finale è %lf", Tf, Af);


    double mui = sist.w;
    FILE* filer1 = fopen("mu_Tf_Af6.dat", "w");
    for(j = 0; j < 10; j++) {

      t = 0;
      
      for(i = 0; i < nsteps; i++) {
	
	passo(t, q, p_sist, p_vanderpol);
	
	prev_v = current_v;
	current_v = p.v;
	
	if ((prev_v * current_v) < 0) {
	  
	  lowahh_x = peak_x;
	  peak_x = p.x;
	  
	  lowahh_t = peak_t;
	  peak_t = t;
	  
	}

	t = t + sist.dt;
	
      }
      	
	Tf = fabs(peak_t - lowahh_t);
	Af = fabs(peak_x - lowahh_x);
	printf("%lf %lf %lf\n", p_sist->w, Tf, Af);
	fprintf(filer1, "%lf %lf %lf\n", p_sist->w, Tf, Af);
	
	sist.w = sist.w + mui * 0.5;
    }
    
    fclose(filer1);
    
    return 0;
}
