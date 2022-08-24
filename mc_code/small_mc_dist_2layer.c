char   t1[80] = "Small Monte Carlo by Scott Prahl (http://omlc.ogi.edu)";
char   t2[80] = "modified to keep track of distance and for two layers by Ben Smith";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define BINS 101
int DEBUG=0;
int ASCII=1; 
double g = 0.0;	/* Scattering Anisotropy -1<=g<=1 */
double z0=0; /* start photons at this level */
long this_layer;
long   i, photons = 10000000;
double d,x,y,z,u,v,w,weight,g2,z_max;
double t_layer;
double dist[2], z_top[2], z_bot[2];

char *outfile;
FILE *outPtr;
double outbuf[11];

double rs, rd;
long count;
long maxrand;
double d;
long N_max=1000;
double layer_number=0;

void print_help() {
    fprintf(stderr,"small_mc -g asymmetry_parameter -p n_Photons -z z0 -l layer_thickness -N N_max\n");
    fprintf(stderr,"for binary output: x, y, z, u, v, w, dist0, dist1, z_max, event count, photon\n");
}

void dump() {
  if (ASCII==1) {
    printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %10.6f %10.6f %10.6f %10ld %16ld \n", x, y, z, u, v, w, dist[0], dist[1], z_max, count, i);
  } else {
  outbuf[0]=x;
  outbuf[1]=y;
  outbuf[2]=z;
  outbuf[3]=u;
  outbuf[4]=v;
  outbuf[5]=w;
  outbuf[6]=dist[0];
  outbuf[7]=dist[1];
  outbuf[8]=z_max;
  outbuf[9]=(double) i;
  outbuf[10]=(double) count ;
  fwrite(outbuf, sizeof(double), 11, outPtr);
  }
}

void launch() /* Start the photon */
{
  if(DEBUG==1) {  printf("launch\n"); }
  x = 0.0; y = 0.0; z = z0;		  
  u = 0.0; v = 0.0; w = -1.0;
  z_max=0;
  count=0;
  dist[0]=0.;
  dist[1]=0.;
  this_layer=0;
}

void move() /* move to next scattering or absorption event */
{
  double next_z;  /*next layer boundary that might be encountered on the current path */
  double dist_this;
  /*d = -log((random()+1.0)/(maxrand+1.0)); */
  if (w>0) {
    next_z = z_top[this_layer];
  }else{
    next_z = z_bot[this_layer];
  }
  dist_this=(next_z-z)/w;  /*distance to the next layer boundary */
  if (dist_this > d) {
    /* trajectory ends in the current layer */  
    dist[this_layer]+=d;
    x += d * u;
    y += d * v;
    z += d * w; 
    d=0;
    if (z_max > z) {
      z_max=z;
    }
    if (DEBUG==1){
      printf("stayed in layer %ld\n", this_layer);
      printf("\t--%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %10.6f %10.6f %10ld %16ld \n", x, y, z, u, v, w, dist[0], dist[1], count, i);
    }
    return;  
  }
  if (DEBUG==1) { 
    printf("\tbefore move, %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %10.6f %10.6f %10ld %16ld \n", x, y, z, u, v, w, dist[0], dist[1], count, i);
    printf("\t--exiting layer %ld, w=%f, next_z=%f\n", this_layer, w, next_z);
    printf("\t--z_top=%f, z_bot=%f\n", z_top[this_layer], z_bot[this_layer]);
  }
  /* OTW record the distance in this layer, continue on to the next layer */
  dist[this_layer]+=dist_this;
  x += dist_this * u;
  y += dist_this * v;
  z += dist_this * w;
  d -= dist_this;
  if (DEBUG==1) {
    printf("\tafter move, %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %10.6f %10.6f %10ld %16ld \n", x, y, z, u, v, w, dist[0], dist[1], count, i);
  }
  if (z < z_max){
    z_max=z;
  }
  if (this_layer<=0 && w>0) {
    /* exit location at the top*/
    dump();
    count=2e9+1;
    d=0;
    return;
  } 
  if (this_layer>=1 && w<0) { 
    /* exit location on the bottom */
    dump();
    count=2e9+1;
    d=0;
    return;
  }
  /* from this point, we are crossing into the next layer */
  /* increment the layer number */
  if (w>0) { this_layer -=1;  }
  if (w<=0) { this_layer +=1; }

  if ((this_layer > 1) || (this_layer < 0)) { printf("layer=%16ld count=%16ld, i=%16ld\n", this_layer, count, i); }
  /* escape out the bottom 
  if (w<0 && this_layer ==1) {
    printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %10.6f %10.6f %10ld %16ld \n", x, y, z, u, v, w, dist[0], dist[1], count, i);
    count=2e9+1;
  }
  */
}

void scatter() /* Scatter photon and establish new direction */
{ 
  double x1, x2, x3, t, mu, mu2, g2, temp1; 

  for(;;) {								/*new direction*/
    x1=2.0*random()/maxrand - 1.0; 
    x2=2.0*random()/maxrand - 1.0; 
    if ((x3=x1*x1+x2*x2)<1.0) break;
  }	
  if (g==0) {  /* isotropic */
    u = 2.0 * x3 -1.0;
    v = x1 * sqrt((1.0-u*u)/x3);
    w = x2 * sqrt((1.0-u*u)/x3);
    return;
  } 
  mu = (1.0-g2)/(1.0-g+2.0*g*random()/maxrand);
  mu = (1.0 + g2-mu*mu)/2.0/g; 
  mu2 = mu*mu;
  if ( fabs(w) < 0.9 ) {	/* avoid round-off errors in w */
    temp1=sqrt((1-mu2)/(1-w*w)/x3);
    t = mu * u + temp1 * (x1*u*w-x2*v);
    v = mu * v + temp1 * (x1*v*w+x2*u);
    w = mu * w - sqrt((1-mu2)*(1-w*w)/x3) * x1;
  } else {
    temp1=sqrt((1-mu2)/(1-v*v)/x3);
    t = mu * u + temp1 * (x1*u*v + x2*w);
    w = mu * w + temp1 * (x1*v*w - x2*u);
    v = mu * v - sqrt((1-mu2)*(1-v*v)/x3) * x1;
  }
  u = t;
  return;
}
 
int main (int argc, char *argv[] )
{
  int index;
  char c; 
  int  last_arg_processed=0;
  opterr = 0;

  if (argc==0) {
    print_help();
    return 1;
  }
  while ((c = getopt (argc, argv, "hg:p:z:l:m:o:")) != -1) {
    outfile="";
    last_arg_processed+=2;
    switch (c) { 
    case 'g':
      g = (double) atof(optarg);
      break;
    case 'p':
      photons = (long) atol(optarg);
      break;
    case 'l':
      t_layer = (double) atof(optarg);
      break;
    case 'm':
      N_max = (long) atol(optarg);
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'h':
      print_help();
      return 1;
    default:
      abort();
    }
  }
  if (outfile[0]=='\0') {
    ASCII=1;
  } else {
    ASCII=0;
    printf("outfile is %s\n", outfile);
    outPtr=fopen(outfile,"wb");
  }
  z_top[0]=0.0;
  z_bot[0]=-1.0*t_layer;
  z_top[1]=-1.0*t_layer;
  z_bot[1]=-1.0*5000.0/(1.0-g);  /*effectively infinite*/

  g2=g*g;
  srandom(time(0));
  maxrand=(exp(log(2.)*31.))-1;
  
  double g_0=g;
  double g2_0=g2*g2;

  for (i = 1; i <= photons; i++){
    launch ();
    while (count < N_max ) { 
      d = -log((random()+1.0)/(maxrand+1.0)); 
      while (d>0 && count < N_max) {
	move ();  /* move through as many layers as d allows */
      }
      count++;
      scatter ();
    }
  }	
  return 0;
}
