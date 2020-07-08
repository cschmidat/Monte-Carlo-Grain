#include "consts.h"

#define POSORNE(a,b) ((pts[a].nb[0]==b)?0:((pts[a].nb[1]==b)?1:2))
#define MEASP 1000
extern int n;
extern double ars[N];
extern char ck[N];
extern int *cs;



struct pt{
	int nb[3];
	unsigned int x;
	unsigned int y;
};

struct polydata{
	double ar[MEASP];
	int sides[MEASP];
};

//extern int write(char *suffix, struct pt *pts);
double dist(struct pt *,int,int);
static inline unsigned int mein_rand(void);
double polyar(int p1, char pos, struct pt *pts, int polyind);
int meas(struct pt *pts);
double meanar(struct pt *pts);
int intck(struct pt *pts,int p1,int p2,int q1,int q2);
/*
double polyar(int p1, char pos, struct pt *pts, int *cs){	//Each polygon is uniquely specified by directed edge (p1,p2) (not injectively)
	int cont1,cont2;
	unsigned piv[2];
	double ar=0;
	cont1=p1;
	cont2=pts[cont1].nb[pos];
	do{
		cs[3*cont1+pos]=1;
		ar+=((double)(int)(pts[cont1].y-pts[p1].y)*(int)(pts[cont2].x-pts[p1].x)-(double)(int)(pts[cont1].x-pts[p1].x)*(int)(pts[cont2].y-pts[p1].y))*0.5;
		pos=(pts[cont2].nb[0]==cont1)?1:((pts[cont2].nb[1]==cont1)?2:0);
		cont1=cont2;
		cont2=pts[cont2].nb[pos];
//		printf("%d %d\n",cont1,cont2);
	}
	while(cont1!=p1);
//	printf("%d %d %d\n%d %d %d\n%d %d %d\n%d %d %d\n%d %d %d\n%d %d %d\n",cs[3*1+0],cs[3*1+1],cs[3*1+2],cs[3*999215+0],cs[3*999215+1],cs[3*999215+2],cs[3*730+0],cs[3*730+1],cs[3*730+2],cs[3*1005+0],cs[3*1005+1],cs[3*1005+2],cs[3*1552+0],cs[3*1552+1],cs[3*1552+2],cs[3*1652+0],cs[3*1652+1],cs[3*1652+2]);
	return ar;
}
*/


int meas(struct pt *pts){
	int i,j,arind;
	int vis[N]={0};
	for(i=0;i<3*n;i++){
		j=(i-i%3)/3;
		arind=cs[i];
		if(vis[arind]==0){
			polyar(j,i%3,pts,arind);
			vis[arind]=1;
		}
		if(cs[i]!=arind){
			printf("Error\n"); exit(1);}
		}
		for(i=0;i<N;i++){
			if(vis[i]==0){	ck[i]=0; ars[i]=0;}
		}
	return 1;
}

double meanar(struct pt *pts){//algo copied from init(), writes to ck[] and ars[]
	double allar;
	int i,polcount=0;
	meas(pts);
	for(i=0;i<N;i++){
		if(ck[i]!=0){	allar+=ars[i]; polcount++;}
	}
	return(allar/polcount);
}


double polyar(int p1, char pos, struct pt *pts, int polyind){
	int cont1,cont2;
	char sides=0;
	unsigned piv[2];
	double ar=0;
	cont1=p1;
	cont2=pts[cont1].nb[pos];
//	if(p1==42765){printf("p1: %d\n",p1);}
	do{
		if(sides>50){
			printf("Too many sides starting at p1: %d pos: %d n: %d ind: %d\n",p1,pos,n,polyind);
			write("toomuch",pts);
			exit(1);
		}
		cs[3*cont1+pos]=polyind;
		ar+=((double)(int)(pts[cont1].y-pts[p1].y)*(int)(pts[cont2].x-pts[p1].x)-(double)(int)(pts[cont1].x-pts[p1].x)*(int)(pts[cont2].y-pts[p1].y))*0.5;
	pos=(pts[cont2].nb[0]==cont1)?1:((pts[cont2].nb[1]==cont1)?2:0);
		cont1=cont2;
		cont2=pts[cont2].nb[pos];
//		printf("%d %d\n",cont1,cont2);

		sides++;
	}
	while(cont1!=p1);
	ck[polyind]=sides;
	ars[polyind]=ar;
	return ar;
}

double meansides(struct pt *pts){
	int i, sidesnum=0,polanz=0;

	for(i=0;i<n;i++){
		sidesnum+=ck[i];
		if(ck[i]!=0) polanz++;
	}

	return((double)sidesnum/polanz);
}
int gethist(struct pt *pts,int *histp, double *meansizep){ //30 points built in!
	int i,j;
	double meansum[30]={0};
	for(i=0;i<30;i++){
		for(j=0;j<N;j++){
			if(ck[j]==i){histp[i]++;
			meansum[i]+=ars[j];}
		}
		if(meansum[i]!=0) meansizep[i]=meansum[i]/histp[i];
	}
}
double meandist(struct pt *pts){
	int i,j, distnum=0;
	double distmean,distsum=0;
	for(i=0;i<n;i++){
		for(j=0;j<3;j++){
			distsum+=dist(pts,i,pts[i].nb[j]);
			distnum++;
		}
	}
	distmean=distsum/distnum;
//	printf("Distnum: %d Distsum: %lf Distmean: %lf\n",distnum,distsum,distmean);
	return(distmean);
}
/*
	Reads polygon areas, cons and sides, returns #(polygons)
*/
double init(struct pt *pts){
	int i,polcount=1;
	double ar=0;
	double sl=0;
	for(i=0;i<3*n;i++){
		if(cs[i]==0){
			ar+=polyar((i-i%3)/3,i%3,pts,polcount);
			polcount++;
//			printf("i: %d\n",i);
		}
	}
	for(i=1;i<polcount;i++){
		sl+=ck[i];
	}
	ck[polcount]=-1;
	printf("%lf\n",sl/(polcount-1));
	return(polcount-1);
		
}

/*
Resolves intersections
*/

int cross(struct pt *pts){
	int i,j,k=0;
	printf("Starting cross\n");
	for(i=0;i<n;i++){
		for(j=0;j<3;j++){
//			printf("i: %d j: %d\n",i,j);
			if(dist(pts,i,pts[i].nb[j])<200000&&intersects(pts,i,pts[i].nb[j])){
				printf("Dist: %lf ",dist(pts,i,pts[i].nb[j]));
				pts[i].x+=2*(pts[pts[i].nb[j]].x-pts[i].x)+(int)(mein_rand())%5;
				pts[i].y+=2*(pts[pts[i].nb[j]].y-pts[i].y)+(int)(mein_rand())%5;
				printf("Resolved %d with %d\n",i,pts[i].nb[j]);
				printf("Koord nachher: %u %u\n",pts[i].x,pts[i].y);
				j=0;
				k++;
//				printf("Anzahl: %d\n",k);
			}
		}
	}
}

int intersects(struct pt *pts,int i,int j){
	char posi,posj;
	posi=POSORNE(i,j);
	posj=POSORNE(j,i);
	return(intck(pts,i,pts[i].nb[(posi+1)%3],j,pts[j].nb[(posj+1)%3])||intck(pts,i,pts[i].nb[(posi+1)%3],j,pts[j].nb[(posj+2)%3])||intck(pts,i,pts[i].nb[(posi+2)%3],j,pts[j].nb[(posj+1)%3])||intck(pts,i,pts[i].nb[(posi+2)%3],j,pts[j].nb[(posj+2)%3]));
}

int intck(struct pt *pts,int p1,int p2,int q1,int q2){
	double a1,a2,b1,b2;
	double xi,yi;
	b1=(double)(int)(pts[p2].y-pts[p1].y)/((int)(pts[p2].x-pts[p1].x));
	b2=(double)(int)(pts[q2].y-pts[q1].y)/((int)(pts[q2].x-pts[q1].x));
	a1=pts[p1].y-b1*pts[p1].x;
	a2=pts[q1].y-b2*pts[q1].x;
	xi=-(a1-a2)/(b1-b2);
	yi=a1+b1*xi;
//	if(p1==518606||q1==518606)
//		printf("\n\np1: %d q1: %d\np1: %u %u\np2: %u %u\n q1: %u %u\n q2: %u %u\nb1: %lf\nb2: %lf\na1: %lf\na2: %lf\nxi: %lf\nyi:%lf\n",p1,q1,pts[p1].x,pts[p1].y,pts[p2].x,pts[p2].y,pts[q1].x,pts[q1].y,pts[q2].x,pts[q2].y,b1,b2,a1,a2,xi,yi);
	return((pts[p1].x-xi)*(xi-pts[p2].x)>=0&&(pts[q1].x-xi)*(xi-pts[q2].x)>=0&&(pts[p1].y-yi)*(yi-pts[p2].y)>=0&&(pts[q1].y-yi)*(yi-pts[q2].y)>=0);
}
/*int meas(struct pt *pts, int *allcs){
	int i;
	double arges;
	for(i=0;i<3*n;i++){
	if(allcs[i]==0)
		arges+=polyar((i-i%3)/3,i%3,pts,allcs);
//		printf("i: %d\n",i);
	}
	printf("Areatest: %lf\n",arges);
}*/

