#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#define STEPS 20000000000
#define JUMPAR 1.8e11//initially 1.8e11
#define KTEMPAN 10e3 //initially 7e3
#define MINAR 5e11 //initially 5e12
#define T1TRIES 1
#define MEASP 1000 //measurement points
#define MEASF 1000
#define RATEP 10000 //#(polygons whose growth rate we'll calculate)
#define ABST 10000
#define BASEDIR "data/"
#include "meas.c"

#define unlikely(expr) __builtin_expect(!!(expr),0)
#define likely(expr) __builtin_expect(!!(expr),1)
#define POSORNE(a,b) ((pts[a].nb[0]==b)?0:((pts[a].nb[1]==b)?1:2))
#define DOUBCON(a) (pts[a].nb[0]==pts[a].nb[1]||pts[a].nb[0]==pts[a].nb[2]||pts[a].nb[1]==pts[a].nb[2])

/*
Simulation: 
Events: T1...Swap
T2...Disappearance of small grain
Jump length is chosen s.t. 30% of jumps pass without considering Boltzmann factor
Temperatur is chosen s.t. points fluctuate an average of 1% of the typical initial length.
*/

int n=N;
double jumpar=JUMPAR,minar=MINAR;//jumpar: maximal allowed "jump area". minar: minimal area for triangles; sqrt(minar): minimal t1 length
double ktemp=KTEMPAN;
long long trisucc=0, succ=0, expsucc=0;
double distall=0;
unsigned int MEIN_RAND_STATE=0;
double ars[N];	//areas of polygons
char ck[N];	//numbers of sides
int *cs;	//connection
double mdists[MEASP],msides[MEASP],ratedata[RATEP][4],timemeas[MEASP],rate500000[N/2][14];
int hist[MEASP][30];
double meansize[MEASP][30];
int spolylist[]={10000,4000,1000,2001,3002,4003,5004,6005,1762,2434,5339,70927,1123,9003,8327,90473,308328};
struct polydata *spolydata;
double t=0;


static inline unsigned int mein_rand(void){
    MEIN_RAND_STATE=1664525L*MEIN_RAND_STATE+1013904223L;
    return MEIN_RAND_STATE;
}
inline int setvars(void){
	minar=MINAR*N/n;
	jumpar=JUMPAR*N/n;
	ktemp=KTEMPAN*sqrt((double)n/N);
}

inline char jumptest(double dx, double dy, int dist1x, int dist1y, int dist2x, int dist2y, int dist3x, int dist3y){//test whether jump shall be made
	double x=(fabs(dy*dist1x-dx*dist1y)+fabs(dy*dist2x-dx*dist2y)+fabs(dy*dist3x-dx*dist3y))*0.5; //Area (3 triangles)
	double distbef,distaf;
//	printf("%lf\n",x);
	if(x<=jumpar){
	trisucc++;
//	distall+=sqrt(dx*dx+dy*dy);
	distbef=sqrt((double)dist1x*dist1x+(double)dist1y*dist1y)+sqrt((double)dist2x*dist2x+(double)dist2y*dist2y)+sqrt((double)dist3x*dist3x+(double)dist3y*dist3y);
	distaf=sqrt((dist1x+dx)*(dist1x+dx)+(dist1y+dy)*(dist1y+dy))+sqrt((dist2x+dx)*(dist2x+dx)+(dist2y+dy)*(dist2y+dy))+sqrt((dist3x+dx)*(dist3x+dx)+(dist3y+dy)*(dist3y+dy));
//	printf("distbef: %lf, distaf: %lf\n",distbef,distaf);
	if(exp((distbef-distaf)/ktemp)>=((double) mein_rand())/4294967295){
//	if(distbef>distaf){	//T=0
		expsucc++;
		if(distbef>distaf) succ++;
		return 1;
//		printf("dist: %lf\n",sqrt(dx*dx+dy*dy));
	}
	else return 0;
	}
	else return 0;
}

int conswap(int pfix,int por,int pnew,struct pt *pts){ //pfix's connection to por->connection to pnew, can be extended so pnew->pfix, but instead of what?
	int buf,pos;
	pos=(por==pts[pfix].nb[0])?0:((por==pts[pfix].nb[1])?1:2);
	pts[pfix].nb[pos]=pnew;
}
struct pt * rm(int i, struct pt *pts){
	if(n%100000==0)
	pts=(struct pt *) realloc(pts, (n-1)*sizeof(struct pt));
	return pts;
}

//double polyar(int p1, char pos, struct pt *pts, int *cs);

int meas(struct pt *pts);

int tabrep(int inew, int iold, struct pt *pts){ //inew->iold, iold's data is overwritten!
	int j;
	for(j=0;j<3;j++){
		pts[iold].nb[j]=pts[inew].nb[j];
		conswap(pts[iold].nb[j],inew,iold,pts);
		cs[3*iold+j]=cs[3*inew+j];
	}
	pts[iold].x=pts[inew].x;
	pts[iold].y=pts[inew].y;
}

double dist(struct pt *pts,int p1,int p2){
	return sqrt(pow((int)(pts[p1].x-pts[p2].x),2)+pow((int)(pts[p1].y-pts[p2].y),2));
}

int deltri(int i, int pt1, int pt2, int noti, int notpt1, int notpt2, struct pt *pts, int leavecs){
	int ni,np1,np2; //real indices of points above
	ni=pts[i].nb[noti];
	np1=pts[pt1].nb[notpt1];
	np2=pts[pt2].nb[notpt2];

	//Zuerst übernimmt i neue Nachbarn, diese übernehmen i:
	pts[i].nb[(noti+1)%3]=pts[pt1].nb[notpt1];
	pts[i].nb[(noti+2)%3]=pts[pt2].nb[notpt2];
	conswap(pts[pt1].nb[notpt1],pt1,i,pts);
	conswap(pts[pt2].nb[notpt2],pt2,i,pts);
	
	if(pts[ni].nb[POSORNE(ni,i)]!=i||pts[np1].nb[POSORNE(np1,i)]!=i||pts[np2].nb[POSORNE(np2,i)]!=i){
		printf("Achtung, Posorne passt ned! i: %d p1: %d p2: %d ni: %d n1: %d n2: %d\n",i,pt1,pt2,ni,np1,np2);
		exit(1);
	}
	if(leavecs%2==0) polyar(ni,POSORNE(ni,i),pts,cs[3*ni+POSORNE(ni,i)]);//Achtung: Was passiert beim Vertauschen?
	if((leavecs-leavecs%2)%4!=2) polyar(np1,POSORNE(np1,i),pts,cs[3*np1+POSORNE(np1,i)]);
	if((leavecs-(leavecs-leavecs%2)%4)%8!=4) polyar(np2,POSORNE(np2,i),pts,cs[3*np2+POSORNE(np2,i)]);
	
	if((pt1==n-1&&pt2==n-2)||(pt1==n-2&&pt2==n-1))
		;
	else if(pt1==n-1)
		tabrep(n-2,pt2,pts);
	else if(pt2==n-2)
		tabrep(n-1,pt1,pts);
	else if(pt1==n-2)
		tabrep(n-1,pt2,pts);
	else if(pt2==n-1)
		tabrep(n-2,pt1,pts);
	else{
		tabrep(n-1,pt1,pts);
		tabrep(n-2,pt2,pts);
	}
	pts=rm(pt1,pts);
	pts=rm(pt2,pts);
	n-=2;
	setvars();
	
//	printf("%lf %lf %lf\n",ktemp, minar, jumpar);
}
int sort(int i, struct pt *pts){
	int dist1x, dist1y, dist2x, dist2y,dist3x,dist3y;
	double spo1, spo2, sp1, sp2;	//spoi: angle bt. d1^t and d(i+1), >0=>first side
	int buf;
	
	dist1x=pts[i].x-pts[pts[i].nb[0]].x;
	dist1y=pts[i].y-pts[pts[i].nb[0]].y;
	dist2x=pts[i].x-pts[pts[i].nb[1]].x;
	dist2y=pts[i].y-pts[pts[i].nb[1]].y;
	dist3x=pts[i].x-pts[pts[i].nb[2]].x;
	dist3y=pts[i].y-pts[pts[i].nb[2]].y;
	
	spo1=(-(double)dist2x*dist1y+(double)dist1x*dist2y)/(sqrt(pow(dist1x,2)+pow(dist1y,2))*sqrt(pow(dist2x,2)+pow(dist2y,2)));
	spo2=(-(double)dist3x*dist1y+(double)dist1x*dist3y)/(sqrt(pow(dist1x,2)+pow(dist1y,2))*sqrt(pow(dist3x,2)+pow(dist3y,2)));

	if(spo1>=0&&spo2<=0)	return 0;
	else if(spo1<0&&spo2>0){
		buf=pts[i].nb[1];
		pts[i].nb[1]=pts[i].nb[2];
		pts[i].nb[2]=buf;
		return 1;
	}
	else if(spo1>=0&&spo2>=0){
		sp1=((double)dist2x*dist1x+(double)dist2y*dist1y)/(sqrt(pow(dist1x,2)+pow(dist1y,2))*sqrt(pow(dist2x,2)+pow(dist2y,2)));
		sp2=((double)dist3x*dist1x+(double)dist3y*dist1y)/(sqrt(pow(dist1x,2)+pow(dist1y,2))*sqrt(pow(dist3x,2)+pow(dist3y,2)));
		if(sp1<sp2){
			buf=pts[i].nb[1];
			pts[i].nb[1]=pts[i].nb[2];
			pts[i].nb[2]=buf;
			return 1;
		}
		else return 0;
	}
	else if(spo1<=0&&spo2<=0){
		sp1=((double)dist2x*dist1x+(double)dist2y*dist1y)/(sqrt(pow(dist1x,2)+pow(dist1y,2))*sqrt(pow(dist2x,2)+pow(dist2y,2)));
		sp2=((double)dist3x*dist1x+(double)dist3y*dist1y)/(sqrt(pow(dist1x,2)+pow(dist1y,2))*sqrt(pow(dist3x,2)+pow(dist3y,2)));
		if(sp1>sp2){
			buf=pts[i].nb[1];
			pts[i].nb[1]=pts[i].nb[2];
			pts[i].nb[2]=buf;
			return 1;
		}
		else return 0;
	}
	
}
int fixdoubconn(int left1, int left2, struct pt *pts){
	int nleft1, nleft2,dleft1,dleft2;
	nleft1=(pts[left1].nb[0]==pts[left1].nb[1])?pts[left1].nb[2]:((pts[left1].nb[0]==pts[left1].nb[2])?pts[left1].nb[1]:pts[left1].nb[0]);
	nleft2=(pts[left2].nb[0]==pts[left2].nb[1])?pts[left2].nb[2]:((pts[left2].nb[0]==pts[left2].nb[2])?pts[left2].nb[1]:pts[left2].nb[0]);
	if(left1==left2){ //nleft1==nleft2: Single one left, needs to be deleted
		dleft1=(pts[nleft1].nb[0]==left1)?1:((pts[nleft1].nb[1]==left1)?2:0);
		dleft2=pts[nleft1].nb[(dleft1+1)%3];
		dleft1=pts[nleft1].nb[dleft1];
		printf("\n\nACHTUNG: left1: %d left2: %d nleft1: %d nleft2: %d dleft1: %d dleft2: %d\n\n",left1,left2,nleft1,nleft2,dleft1,dleft2);
		conswap(dleft1,nleft1,dleft2,pts);
		conswap(dleft2,nleft1,dleft1,pts);
		tabrep(n-1,nleft1,pts);
		tabrep(n-2,left1,pts);
		n-=2;
		setvars();
		if(DOUBCON(dleft1)||DOUBCON(dleft2)){	//just checked for "end-case" 2
			printf("Error at three triangs: left1: %d left2: %d\n",left1,left2);
			fixdoubconn(dleft1,dleft2,pts);
		}
	}
	else{
		printf("\n\nACHTUNG2: left1: %d left2: %d nleft1: %d nleft2: %d\n\n",left1,left2,nleft1,nleft2);
		printf("Polygon: %d %d\n", cs[3*left1+(POSORNE(left1,nleft1)+2)%3],cs[3*left2+(POSORNE(left2,nleft2)+2)%3]);
		ck[cs[3*left1+(POSORNE(left1,nleft1)+2)%3]]=0;
		pts=rm(left1,pts);
		pts=rm(left2,pts);
		conswap(nleft1,left1,nleft2,pts);
		conswap(nleft2,left2,nleft1,pts);
		tabrep(n-1,left1,pts);
		tabrep(n-2,left2,pts);
		n-=2;
		setvars();
		if(DOUBCON(nleft1)||DOUBCON(nleft2)){
			printf("Error at double triangs's else: left1: %d left2: %d\n",left1,left2);
			fixdoubconn(nleft1,nleft2,pts);
		}
		else{
			polyar(nleft1,POSORNE(nleft1,nleft2),pts,cs[3*nleft1+POSORNE(nleft1,nleft2)]);
			polyar(nleft2,POSORNE(nleft2,nleft1),pts,cs[3*nleft2+POSORNE(nleft2,nleft1)]);
		}
	}
	return 0;
}


/* 
	Implemented: "Normal" double triangles, Three triangles, double triangle & quadrilateral
	Points in triangles: i,pt1,pt2,pt3
	Points outside: left1, left2
*/
int doubtri(struct pt *pts,int i,int pt1,int pt2,int noti,int notpt1,int notpt2){
	char check;
	double tri2ar;
	int left1,left2,pt3;
	int nleft1, nleft2,buf2,dleft1,dleft2;
//	printf("Error: Doppeldreieck!\nKoordinaten: %d %d %d\n",i,pt1,pt2);
//	write("midsneu","nachneu",points,nach);
//	exit(1);
	check=(pts[pt1].nb[notpt1]==pts[pt2].nb[notpt2])?0:((pts[i].nb[noti]==pts[pt1].nb[notpt1])?2:1);//which point is not in the other triangle?
	switch (check){
		case 0://i is not in other triangle
			{
			tri2ar=fabs((double)(int)(pts[pt2].x-pts[pts[pt2].nb[notpt2]].x)*(int)(pts[pt1].y-pts[pts[pt2].nb[notpt2]].y)-(double)(int)(pts[pt1].x-pts[pts[pt2].nb[notpt2]].x)*(int)(pts[pt2].y-pts[pts[pt2].nb[notpt2]].y))*0.5;
			left1=pts[i].nb[noti];
			left2=!(pts[pts[pt2].nb[notpt2]].nb[0]==pt1||pts[pts[pt2].nb[notpt2]].nb[0]==pt2)?pts[pts[pt2].nb[notpt2]].nb[0]:(!(pts[pts[pt2].nb[notpt2]].nb[1]==pt1||pts[pts[pt2].nb[notpt2]].nb[1]==pt2)?pts[pts[pt2].nb[notpt2]].nb[1]:pts[pts[pt2].nb[notpt2]].nb[2]);
			break;
			}
		case 1://pt1 is not in triangle
			{
			tri2ar=fabs((double)(int)(pts[pt2].x-pts[i].x)*(int)(pts[pts[pt2].nb[notpt2]].y-pts[i].y)-(double)(int)(pts[pts[pt2].nb[notpt2]].x-pts[i].x)*(int)(pts[pt2].y-pts[i].y))*0.5;
			left1=pts[pt1].nb[notpt1];
			left2=!(pts[pts[pt2].nb[notpt2]].nb[0]==i||pts[pts[pt2].nb[notpt2]].nb[0]==pt2)?pts[pts[pt2].nb[notpt2]].nb[0]:(!(pts[pts[pt2].nb[notpt2]].nb[1]==i||pts[pts[pt2].nb[notpt2]].nb[1]==pt2)?pts[pts[pt2].nb[notpt2]].nb[1]:pts[pts[pt2].nb[notpt2]].nb[2]);
			break;
			}
		case 2://pt2 is not in triangle
			{
			tri2ar=fabs((double)(int)(pts[pts[pt1].nb[notpt1]].x-pts[i].x)*(int)(pts[pt1].y-pts[i].y)-(double)(int)(pts[pt1].x-pts[i].x)*(int)(pts[pts[pt1].nb[notpt1]].y-pts[i].y))*0.5;
			left1=pts[pt2].nb[notpt2];
			left2=!(pts[pts[pt1].nb[notpt1]].nb[0]==pt1||pts[pts[pt1].nb[notpt1]].nb[0]==i)?pts[pts[pt1].nb[notpt1]].nb[0]:(!(pts[pts[pt1].nb[notpt1]].nb[1]==pt1||pts[pts[pt1].nb[notpt1]].nb[1]==i)?pts[pts[pt1].nb[notpt1]].nb[1]:pts[pts[pt1].nb[notpt1]].nb[2]);
			break;
			}
		}
			
//		printf("Area of second triangle: %lf\n",tri2ar);
		if(tri2ar<minar){
			switch (check){
				case 0://i is not in other triangle
					{
					pt3=pts[pt2].nb[notpt2];
					ck[cs[3*i+(noti+2)%3]]=0;
					ck[cs[3*pt3+(POSORNE(pt3,left2)+2)%3]]=0;
					conswap(left1,i,left2,pts);
					conswap(left2,pts[pt2].nb[notpt2],left1,pts);
					break;
					}
				case 1://pt1 is not in triangle
					{
					pt3=pts[pt2].nb[notpt2];
					ck[cs[3*pt1+(notpt1+2)%3]]=0;
					ck[cs[3*pt3+(POSORNE(pt3,left2)+2)%3]]=0;
					conswap(left1,pt1,left2,pts);
					conswap(left2,pts[pt2].nb[notpt2],left1,pts);
					break;
					}
				case 2://pt2 is not in trianglew
					{
					pt3=pts[pt1].nb[notpt1];
					ck[cs[3*pt2+(notpt2+2)%3]]=0;
					ck[cs[3*pt3+(POSORNE(pt3,left2)+2)%3]]=0;
					conswap(left1,pt2,left2,pts);
					conswap(left2,pts[pt1].nb[notpt1],left1,pts);
					break;
					}
			}
			if(i>n-5||pt1>n-5||pt2>n-5||pt3>n-5){
			printf("Oh oh (not that bad ;) ) i: %d pt1: %d pt2: %d pt3: %d",i,pt1,pt2,pt3);
			write("ohoh",pts);
			exit(1);
			}

			else{
				printf("left1: %d left2: %d pt3: %d check: %d noti: %d notpt1: %d notpt2: %d\n",left1,left2,pt3,check,noti,notpt1,notpt2);
				tabrep(n-1,pt1,pts);
				tabrep(n-2,pt2,pts);
				tabrep(n-3,pt3,pts);
				tabrep(n-4,i,pts);
				pts=rm(pt1,pts);
				pts=rm(pt2,pts);
				pts=rm(pt3,pts);
				pts=rm(i,pts);
				n-=4;
				setvars();
				
				if(pts[left1].nb[0]==pts[left1].nb[1]||pts[left1].nb[0]==pts[left1].nb[2]||pts[left1].nb[1]==pts[left1].nb[2]||pts[left2].nb[0]==pts[left2].nb[1]||pts[left2].nb[0]==pts[left2].nb[2]||pts[left2].nb[1]==pts[left2].nb[2]){
					fixdoubconn(left1,left2,pts);			
				}
				
				else{
					polyar(left1,POSORNE(left1,left2),pts,cs[3*left1+POSORNE(left1,left2)]);
					polyar(left2,POSORNE(left2,left1),pts,cs[3*left2+POSORNE(left2,left1)]);
					return 2;
				}
			}
		}
}
		
					
			
//	if(
//		tri2ar=fabs((double)(int)(points[pt2][0]-points[i][0])*(int)(points[pt1][1]-points[i][1])-(double)(int)(points[pt1][0]-points[i][0])*(int)(points[pt2][1]-points[i][1]))*0.5;

int t1(int i, double dist1, double dist2, double dist3, struct pt *pts){ //Only check in special directions (dependent on dx,dy)?
//	double minlen=sqrt(minar)/1.52;
//T1-Nachbar zufällig auswählen
	double minlen=minar/2.31;//*2.31
//	if(i==1) printf("Minlen: %lf\n",sqrt(minlen));
	int ran,j,k,t1not[2];
	double dst[3];
	int inpos,buf,inapos,incpos,pt2,badpos,sneigh;
	if(i==57901&&n==500002)	printf("n1: %d dist1: %lf n2: %d dist2: %lf n3: %d dist3: %lf minlen: %lf\n",pts[i].nb[0],dist1,pts[i].nb[1],dist2,pts[i].nb[2],dist3,minlen);
	if(dist1<minlen||dist2<minlen||dist3<minlen){
		dst[0]=dist1;dst[1]=dist2;dst[2]=dist3;
		while(dst[sneigh=(mein_rand()%3)]>minlen)
		;
//	int sneigh=(dist1<minlen) ? 0 : ((dist2<minlen)? 1 : ((dist3<minlen) ? 2:-1));	//Sry about that^^
	
//	if(sneigh!=-1){ //For now, I won't care about the directions of the swap
//		printf("%d, i:%d, Nachbarn: %d %d %d\n", sneigh,i,nach[i][0], nach[i][1], nach[i][2]);
		pt2=pts[i].nb[sneigh];
		inpos=(i==pts[pt2].nb[0])? 0:((i==pts[pt2].nb[1])?1:2); //nach[i][sneigh] is point2, i is point1, nach[point2][inpos] is point1
		if(pts[i].nb[(sneigh+1)%3]==pts[pt2].nb[(inpos+1)%3]||pts[i].nb[(sneigh+2)%3]==pts[pt2].nb[(inpos+1)%3]||pts[i].nb[(sneigh+1)%3]==pts[pt2].nb[(inpos+2)%3]||pts[i].nb[(sneigh+2)%3]==pts[pt2].nb[(inpos+2)%3])
			return -1; //Bei Dreieck nix machen!
			
			//prevent double triangles
			//check inpos+1 and sneigh+1...erledigt
		if(pts[pts[pt2].nb[(inpos+1)%3]].nb[0]==pts[i].nb[(sneigh+1)%3]||pts[pts[pt2].nb[(inpos+1)%3]].nb[1]==pts[i].nb[(sneigh+1)%3]||pts[pts[pt2].nb[(inpos+1)%3]].nb[2]==pts[i].nb[(sneigh+1)%3]){//Is da Viereck
			badpos=(pts[pts[pt2].nb[(inpos+1)%3]].nb[0]==pts[i].nb[(sneigh+1)%3])?0:((pts[pts[pt2].nb[(inpos+1)%3]].nb[1]==pts[i].nb[(sneigh+1)%3])?1:2);
			if(pts[pts[pt2].nb[(inpos+1)%3]].nb[(badpos+1)%3]==pts[pts[i].nb[(sneigh+1)%3]].nb[0]||pts[pts[pt2].nb[(inpos+1)%3]].nb[(badpos+1)%3]==pts[pts[i].nb[(sneigh+1)%3]].nb[1]||pts[pts[pt2].nb[(inpos+1)%3]].nb[(badpos+1)%3]==pts[pts[i].nb[(sneigh+1)%3]].nb[2]||pts[pts[pt2].nb[(inpos+1)%3]].nb[(badpos+2)%3]==pts[pts[i].nb[(sneigh+1)%3]].nb[0]||pts[pts[pt2].nb[(inpos+1)%3]].nb[(badpos+2)%3]==pts[pts[i].nb[(sneigh+1)%3]].nb[1]||pts[pts[pt2].nb[(inpos+1)%3]].nb[(badpos+2)%3]==pts[pts[i].nb[(sneigh+1)%3]].nb[2])
				return -1;
		}
		
		//check inpos+1 and sneigh+2
		if(pts[pts[pt2].nb[(inpos+1)%3]].nb[0]==pts[i].nb[(sneigh+2)%3]||pts[pts[pt2].nb[(inpos+1)%3]].nb[1]==pts[i].nb[(sneigh+2)%3]||pts[pts[pt2].nb[(inpos+1)%3]].nb[2]==pts[i].nb[(sneigh+2)%3]){
			badpos=(pts[pts[pt2].nb[(inpos+1)%3]].nb[0]==pts[i].nb[(sneigh+2)%3])?0:((pts[pts[pt2].nb[(inpos+1)%3]].nb[1]==pts[i].nb[(sneigh+2)%3])?1:2);
			if(pts[pts[pt2].nb[(inpos+1)%3]].nb[(badpos+1)%3]==pts[pts[i].nb[(sneigh+2)%3]].nb[0]||pts[pts[pt2].nb[(inpos+1)%3]].nb[(badpos+1)%3]==pts[pts[i].nb[(sneigh+2)%3]].nb[1]||pts[pts[pt2].nb[(inpos+1)%3]].nb[(badpos+1)%3]==pts[pts[i].nb[(sneigh+2)%3]].nb[2]||pts[pts[pt2].nb[(inpos+1)%3]].nb[(badpos+2)%3]==pts[pts[i].nb[(sneigh+2)%3]].nb[0]||pts[pts[pt2].nb[(inpos+1)%3]].nb[(badpos+2)%3]==pts[pts[i].nb[(sneigh+2)%3]].nb[1]||pts[pts[pt2].nb[(inpos+1)%3]].nb[(badpos+2)%3]==pts[pts[i].nb[(sneigh+2)%3]].nb[2])
				return -1;
		
		}
					//check inpos+2 and sneigh+1
		if(pts[pts[pt2].nb[(inpos+2)%3]].nb[0]==pts[i].nb[(sneigh+1)%3]||pts[pts[pt2].nb[(inpos+2)%3]].nb[1]==pts[i].nb[(sneigh+1)%3]||pts[pts[pt2].nb[(inpos+2)%3]].nb[2]==pts[i].nb[(sneigh+1)%3]){
			badpos=(pts[pts[pt2].nb[(inpos+2)%3]].nb[0]==pts[i].nb[(sneigh+1)%3])?0:((pts[pts[pt2].nb[(inpos+2)%3]].nb[1]==pts[i].nb[(sneigh+1)%3])?1:2);
			if(pts[pts[pt2].nb[(inpos+2)%3]].nb[(badpos+1)%3]==pts[pts[i].nb[(sneigh+1)%3]].nb[0]||pts[pts[pt2].nb[(inpos+2)%3]].nb[(badpos+1)%3]==pts[pts[i].nb[(sneigh+1)%3]].nb[1]||pts[pts[pt2].nb[(inpos+2)%3]].nb[(badpos+1)%3]==pts[pts[i].nb[(sneigh+1)%3]].nb[2]||pts[pts[pt2].nb[(inpos+2)%3]].nb[(badpos+2)%3]==pts[pts[i].nb[(sneigh+1)%3]].nb[0]||pts[pts[pt2].nb[(inpos+2)%3]].nb[(badpos+2)%3]==pts[pts[i].nb[(sneigh+1)%3]].nb[1]||pts[pts[pt2].nb[(inpos+2)%3]].nb[(badpos+2)%3]==pts[pts[i].nb[(sneigh+1)%3]].nb[2])
				return -1;
		}
		
		//check inpos+2 and sneigh+2
		if(pts[pts[pt2].nb[(inpos+2)%3]].nb[0]==pts[i].nb[(sneigh+2)%3]||pts[pts[pt2].nb[(inpos+2)%3]].nb[1]==pts[i].nb[(sneigh+2)%3]||pts[pts[pt2].nb[(inpos+2)%3]].nb[2]==pts[i].nb[(sneigh+2)%3]){
			badpos=(pts[pts[pt2].nb[(inpos+2)%3]].nb[0]==pts[i].nb[(sneigh+2)%3])?0:((pts[pts[pt2].nb[(inpos+2)%3]].nb[1]==pts[i].nb[(sneigh+2)%3])?1:2);
			if(pts[pts[pt2].nb[(inpos+2)%3]].nb[(badpos+1)%3]==pts[pts[i].nb[(sneigh+2)%3]].nb[0]||pts[pts[pt2].nb[(inpos+2)%3]].nb[(badpos+1)%3]==pts[pts[i].nb[(sneigh+2)%3]].nb[1]||pts[pts[pt2].nb[(inpos+2)%3]].nb[(badpos+1)%3]==pts[pts[i].nb[(sneigh+2)%3]].nb[2]||pts[pts[pt2].nb[(inpos+2)%3]].nb[(badpos+2)%3]==pts[pts[i].nb[(sneigh+2)%3]].nb[0]||pts[pts[pt2].nb[(inpos+2)%3]].nb[(badpos+2)%3]==pts[pts[i].nb[(sneigh+2)%3]].nb[1]||pts[pts[pt2].nb[(inpos+2)%3]].nb[(badpos+2)%3]==pts[pts[i].nb[(sneigh+2)%3]].nb[2])
				return -1;
		
		}

		ran=1+(mein_rand())%2;
		
		if(exp((dist(pts,i,pts[i].nb[(ran+sneigh)%3])+dist(pts,pt2,pts[pt2].nb[(ran+inpos)%3])-dist(pts,i,pts[pt2].nb[(ran+inpos)%3])-dist(pts,pt2,pts[i].nb[(ran+sneigh)%3]))/(ktemp))>=((double) mein_rand())/4294967295){
//			printf("dist1: %lf dist2: %lf dist3: %lf dist4: %lf\n i: %d pt2: %d A: %d C: %d\n",dist(points,i,nach[i][(ran1+sneigh)%3]),dist(points,pt2,nach[pt2][(ran2+inpos)%3]),dist(points,i,nach[pt2][(ran2+inpos)%3]),dist(points,pt2,nach[i][(ran1+sneigh)%3]),i,pt2,nach[i][(ran1+sneigh)%3],nach[pt2][(ran2+inpos)%3]);
			buf=pts[i].nb[(sneigh+ran)%3];
			pts[i].nb[(sneigh+ran)%3]=pts[pt2].nb[(inpos+ran)%3];
//			if(n<305000&&pts[i].x<0.21e9&&pts[i].y>3.3e9&&pts[i].x>2e8)
//				printf("pt2(inpos+ran)mod3: %d\n",pts[pt2].nb[(inpos+ran)%3]);
			pts[pt2].nb[(inpos+ran)%3]=buf;	//First neighbour changed for pts 1 and 2, we still have to change their neighbours' neighbours. C is now nach[point1][sneigh+ran], A is now nach[point2][inpos+ran] and buf.
			inapos=(i==pts[buf].nb[0])? 0:((i==pts[buf].nb[1])?1:2); //Search point1 in A
			pts[buf].nb[inapos]=pt2; //A now points to point2
			incpos=(pt2==pts[pts[i].nb[(sneigh+ran)%3]].nb[0])? 0:((pt2==pts[pts[i].nb[(sneigh+ran)%3]].nb[1])?1:2); //Search point2 in C
			pts[pts[i].nb[(sneigh+ran)%3]].nb[incpos]=i; //C now points to point1=i

//			if(n>280000&&n<305000&&pts[i].x<0.21e9&&pts[i].y>3.3e9&&pts[i].x>2e8)
//				write("badbfsrt",pts);

			//Sort:
			buf=pts[i].nb[(sneigh+1)%3];
			pts[i].nb[(sneigh+1)%3]=pts[i].nb[(sneigh+2)%3];
			pts[i].nb[(sneigh+2)%3]=buf;
			buf=pts[pt2].nb[(inpos+1)%3];
			pts[pt2].nb[(inpos+1)%3]=pts[pt2].nb[(inpos+2)%3];
			pts[pt2].nb[(inpos+2)%3]=buf;

			for(j=0;j<3;j++)
				polyar(pts[pt2].nb[(inpos+1)%3],j,pts,cs[3*pts[pt2].nb[(inpos+1)%3]+j]);
			for(j=0;j<3;j++)
				polyar(pts[pt2].nb[(inpos+2)%3],j,pts,cs[3*pts[pt2].nb[(inpos+2)%3]+j]);
			for(j=0;j<3;j++)
				polyar(pts[i].nb[(sneigh+1)%3],j,pts,cs[3*pts[i].nb[(sneigh+1)%3]+j]);
			for(j=0;j<3;j++)
				polyar(pts[i].nb[(sneigh+2)%3],j,pts,cs[3*pts[i].nb[(sneigh+2)%3]+j]);
			if(n<305000&&pts[i].x<0.21e9&&pts[i].y>3.3e9&&pts[i].x>2e8)
				printf("i: %d pt2: %d ran: %d sneigh: %d inpos: %d n: %d buf: %d\n",i,pt2,ran,sneigh,inpos,n,buf);//i: 177515 pt2: 58286 ran: 1 sneigh: 2 n: 299594
//Question: does this occur bf. or after swapping orders?
//		printf("Nachbarn von i: %d %d %d\n", nach[i][0], nach[i][1], nach[i][2]);
			return 1;
		}
		else return 0;
	}
	else return -2;
}

char tricheck(int i, struct pt *pts){ //auf return Condition umschreiben? Bleiben nach[i][0,1,2] im Cache?
	/*
	i has 3 neighbours, question: does one of them have a neighbour of i as neighbour? ^^ 
	Beware 'double triangles'
	*/
	char pos1,pos2;
	//Neighbour 0: (Are multiple if-conditions more economic?)
	if(pts[pts[i].nb[0]].nb[0]==pts[i].nb[1]||pts[pts[i].nb[0]].nb[0]==pts[i].nb[2]||pts[pts[i].nb[0]].nb[1]==pts[i].nb[1]||pts[pts[i].nb[0]].nb[1]==pts[i].nb[2]||pts[pts[i].nb[0]].nb[2]==pts[i].nb[1]||pts[pts[i].nb[0]].nb[2]==pts[i].nb[2]){
		pos1=0;
		if(pts[pts[i].nb[0]].nb[0]==pts[i].nb[1]||pts[pts[i].nb[0]].nb[1]==pts[i].nb[1]||pts[pts[i].nb[0]].nb[2]==pts[i].nb[1])
			return(2);
		else return(1);
	}
	else if(pts[pts[i].nb[1]].nb[0]==pts[i].nb[2]||pts[pts[i].nb[1]].nb[1]==pts[i].nb[2]||pts[pts[i].nb[1]].nb[2]==pts[i].nb[2]){
		pos1=1;pos2=2;
		return(0);
	}
	else return -1;
}

int t2(int i, struct pt *pts){
	int pt1,pt2;
	double triar;
	char noti, notpt1, notpt2;
	char ck2;
	if((noti=tricheck(i,pts))+1){//noti is in brackets ;)
		pt1=pts[i].nb[(noti+1)%3];
		pt2=pts[i].nb[(noti+2)%3];
//			printf("Triang found: %d %d %d\n",i,pt1,pt2);
		triar=fabs((double)(int)(pts[pt2].x-pts[i].x)*(int)(pts[pt1].y-pts[i].y)-(double)(int)(pts[pt1].x-pts[i].x)*(int)(pts[pt2].y-pts[i].y))*0.5;
//			printf("Fläche: %lf\n", triar);
		if(unlikely(triar<minar)){
//				printf("Kleines Dreieck! Fläche: %lf\n",triar);
			//Punkte: i, pt1, pt2. Zuerst schaun, wo welche Nachbarn gespeichert:
			notpt1=((pts[pt1].nb[0]==pt2&&pts[pt1].nb[1]==i)||(pts[pt1].nb[0]==i&&pts[pt1].nb[1]==pt2))?2:((pts[pt1].nb[0]==i||pts[pt1].nb[0]==pt2)?1:0);
			notpt2=((pts[pt2].nb[0]==pt1&&pts[pt2].nb[1]==i)||(pts[pt2].nb[0]==i&&pts[pt2].nb[1]==pt1))?2:((pts[pt2].nb[0]==i||pts[pt2].nb[0]==pt1)?1:0);
//				printf("Koordinaten der \"Außenpunkte\": %d %d %d\n",nach[i][noti],nach[pt1][notpt1],nach[pt2][notpt2]);
			if(pts[i].nb[noti]==pts[pt1].nb[notpt1]||pts[i].nb[noti]==pts[pt2].nb[notpt2]||pts[pt1].nb[notpt1]==pts[pt2].nb[notpt2]){
//				printf("Achtung Doppeldreieck, n: %d!\n",n);
//				write("doubpre",pts);
				doubtri(pts,i,pt1,pt2,noti,notpt1,notpt2);
//				write("doubnach",pts);
				return 1;
			}
			else{
				if(cs[3*i+(noti+2)%3]!=cs[3*pt1+(notpt1+2)%3]||cs[3*pt2+(notpt2+2)%3]!=cs[3*pt1+(notpt1+2)%3]||cs[3*i+(noti+2)%3]!=cs[3*pt2+(notpt2+2)%3]){
					printf("\n\n\nWarning!\n\n\n\ncs1: %d cs2: %d i: %d pt1: %d n: %d noti: %d notpt1: %d\n",cs[3*i+(noti+2)%3],cs[3*pt1+(notpt1+2)%3],i,pt1,n,noti,notpt1);
					if(cs[3*i+noti%3]==cs[3*pts[i].nb[noti]+(POSORNE(pts[i].nb[noti],i))%3]){
						if(cs[3*i+noti%3]==cs[3*i+(noti+1)%3]){
							printf("Passt: cs(noti): %d cs(noti+1): %d cs(noti+2): %d\n",cs[3*i+noti%3],cs[3*i+(noti+1)%3],cs[3*i+(noti+2)%3]);
							deltri(i,pt1,pt2,noti,notpt1,notpt2,pts,1);
							write("trinach",pts);
						}
						else exit(1);
					}
					else if(cs[3*pt1+notpt1%3]==cs[3*pts[pt1].nb[notpt1]+(POSORNE(pts[pt1].nb[notpt1],pt1))%3]){
						if(cs[3*pt1+notpt1%3]==cs[3*pt1+(notpt1+1)%3]){
							printf("Passt: cs(notpt1): %d cs(notpt1+1): %d cs(notpt1+2): %d\n",cs[3*pt1+notpt1%3],cs[3*pt1+(notpt1+1)%3],cs[3*pt1+(notpt1+2)%3]);
							deltri(i,pt1,pt2,noti,notpt1,notpt2,pts,2);
							write("trinach",pts);
						}
						else exit(1);
					}
					else if(cs[3*pt2+notpt2%3]==cs[3*pts[pt2].nb[notpt2]+(POSORNE(pts[pt2].nb[notpt2],pt2))%3]){
						if(cs[3*pt2+notpt2%3]==cs[3*pt2+(notpt2+1)%3]){
							printf("Passt: cs(notpt2): %d cs(notpt2+1): %d cs(notpt2+2): %d\n",cs[3*pt2+notpt2%3],cs[3*pt2+(notpt2+1)%3],cs[3*pt2+(notpt2+2)%3]);
							deltri(i,pt1,pt2,noti,notpt1,notpt2,pts,4);
							write("trinach",pts);
						}
						else exit(1);
					}
					else exit(1);
					return 0;
				}
				
				else deltri(i,pt1,pt2,noti,notpt1,notpt2,pts,0);
	
			}
			return 1;
		}
	}
	return 0;
		
}

int read(char *, char *, struct pt *);
int jump(int i,struct pt *points);
int write(char *suffix, struct pt *pts);
int main(int argc, char *argv[]){
	int test1,j,k,r;
	long long i;
	int dist1x,dist1y,dist2x,dist2y,dist3x,dist3y;
	char made[N/ABST+20]={0};
	struct pt *pts;
	char fstr[10];

	pts=(struct pt *) calloc(n,sizeof(struct pt));
	spolydata=(struct polydata *) calloc(sizeof(spolylist)/sizeof(int),sizeof(struct polydata));
	cs=(int *) calloc(3*n,sizeof(int));
	if(pts==NULL||cs==NULL){
		error(0,0,"Not enough memory!");
		exit(1);
	}
	if(read("mids","nach",pts)!=n){
		error(0,0,"Indizes stimmen ned ueberein!");
		exit(1);
	}
	j=mein_rand()%n;
	printf("j: %u\n",j);
	j=mein_rand()%n;
	cross(pts);
//	printf("Intersects? %d\ndist: %lf\n",intersects(pts,518606,518605),dist(pts,518605,518606));
	for(k=0;k<n;k++) sort(k,pts);
	i=0;
//	exit(1);
	init(pts);
	write("neu",pts);
	
	while(n>100000){ //Ursprünglich: 61400
//	while(i<3100000000){ //zerst 317
		j=pts[j].nb[mein_rand()%3];
		if(i%500==0)	j=mein_rand()%n;
		jump(j,pts);
//		if(sort(j,pts)==1&&j!=147980){ printf("Achtung bei %d, i: %lld\n",j,i); write("midstest","nachtest",pts); exit(1);}
		if(i%100000000==0){ printf("i: %lld n: %d\nsucc: %lld, trisucc: %lld, ratioall: %lf ratioot: %lf\n",i,n,succ,trisucc,(double)(expsucc)/trisucc,(double)(succ)/trisucc);}
		if(mdists[MEASP-n/MEASF]==0&&n<=(n-n%MEASF)==0){getdata(pts,spolydata);}
		
//		if((n==999000||n==999002)&&made[6]==0){ write("999000",pts); ++made[6]; printf("Writing 9.97\n");}
		if(n<=500000&&made[N/ABST+1]==0){ 
			for(r=0;r<RATEP;r++){
				ratedata[r][2]=mein_rand()%(N/2);
				ratedata[r][0]=ars[(int)ratedata[r][2]];
				ratedata[r][1]=ratedata[r][0];
				ratedata[r][3]=ck[(int)ratedata[r][2]];
			}
			printf("Meanar: %lf\n",meanar(pts));
			++made[N/ABST+1];}
		if(n<=480000&&made[N/ABST+2]==0){for(r=0;r<RATEP;r++){ratedata[r][1]=(ars[(int)(ratedata[r][2])]-ratedata[r][1]);} ++made[N/ABST+2];}

		if(n<=500000&&made[N/ABST+3]==0){
			for(r=0;r<N/2;r++){
				rate500000[r][0]=ars[r];
				rate500000[r][1]=ck[r];
			}
			++made[N/ABST+3];
			printf("Suc. jumps at %d: %lld\n",n,expsucc);
		}
		if(n<=490000&&made[N/ABST+4]==0){for(r=0;r<N/2;r++){rate500000[r][2]=(ars[r]-rate500000[r][0]);} ++made[N/ABST+4];printf("Suc. jumps at %d: %lld\n",n,expsucc);}
		if(n<=480000&&made[N/ABST+5]==0){for(r=0;r<N/2;r++){rate500000[r][3]=(ars[r]-rate500000[r][0]);} ++made[N/ABST+5];printf("Suc. jumps at %d: %lld\n",n,expsucc);}
		if(n<=470000&&made[N/ABST+6]==0){for(r=0;r<N/2;r++){rate500000[r][4]=(ars[r]-rate500000[r][0]);} ++made[N/ABST+6];printf("Suc. jumps at %d: %lld\n",n,expsucc);}
		if(n<=460000&&made[N/ABST+7]==0){for(r=0;r<N/2;r++){rate500000[r][5]=(ars[r]-rate500000[r][0]);} ++made[N/ABST+7];printf("Suc. jumps at %d: %lld\n",n,expsucc);}
		if(n<=450000&&made[N/ABST+8]==0){for(r=0;r<N/2;r++){rate500000[r][6]=(ars[r]-rate500000[r][0]);} ++made[N/ABST+8];printf("Suc. jumps at %d: %lld\n",n,expsucc);}
		if(n<=440000&&made[N/ABST+9]==0){for(r=0;r<N/2;r++){rate500000[r][7]=(ars[r]-rate500000[r][0]);} ++made[N/ABST+9];printf("Suc. jumps at %d: %lld\n",n,expsucc);}
		if(n<=430000&&made[N/ABST+10]==0){for(r=0;r<N/2;r++){rate500000[r][8]=(ars[r]-rate500000[r][0]);} ++made[N/ABST+10];printf("Suc. jumps at %d: %lld\n",n,expsucc);}
		if(n<=420000&&made[N/ABST+11]==0){for(r=0;r<N/2;r++){rate500000[r][9]=(ars[r]-rate500000[r][0]);} ++made[N/ABST+11];printf("Suc. jumps at %d: %lld\n",n,expsucc);}
		if(n<=410000&&made[N/ABST+12]==0){for(r=0;r<N/2;r++){rate500000[r][10]=(ars[r]-rate500000[r][0]);} ++made[N/ABST+12];printf("Suc. jumps at %d: %lld\n",n,expsucc);}
		if(n<=400000&&made[N/ABST+13]==0){for(r=0;r<N/2;r++){rate500000[r][11]=(ars[r]-rate500000[r][0]);} ++made[N/ABST+13];printf("Suc. jumps at %d: %lld\n",n,expsucc);}
		if(n<=495000&&made[N/ABST+14]==0){for(r=0;r<N/2;r++){rate500000[r][12]=(ars[r]-rate500000[r][0]);} ++made[N/ABST+14];printf("Suc. jumps at %d: %lld\n",n,expsucc);}
		if(n<=497500&&made[N/ABST+15]==0){for(r=0;r<N/2;r++){rate500000[r][13]=(ars[r]-rate500000[r][0]);} ++made[N/ABST+15];printf("Suc. jumps at %d: %lld\n",n,expsucc);}
		
		if(made[(n-n%ABST)/ABST]==0&&n<=n-n%ABST){
			sprintf(fstr,"%d",n-n%ABST);
			write(fstr,pts);
			++made[(n-n%ABST)/ABST];
			if(n%20000==0) sleep(3);
		}
//		if(test1==-1) break;
		i++;
	}
	write("neu",pts);
	printf("succ: %lld, trisucc: %lld, ratio: %lf\n",succ,trisucc,(double)(expsucc)/trisucc);
	printf("distdurch: %lf\n",distall/trisucc);
	printf("n: %d\n",n);
	return 0;
}

int read(char *filenamept, char *filenamenach, struct pt *pts){
	int i=0,j=0;
	FILE *fp;
	fp=fopen(filenamept,"r");
	while(!feof(fp)){
		fscanf(fp,"%u %u \n", &pts[i].x,&pts[i].y);
		i++;
	}
	
	fclose(fp);
	fp=fopen(filenamenach,"r");
	
	while(!feof(fp)){
		fscanf(fp,"%d %d %d \n",&pts[j].nb[0],&pts[j].nb[1],&pts[j].nb[2]);
		j++;
	}
	fclose(fp);
	if(i!=j){
		error(0,0,"Indizes stimmen ned ueberein!");
		exit(1);
		}
	return i;
}

int write(char *suffix, struct pt *pts){
	FILE *fp,*fp2;
	char *name;
	name=(char *) calloc(50,sizeof(char));
	int i,j;
	strcpy(name,BASEDIR "mids");
	strcat(name,suffix);
	fp=fopen(name,"w");
	for(i=0;i<n;i++){
		fprintf(fp,"%u %u\n", pts[i].x,pts[i].y);
	}
	fclose(fp);
	
	strcpy(name,BASEDIR "nach");
	strcat(name,suffix);
	fp=fopen(name,"w");
	
	for(i=0;i<n;i++){
		fprintf(fp,"%d %d %d\n",pts[i].nb[0],pts[i].nb[1],pts[i].nb[2]);
	}
	fclose(fp);
	
	
	strcpy(name,BASEDIR "cons");
	strcat(name,suffix);
	fp=fopen(name,"w");
	
	for(i=0;i<n;i++){
		fprintf(fp,"%d %d %d\n",cs[3*i],cs[3*i+1],cs[3*i+2]);
	}
	fclose(fp);
	
	strcpy(name,BASEDIR "polys");
	strcat(name,suffix);
	fp=fopen(name,"w");
	
	i=0;
	for(i=0;i<N;i++){
		fprintf(fp,"%d %lf\n",ck[i],ars[i]);
	}
	fclose(fp);
	fp=fopen(BASEDIR "mdists","w");
	
	i=0;
	for(i=0;i<MEASP;i++){
		if(mdists[i]==0){fprintf(fp,"NaN %d\n",N-N*i/MEASP);}
		else fprintf(fp,"%lf %lf %lf %d\n",mdists[i],msides[i],timemeas[i],N-N*i/MEASP);
	}
	fclose(fp);
	
	for(i=0;i<sizeof(spolylist)/sizeof(int);i++){
		sprintf(name,BASEDIR "spoly%d",spolylist[i]);
		fp=fopen(name,"w");
		for(j=0;j<MEASP;j++){
			fprintf(fp,"%lf %d %lf\n",spolydata[i].ar[j],spolydata[i].sides[j],timemeas[j]);
		}
		fclose(fp);
	}

	fp=fopen(BASEDIR "rate","w");
	for(i=0;i<RATEP;i++){
		fprintf(fp,"%lf %lf %lf %d\n",ratedata[i][0],ratedata[i][1]/ratedata[i][0],ratedata[i][1],(int)ratedata[i][3]);	//A(t1),(A(t2)-A(t1))/A(t1),Absdiff.,sidesnum
	}
	fclose(fp);

	fp=fopen(BASEDIR "hist","w");
	for(i=0;i<MEASP;i++){
		for(j=0;j<30;j++){
			fprintf(fp,"%d ",hist[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	
	fp=fopen(BASEDIR "meansize","w");
	for(i=0;i<MEASP;i++){
		for(j=0;j<30;j++){
			fprintf(fp,"%lf ",meansize[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	
	fp=fopen(BASEDIR "rate500000","w");
	for(i=0;i<N/2;i++){
		fprintf(fp,"%lf ",rate500000[i][0]);
		for(j=2;j<14;j++){
			fprintf(fp,"%lf %lf ",rate500000[i][j]/rate500000[i][0],rate500000[i][j]);
		}
		fprintf(fp,"%d\n",(int) rate500000[i][1]);
	}
	fclose(fp);
	return 0;
}

int jump(int i,struct pt *pts){ /*Choose random vector (dx,dy)*/
//	printf("%d\n",i);
	int dist1x, dist1y, dist2x,dist2y, dist3x, dist3y,dx,dy,count=0;
	double distxgq, distygq, scal;
	double maxlen;

	dist1x=pts[i].x-pts[pts[i].nb[0]].x;
	dist1y=pts[i].y-pts[pts[i].nb[0]].y;
	dist2x=pts[i].x-pts[pts[i].nb[1]].x;
	dist2y=pts[i].y-pts[pts[i].nb[1]].y;
	dist3x=pts[i].x-pts[pts[i].nb[2]].x;
	dist3y=pts[i].y-pts[pts[i].nb[2]].y;

	distxgq=pow(dist1x,2)+pow(dist2x,2)+pow(dist3x,2);
	distygq=pow(dist1y,2)+pow(dist2y,2)+pow(dist3y,2);
	scal=(double)dist1x*dist1y+(double)dist2x*dist2y+(double)dist3x*dist3y;
	if(scal*scal>=distxgq*distygq)	return 0;//prevent division by sqrt of number <=0
	maxlen=2.83*jumpar/sqrt(distxgq+distygq-sqrt(pow(distxgq-distygq,2)+4*pow(scal,2)));//This was a nice math exercise :)
//	printf("distxgq: %lf distygq: %lf scal: %lf\n dist1x: %d dist1y: %d dist2x: %d dist2y: %d dist3x: %d dist3y: %d\nmaxlen: %lf",distxgq, distygq,scal,dist1x,dist1y,dist2x,dist2y,dist3x,dist3y,maxlen);
//	printf("%lf\n", maxlen);
	if(maxlen>2147483647)	return 0; //too big
	dx=((int) mein_rand())%(int)maxlen;
	dy=((int) mein_rand())%(int)maxlen;
	if(jumptest(dx,dy,(int)dist1x, (int)dist1y, (int)dist2x,(int)dist2y, (int)dist3x,(int) dist3y)){
//	printf("Geht bei %d, gibt aus %d\n",i,j);
//	printf("dx: %d dy: %d\n",dx,dy);
//	if(j==2)
//		printf("Index: %d\n",i);
	pts[i].x+=dx;
	pts[i].y+=dy;
	t+=(double)sqrt(N)/pow(n,1.5);
	while(t1(i,pow(dist1x,2)+pow(dist1y,2),pow(dist2x,2)+pow(dist2y,2),pow(dist3x,2)+pow(dist3y,2),pts)==1&&count<T1TRIES){
		dist1x=pts[i].x-pts[pts[i].nb[0]].x;
		dist1y=pts[i].y-pts[pts[i].nb[0]].y;
		dist2x=pts[i].x-pts[pts[i].nb[1]].x;
		dist2y=pts[i].y-pts[pts[i].nb[1]].y;
		dist3x=pts[i].x-pts[pts[i].nb[2]].x;
		dist3y=pts[i].y-pts[pts[i].nb[2]].y;
//		t1(i,pow(dist1x,2)+pow(dist1y,2),pow(dist2x,2)+pow(dist2y,2),pow(dist3x,2)+pow(dist3y,2),pts);
		count++;
//		if(count>3) printf("Count: %d\n",count);
	}
	t2(i,pts);
	
	}
//	else printf("Geht bei %d nicht \n",i);
//	return i;
	}

int getdata(struct pt *pts, struct polydata *spolydata){
	int i;
	meanar(pts);
	mdists[MEASP-n/MEASF]=meandist(pts);
	msides[MEASP-n/MEASF]=meansides(pts);
	timemeas[MEASP-n/MEASF]=t;
	for(i=0;i<sizeof(spolylist)/sizeof(int);i++){
		spolydata[i].ar[MEASP-n/MEASF]=ars[spolylist[i]];
		spolydata[i].sides[MEASP-n/MEASF]=ck[spolylist[i]];
	}
	gethist(pts,hist[MEASP-n/MEASF],meansize[MEASP-n/MEASF]);
	
}
