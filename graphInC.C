/*  Graph Algorithms in C*/
/* 
 * Author PAURUSH PRAVEEN
 * Algoritms to read and process graph structures and 
 * graph theory applications
 * Objective: Perform consistency mappings of PPI databases
*/


// 1. Read an adjacency matrix file
#include<stdio.h>
#define Vmax 10
int j, x, y, V, E;
int g[Vmax][Vmax];
adjmatrix()
  {
    scanf("%d %d\n" &V, &E);
    for(x=1; x<=V; x++)
	for(y=1; y<=V; y++)
	  g[x][y]=0;
	for(x=1;x<=V;x++)
	  g[x][y]=1;
	for(j=1;j<=E;j++)
	{
	  scanf("%c %c\n", &v1, &v2);
	  x=index(v1);
	  y=index(v2);
	  g[x][y]=1;
	  g[y][x]=1;
	}
  }
  
	  

/*
#include<stdio.h>
//#include<conio.h>
struct n{
  int data;
  struct n *link;
};
typedef struct n NODE;
NODE *getnode(int);
NODE *findlast(NODE *);
void display(NODE *[],int);
void main(){
  NODE *ptr,*temp,*h[10];
  int n,a[10][10],i,j;
  clrscr();
  printf("\n Enter total number of vertices : ");
  scanf("%d",&n);
  printf("\n Enter entries of an adjacency matrix :\n");
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      printf("\n Enter a[%d][%d] : ",i+1,j+1);
      scanf("%d",&a[i][j]);
    }
  }
  printf("\n Entered Adjacency matrix is ... \n");
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      printf("\t %d ",a[i][j]);
    }
    printf("\n");
  }
  for(i=0;i<n;i++)
    h[i]=NULL;
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
	if(a[i][j]==1){
	  temp=getnode(j+1);
	  if(h[i]==NULL)
	    h[i]=temp;
	  else{
	    ptr=findlast(h[i]);
	    ptr->link=temp;
	  }
	}
      }
    }
    printf("\n The Adjacency list looks like ... \n");
    display(h,n);
    getch();
    }
    NODE *getnode(int j){
      NODE * temp;
      temp=(NODE *)malloc(sizeof(NODE));
      temp->data=j;
      temp->link=NULL;
      return(temp);
    }
    NODE *findlast(NODE *h){
      NODE *ptr;
      for(ptr=h;ptr->link!=NULL;ptr=ptr->link);
	return(ptr);
    }
    void display(NODE *h[10],int n){
      NODE *ptr;
      int i;
      for(i=0;i<n;i++){
	printf("\n V%d ",i+1);
	ptr=h[i];
	if(ptr==NULL)
	  printf(" NULL");
	while(ptr!=NULL){
	  printf(" ->V%d",ptr->data);
	  ptr=ptr->link;
	}
	printf("\n");
    }
  }
*/
// 2. Read an sif file file

// 3. Read an adjacency list file


// 4. Read an XGML file 

// 5. Dijkastra algorithm to find shortest path between two vertices

// 6. Steiner tree implementation based on Sadeghi et al.

// 7. 