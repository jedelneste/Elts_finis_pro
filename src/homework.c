# include "tsunami.h"


double interpolate(double *phi, double *U, int *map, int n)
{
    double u = 0.0; int i;
    for (i=0; i <n; i++)
        u += phi[i]*U[map[i]];
    return u;
}


void stereoCoordonnee(double* X, double* Y, double* x, double* y, double* z, int nNode){

    for(int i=0; i<nNode; i++){
        X[i] = (2*R*x[i])/(R+z[i]);
        Y[i] = (2*R*y[i])/(R+z[i]);
    }

}


void tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName)
{ 

//
// A modifer :-)
// Ici, on se contente d'initialiser le champs d'�l�vation avec la bathym�trie
// pour avoir un joli plot de d�part !
//
// On divise la bathym�trie par un facteur (10*BathMax) pour que l'�chelle de couleur fonctionne
// bien dans la fonction d'animation fournie....
//

    int i,j,nNode,nElem, nEdges, trash;
    double dtrash;
    
    double BathMax = 9368;
    
    FILE* file = fopen(meshFileName,"r");
    fscanf(file, "Number of nodes %d \n",&nNode);   
    double *bath = malloc(sizeof(double)*nNode);
    double *x = malloc(sizeof(double)*nNode);
    double *y = malloc(sizeof(double)*nNode);
    for (i = 0; i < nNode; i++) 
        fscanf(file,"%d : %le %le %le\n",&trash,&x,&y,&bath[i]); 
    fscanf(file, "Number of triangles %d \n",&nElem); 
    int *elem = malloc(sizeof(int)*3*nElem);
    for (i = 0; i < nElem; i++) 
        fscanf(file,"%d : %d %d %d \n",&trash,&elem[i*3],&elem[i*3+1],&elem[i*3+2]);   
    fscanf(file, "Number of edges %d \n", &nEdges);
    femEdge *edges = malloc(sizeof(femEdge)*nEdges);
    for (i = 0; i<nEdges; i++){
        fscanf(file, "%d : %d %d : %d %d \n", &trash, &(edges[i].node[0]), &(edges[i].node[1]), &(edges[i].elem[0]), &(edges[i].elem[1]));
        fprintf(stdout, "%d %d \n", &(edges[i].node[0]), &(edges[i].node[1]));
    }
    fclose(file); 
    
    double *E  = malloc(sizeof(double)*nElem*3);
    for (i = 0; i < nElem; i++)
        for (j = 0; j < 3; j++)
            E[i*3+j] = bath[elem[i*3+j]]/(10*BathMax);
  
    tsunamiWriteFile(baseResultName,0,E,E,E,nElem,3); 
    
    free(bath);
    free(E);
    free(elem);
 
}








