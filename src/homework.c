# include "tsunami.h"

typedef struct{
    int elem[2];
    int node[2];
} femEdge; 

double interpolate(double *phi, double *U, int *map, int n){

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

void tsunamiTriangleMap(int index, int map[3])
{
    int j;
    for (j=0; j < 3; ++j) 
        map[j] = index*3 + j; 
}

void tsunamiEdgeMap(int index, int map[2][2], femEdge *edges, int nElem, int* elem)
{
    int i,j,k;   
    for (j=0; j < 2; ++j) {
        int node = edges[index].node[j];
        for (k=0; k < 2; k++) {
            int e = edges[index].elem[k];
            map[k][j] = (nElem)*3;
            if (e >= 0) {
                for (i=0; i < 3; i++) {
                    if (elem[e*3 + i] == node) {
                        map[k][j] = e*3 + i;  }}}}}
}

void phiCreate1(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;  
    phi[1] = xsi;
    phi[2] = eta;
}

void phiCreate2(double xsi, double *phi) 
{   
    phi[0] = (1.0 - xsi)/2.0;
    phi[1] = (1.0 + xsi)/2.0;   
    
}

int femEdgesCompare(const void *edgeOne, const void *edgeTwo)
{
    int *nodeOne = ((femEdge*) edgeOne)->node;
    int *nodeTwo = ((femEdge*) edgeTwo)->node;  
    int  diffMin = fmin(nodeOne[0],nodeOne[1]) - fmin(nodeTwo[0],nodeTwo[1]);
    int  diffMax = fmax(nodeOne[0],nodeOne[1]) - fmax(nodeTwo[0],nodeTwo[1]);
    
    if (diffMin < 0)    return  1;
    if (diffMin > 0)    return -1;
    if (diffMax < 0)    return  1;
    if (diffMax > 0)    return -1; 
                        return  0;
}

void InverseMatrix(double* BE, double* BU, double* BV, int* elem, int nElem, double* X, double* Y)
{
         
    double invA[3][3] = {{18.0,-6.0,-6.0},{-6.0,18.0,-6.0},{-6.0,-6.0,18.0}};
    double BEloc[3],BUloc[3],BVloc[3];
 
    double xLoc[3],yLoc[3],jac;
    int    ielem,i,j,mapElem[3];
    
    for (ielem=0; ielem < nElem; ielem++) {
        tsunamiTriangleMap(ielem,mapElem);
        int *mapCoord = &(elem[ielem*3]);
        for (j=0; j < 3; ++j) {
        	  xLoc[j] = X[mapCoord[j]];
        	  yLoc[j] = Y[mapCoord[j]]; }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        for (i=0; i < 3; i++) {
            BEloc[i] = BE[mapElem[i]];
            BUloc[i] = BU[mapElem[i]];
            BVloc[i] = BV[mapElem[i]];
            BE[mapElem[i]] = 0.0; 
            BU[mapElem[i]] = 0.0; 
            BV[mapElem[i]] = 0.0; }
        for (i=0; i < 3; i++) { 
            for (j=0; j < 3; j++) {
                BE[mapElem[i]] += invA[i][j] * BEloc[j] / jac; 
                BU[mapElem[i]] += invA[i][j] * BUloc[j] / jac; 
                BV[mapElem[i]] += invA[i][j] * BVloc[j] / jac; }}}

}

void tsunamiAddIntegralsElements(int nElem, int* elem, double* X, double* Y, double* E, double* U, double* V, double* BE, double* BU, double* BV){

    int iElem, i, j, k, mapElem[3];
    double weight,jac,xLoc[3],yLoc[3],phi[3],dphidx[3],dphidy[3];
    double x,y,u,v,e;

    for (iElem=0; iElem<nElem; iElem++){
        tsunamiTriangleMap(iElem, mapElem);
        int *mapCoord = &(elem[iElem*3]);
        for(j=0; j<3; j++){
            xLoc[j] = X[mapCoord[j]];
            yLoc[j] = Y[mapCoord[j]];
        }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        for(i=0; i<3; i++){
            dphidx[i] = (yLoc[(i+1)%3] - yLoc[(i+2)%3])/jac;
            dphidy[i] = (xLoc[(i+1)%3] - xLoc[(i+2)%3])/jac;
        }
        for (k=0; k<3; k++){
            weight = gaussTriangleWeight[k];
            phiCreate1(gaussTriangleXsi[k], gaussTriangleEta[k], phi);
            x = interpolate(phi, X, mapCoord, 3);
            y = interpolate(phi, Y, mapCoord, 3);
            e = interpolate(phi, E, mapElem, 3);
            u = interpolate(phi, U, mapElem, 3); 
            v = interpolate(phi, V, mapElem, 3);
            double z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
            double theta = asin(z3d/R)*180/PI;
            double coriolis = 2*Omega*sin(theta); 
            double h = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
            double toAdd = (4*R*R+x*x+y*y)/(4*R*R);
            for(i=0; i<3; i++){
                BE[mapElem[i]] += ( h*u*dphidx[i] + h*v*dphidy[i] )*toAdd*jac*weight + phi[i]*((h*(x*u+v*y))/R*R)*jac*weight; 
                BU[mapElem[i]] += ((coriolis*v - Gamma*u)*phi[i] + g*e*dphidx[i]*toAdd)*jac*weight + (phi[i]*((g*x*e)/2*R*R))*jac*weight; 
                BV[mapElem[i]] += ((- coriolis*u - Gamma*v)*phi[i] + g*e*dphidy[i]*toAdd)*jac*weight + (phi[i]*((g*y*e)/2*R*R))*jac*weight;}}}         
}

void tsunamiAddIntegralsEdges(int nEdges, femEdge* edges, int nElem, int* elem, double* X, double* Y, double* E, double* U, double* V, double* BE, double* BU, double* BV){

    int iEdge, mapEdge[2][2], j, k;
    double xEdge[2], yEdge[2], phiEdge[2];
    double xsi,weight,jac, eL, eR, uL, uR, vL, vR, unL, unR, qe, qu, qv;

    for(iEdge=0; iEdge<nEdges; iEdge++){
        tsunamiEdgeMap(iEdge, mapEdge, edges, nElem, elem);
        for (j=0; j<2; j++){
            int node = edges[iEdge].node[j];
            xEdge[j] = X[node];
            yEdge[j] = Y[node];
        }

        int boundary = (mapEdge[1][0] == nElem*3);

        double dxdxsi = (xEdge[1] - xEdge[0]);
        double dydxsi = (yEdge[1] - yEdge[0]);
        double norm = sqrt(dxdxsi*dxdxsi + dydxsi*dydxsi);
        double nx =  dydxsi/norm;
        double ny = -dxdxsi/norm;
        jac = norm / 2.0;
        for(k=0; k<2; k++){
            xsi = gaussEdgeXsi[k];
            weight = gaussEdgeWeight[k];
            phiCreate2(xsi, phiEdge);
                        
            //A vérifier
            double h = R*(4*R*R - xEdge[k]*xEdge[k] - yEdge[k]*yEdge[k]) / (4*R*R + xEdge[k]*xEdge[k] + yEdge[k]*yEdge[k]);
            eL = interpolate(phiEdge,E,mapEdge[0],2);
            eR = boundary ? eL : interpolate(phiEdge,E,mapEdge[1],2);
            uL = interpolate(phiEdge,U,mapEdge[0],2);
            uR = interpolate(phiEdge,U,mapEdge[1],2);
            vL = interpolate(phiEdge,V,mapEdge[0],2);
            vR = interpolate(phiEdge,V,mapEdge[1],2);
            unL = uL*nx+ vL*ny;
            unR = boundary ? -unL : uR*nx + vR*ny;
            qe =  0.5*h*   ( (unL+unR) + sqrt(g/h)*( eL-eR ) );
            qu =  0.5*g*nx*( ( eL+eR ) + sqrt(h/g)*(unL-unR) );
            qv =  0.5*g*ny*( ( eL+eR ) + sqrt(h/g)*(unL-unR) );    
            double toAdd = (4*R*R+xEdge[k]*xEdge[k]+yEdge[k]*yEdge[k])/(4*R*R);    
            for (int i=0; i < 2; i++) {
                BE[mapEdge[0][i]] -= qe*phiEdge[i]*toAdd*jac*weight; 
                BU[mapEdge[0][i]] -= qu*phiEdge[i]*toAdd*jac*weight; 
                BV[mapEdge[0][i]] -= qv*phiEdge[i]*toAdd*jac*weight; 
                BE[mapEdge[1][i]] += qe*phiEdge[i]*toAdd*jac*weight;
                BU[mapEdge[1][i]] += qu*phiEdge[i]*toAdd*jac*weight;
                BV[mapEdge[1][i]] += qv*phiEdge[i]*toAdd*jac*weight; }}}
}

void initialCondition(double *U, double *V, int nElem, int *elem, double *X, double *Y){
    int e,j,*node;
    for(int e=0; e<nElem; e++){
        node = &(elem[e*3]);
        for(int j=0; j<3; j++){
            double x = X[node[j]];
            double y = Y[node[j]];
            U[e*3+j] = tsunamiInitialConditionOkada(x,y);
            V[e*3+j] = tsunamiInitialConditionOkada(x,y);
        }
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

    int i,j,nNode,nElem, nEdges, trash, iElem, mapElem[3];
    double dtrash, xLoc[3], yLoc[3], phi[3], dphidx[3], dphidy[3];
    
    double BathMax = 9368;
    
    FILE* file = fopen(meshFileName,"r");
    fscanf(file, "Number of nodes %d \n",&nNode);   
    double *bath = malloc(sizeof(double)*nNode);
    double *x = malloc(sizeof(double)*nNode);
    double *y = malloc(sizeof(double)*nNode);
    for (i = 0; i < nNode; i++) 
        fscanf(file,"%d : %le %le %le\n",&trash,&x[i],&y[i],&bath[i]); 
    fscanf(file, "Number of triangles %d \n",&nElem); 
    int *elem = malloc(sizeof(int)*3*nElem);
    for (i = 0; i < nElem; i++) 
        fscanf(file,"%d : %d %d %d \n",&trash,&elem[i*3],&elem[i*3+1],&elem[i*3+2]);   
    fscanf(file, "Number of edges %d \n", &nEdges);
    femEdge *edges = malloc(sizeof(femEdge)*nEdges);
    for (i = 0; i<nEdges; i++){
        fscanf(file, "%d : %d %d : %d %d \n", &trash, &(edges[i].node[0]), &(edges[i].node[1]), &(edges[i].elem[0]), &(edges[i].elem[1]));
        //fprintf(stdout, "%d %d \n", (edges[i].node[0]), (edges[i].node[1]));
    }
    qsort(edges, nEdges, sizeof(femEdge), femEdgesCompare);
    int index = 0;          
    int nBoundary = 0;
    
    for (i=0; i < nEdges; i++) {
      if (i == nEdges - 1 || femEdgesCompare(&edges[i],&edges[i+1]) != 0) {
              edges[index] = edges[i];
              nBoundary++; }
      else {  edges[index] = edges[i];
              edges[index].elem[1] = edges[i+1].elem[0];
              i = i+1;}
      index++; 
    }
    double *X = malloc(sizeof(double)*nNode);
    double *Y = malloc(sizeof(double)*nNode);
    stereoCoordonnee(X,Y,x,y,bath,nNode);

    fclose(file); 

    //Création des éléments qu'on a besoin 
    double *E  = malloc(sizeof(double)*nElem*3); //peut être mettre +1
    double *U = malloc(sizeof(double)*nElem*3);
    double *V = malloc(sizeof(double)*nElem*3);
    double *BE  = malloc(sizeof(double)*nElem*3);
    double *BU = malloc(sizeof(double)*nElem*3);
    double *BV = malloc(sizeof(double)*nElem*3);
    for (i =0; i< nElem*3; i++){
        E[i] = bath[elem[i]]/(10*BathMax);
        U[i] = 0.0;
        V[i] = 0.0;
        BE[i] = 0.0;
        BU[i] = 0.0;
        BV[i] = 0.0;
    }

    // On applique les conditions initiales
    initialCondition(U,V,nElem, elem,X,Y);
    tsunamiAddIntegralsElements(nElem, elem, X,Y,E,U,V,BE,BU,BV);
    tsunamiAddIntegralsEdges(nEdges,edges,nElem,elem,X,Y,E,U,V,BE,BU,BV);
    InverseMatrix(BE,BU,BV,elem,nElem,X,Y);
    for (i=0; i < nElem*3; i++) {
        E[i] += dt * BE[i];
        U[i] += dt* BU[i];
        V[i] += dt * BV[i]; 
    }
   
   //Calcul avec le nombre d'itération 
   /*
   int iteration = 0;

   while(iteration < nmax){
        for (i=0; i < nElem*3; i++) {
            BE[i] = 0.0;
            BU[i] = 0.0;
            BV[i] = 0.0; 
        }
        tsunamiAddIntegralsElements(nElem, elem, X,Y,E,U,V,BE,BU,BV);
        tsunamiAddIntegralsEdges(nEdges,edges,nElem,elem,X,Y,E,U,V,BE,BU,BV);
        InverseMatrix(BE,BU,BV,elem,nElem,X,Y);
        //Mise à jour de E, U, V
        for (i=0; i < nElem*3; i++) {
            E[i] += dt * BE[i];
            U[i] += dt* BU[i];
            V[i] += dt * BV[i]; }
        //On regarde s'il faut enregistrer dans un fichier ou pas
        if( (iteration+1) % sub == 0){
            tsunamiWriteFile(baseResultName,iteration,U,V,E,nElem,3); 
        }
        iteration ++;    
   }
   */
   //On libère la mémoire utilisée
    free(bath);
    free(E);
    free(elem);
    free(U);
    free(V);
    free(edges);
    free(BE);
    free(BU);
    free(BV);
    free(X);
    free(Y);
    free(x);
    free(y);
}








